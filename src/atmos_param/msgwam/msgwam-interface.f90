module msgwam_mod

! ==============================================================================
! This module implements MS-GWaM as described in Bölöni et al. (2021) and as
! implemented in dsconnelly/python-msgwam on GitHub, though there are some
! differences due to the Fortran translation.
! ==============================================================================

use constants_mod,        only: constants_init
use fms_mod,              only: CLOCK_ROUTINE, fms_init, mpp_clock_begin, &
                                mpp_clock_end, mpp_clock_id, MPP_CLOCK_SYNC, &
                                mpp_pe, write_version_number
use time_manager_mod,     only: time_manager_init, time_type

use msgwam_constants_mod, only: i_max, init_msgwam_constants, j_max, n_max, &
                                n_source, print_prune_diag, q_max
use msgwam_io_mod,        only: init_nc_output, init_ray_state, &
                                save_ray_state, send_nc_output
use msgwam_mean_mod,      only: get_accelerations, project_fluxes, &
                                update_mean_fields
use msgwam_rays_mod,      only: t_ray
use msgwam_RK4_mod,       only: take_RK4_step
use msgwam_sinks_mod,     only: apply_breaking, apply_dissipation, &
                                check_boundaries
use msgwam_source_mod,    only: LONG_KIND, check_source, init_source, &
                                t_prune_diag

implicit none
private

character(len=128) :: version = "msgwam.f90, 2025/07/08"
character(len=128) :: tagname = "cayuga"

! ==============================================================================
! public interfaces for use by damping_driver
! ==============================================================================

public msgwam_calc, msgwam_end, msgwam_init

! ==============================================================================
! module status and timing variables
! ==============================================================================

logical :: is_first_step = .true.
logical :: is_initialized = .false.
integer :: clock_rt, clock_io

! ==============================================================================
! mean state variables
! ==============================================================================

real, dimension(:, :, :), allocatable :: u, v, rho, N2, G2, flux_x, flux_y

! ==============================================================================
! ray volume state variables
! ==============================================================================

type(t_ray), dimension(:, :, :), allocatable :: rays
integer, dimension(:, :, :), allocatable :: ghosts
integer, dimension(:, :), allocatable :: last_meta

contains

subroutine msgwam_init(lon_bounds, lat_bounds, p_ref, Time, axes)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:),    intent(in) :: lon_bounds, lat_bounds, p_ref
    type(time_type),       intent(in) :: Time
    integer, dimension(4), intent(in) :: axes

    ! --------------------------------------------------------------------------

    if (is_initialized) then
        return
    end if

    call fms_init
    call time_manager_init
    call constants_init

    call init_msgwam_constants(lon_bounds, lat_bounds, p_ref)

    allocate(u(q_max, i_max, j_max))
    allocate(v(q_max, i_max, j_max))

    allocate(rho(q_max, i_max, j_max))
    allocate(N2(q_max, i_max, j_max))
    allocate(G2(q_max, i_max, j_max))

    allocate(flux_x(q_max + 1, i_max, j_max))
    allocate(flux_y(q_max + 1, i_max, j_max))

    allocate(rays(n_max, i_max, j_max))
    allocate(ghosts(n_source, i_max, j_max))
    allocate(last_meta(i_max, j_max))

    call init_source(lat_bounds)
    call init_nc_output(axes, Time)
    call init_ray_state(rays, ghosts, last_meta)

    clock_rt = mpp_clock_id( &
        "      MS-GWaM", &
        grain=CLOCK_ROUTINE, &
        flags=MPP_CLOCK_SYNC &
    )

    clock_io = mpp_clock_id( &
        "      MS-GWaM IO", &
        grain=CLOCK_ROUTINE, &
        flags=MPP_CLOCK_SYNC &
    )

    is_initialized = .true.

end subroutine msgwam_init

subroutine msgwam_calc(i_start, j_start, lat, p_full, z_full, temp, uuu, vvv, &
    Time, dt, du_dt, dv_dt)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                              intent(in)  :: i_start, j_start
    real, dimension(:, :),                intent(in)  :: lat
    real, dimension(i_max, j_max, q_max), intent(in)  :: p_full, z_full, temp, &
                                                         uuu, vvv
    type(time_type),                      intent(in)  :: Time
    real,                                 intent(in)  :: dt
    real, dimension(i_max, j_max, q_max), intent(out) :: du_dt, dv_dt

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------    
    integer :: i, j
    real :: dt_rays
    type(t_prune_diag), dimension(i_max, j_max) :: diag

    ! --------------------------------------------------------------------------

    if (is_first_step) then
        is_first_step = .false.
        dt_rays = dt
    else
        dt_rays = dt / 2.
    end if

    call mpp_clock_begin(clock_rt)
    !$OMP PARALLEL DO COLLAPSE(2)
    do j = 1, j_max
        do i = 1, i_max

            call advance_column( &
                j, &
                dt_rays, &
                z_full(i, j, :), &
                p_full(i, j, :), &
                temp(i, j, :), &
                uuu(i, j, :), &
                vvv(i, j, :), &
                u(:, i, j), &
                v(:, i, j), &
                rho(:, i, j), &
                N2(:, i, j), &
                G2(:, i, j), &
                rays(:, i, j), &
                ghosts(:, i, j), &
                last_meta(i, j), &
                diag(i, j), &
                flux_x(:, i, j), &
                flux_y(:, i, j), &
                du_dt(i, j, :), &
                dv_dt(i, j, :) &
            )

        end do
    end do
    !$OMP END PARALLEL DO
    call mpp_clock_end(clock_rt)

    call mpp_clock_begin(clock_io)
    call send_nc_output(i_start, j_start, Time, u, v, rho, N2, G2, &
        flux_x, flux_y, du_dt, dv_dt)
    call mpp_clock_end(clock_io)

    if (print_prune_diag) then
        call print_pruning_diagnostics(diag)
    end if

end subroutine msgwam_calc

subroutine advance_column(j, dt, z_full, p_full, temp, uuu, vvv, &
    u, v, rho, N2, G2, rays, ghosts, last_meta, diag, &
    flux_x, flux_y, du_dt, dv_dt)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                       intent(in)    :: j
    real,                          intent(in)    :: dt
    real, dimension(q_max),        intent(in)    :: z_full, p_full, temp, &
                                                    uuu, vvv
    real, dimension(q_max),        intent(out)   :: u, v, rho, N2, G2
    type(t_ray), dimension(n_max), intent(inout) :: rays
    integer, dimension(n_source),  intent(inout) :: ghosts
    integer,                       intent(inout) :: last_meta
    type(t_prune_diag),            intent(out)   :: diag
    real, dimension(q_max + 1),    intent(out)   :: flux_x, flux_y
    real, dimension(q_max),        intent(out)   :: du_dt, dv_dt

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real, dimension(0:q_max + 1) :: z_centers
    real, dimension(q_max + 1)   :: z_faces

    ! --------------------------------------------------------------------------

    call update_mean_fields(z_full, p_full, temp, uuu, vvv, &
        z_centers, z_faces, u, v, rho, N2, G2)

    call take_RK4_step(j, z_centers(1:q_max), u, v, N2, G2, dt, rays)

    call apply_dissipation(j, z_faces, z_centers(1:q_max), rho, dt, rays)
    call apply_breaking(j, z_faces, rho, N2, rays)
    call check_boundaries(j, z_centers, rays)

    call check_source(j, z_centers(1:q_max), u, v, N2, G2, dt, &
        rays, ghosts, last_meta, diag)

    call project_fluxes(z_centers, rays, flux_x, flux_y)
    call get_accelerations(z_faces, rho, flux_x, flux_y, du_dt, dv_dt)

end subroutine advance_column

subroutine print_pruning_diagnostics(diag)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    type(t_prune_diag), dimension(i_max, j_max), intent(in) :: diag

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: n_pruned
    real :: total_flux, total_r
    integer(kind=LONG_KIND) :: total_age

    ! --------------------------------------------------------------------------

    n_pruned = sum(diag(:, :)%n_pruned)
    total_age = sum(diag(:, :)%total_age)
    total_flux = sum(diag(:, :)%total_flux)
    total_r = sum(diag(:, :)%total_r)

    write(*, "(I2, A, I6, A, I12, A, A, F12.6, A, F14.4, A)") &
        mpp_pe(), " pruned ", n_pruned, &
        " rays with total age ", total_age, " seconds", &
        " and total flux fraction ", total_flux, &
        " and total r ", total_r, " kilometers"

end subroutine print_pruning_diagnostics

subroutine msgwam_end

    call save_ray_state(rays, ghosts, last_meta)
    is_initialized = .false.

end subroutine msgwam_end

end module msgwam_mod