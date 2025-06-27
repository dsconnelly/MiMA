module msgwam_mod

! ==============================================================================
! This module implements MS-GWaM as described in Bölöni et al. (2021) and as
! implemented in dsconnelly/python-msgwam on GitHub, though there are some
! differences due to the Fortran translation.
! ==============================================================================

use constants_mod,        only: constants_init
use fms_mod,              only: CLOCK_ROUTINE, fms_init, mpp_clock_begin, &
                                mpp_clock_end, mpp_clock_id, MPP_CLOCK_SYNC, &
                                write_version_number
use time_manager_mod,     only: time_manager_init, time_type

use msgwam_constants_mod, only: i_max, init_msgwam_constants, j_max, n_max, &
                                n_source, q_max, steady_state
use msgwam_debug_mod,     only: check_rays, track_ray
use msgwam_io_mod,        only: init_nc_output, init_ray_state, &
                                save_ray_state, send_nc_output
use msgwam_mean_mod,      only: get_accelerations, project_fluxes, &
                                update_mean_fields
use msgwam_rays_mod,      only: t_ray
use msgwam_RK4_mod,       only: take_RK4_step
use msgwam_sinks_mod,     only: apply_breaking, apply_dissipation, &
                                check_boundaries
use msgwam_source_mod,    only: check_source, init_source
use msgwam_steady_mod,    only: get_steady_fluxes

implicit none
private

character(len=128) :: version = "msgwam.f90, 2025/04/15"
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
integer, dimension(5) :: clocks

! ==============================================================================
! mean state variables
! ==============================================================================

real, dimension(:, :, :), allocatable :: flux_x, flux_y, rho, N2, G2, u_bar, &
                                         v_bar, z_centers, z_faces

! ==============================================================================
! ray volume state variables
! ==============================================================================

type(t_ray), dimension(:, :, :), allocatable :: rays
integer, dimension(:, :, :), allocatable :: ghosts
integer, dimension(:, :), allocatable :: last_meta

! ==============================================================================

contains

subroutine init_clocks(clocks)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer, dimension(5), intent(out) :: clocks

    ! --------------------------------------------------------------------------

    clocks(1) = mpp_clock_id("      MS-GWaM mean state", grain=CLOCK_ROUTINE, &
        flags=MPP_CLOCK_SYNC)
    clocks(2) = mpp_clock_id("      MS-GWaM RK", grain=CLOCK_ROUTINE, &
        flags=MPP_CLOCK_SYNC)
    clocks(3) = mpp_clock_id("      MS-GWaM sinks", grain=CLOCK_ROUTINE, &
        flags=MPP_CLOCK_SYNC)
    clocks(4) = mpp_clock_id("      MS-GWaM source", grain=CLOCK_ROUTINE, &
        flags=MPP_CLOCK_SYNC)
    clocks(5) = mpp_clock_id("      MS-GWaM fluxes", grain=CLOCK_ROUTINE, &
        flags=MPP_CLOCK_SYNC)

end subroutine init_clocks

subroutine msgwam_init(lon_bounds, lat_bounds, p_ref, Time, axes)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:),    intent(in) :: lon_bounds, lat_bounds, p_ref
    type(time_type),       intent(in) :: Time
    integer, dimension(4), intent(in) :: axes

    if (is_initialized) then
        return
    end if

    call fms_init
    call time_manager_init
    call constants_init

    call init_msgwam_constants(lon_bounds, lat_bounds, p_ref)

    allocate(rho(q_max, i_max, j_max))
    allocate(N2(q_max, i_max, j_max))
    allocate(G2(q_max, i_max, j_max))

    allocate(u_bar(q_max, i_max, j_max))
    allocate(v_bar(q_max, i_max, j_max))

    allocate(z_centers(q_max + 2, i_max, j_max))
    allocate(z_faces(q_max + 1, i_max, j_max))

    allocate(flux_x(q_max + 1, i_max, j_max))
    allocate(flux_y(q_max + 1, i_max, j_max))

    allocate(rays(n_max, i_max, j_max))
    allocate(ghosts(n_source, i_max, j_max))
    allocate(last_meta(i_max, j_max))

    call init_source(p_ref, lat_bounds)
    call init_clocks(clocks)
    call init_nc_output(axes, Time)
    call init_ray_state(rays, ghosts, last_meta)

end subroutine msgwam_init

subroutine msgwam_calc(i_start, j_start, lat, &
    p_full, z_full, temp, uuu, vvv, &
    Time, dt, du_dt, dv_dt)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer, intent(in) :: i_start, j_start
    real, dimension(:, :), intent(in) :: lat
    real, dimension(i_max, j_max, q_max), intent(in)  :: p_full, z_full, temp, &
                                                         uuu, vvv
    type(time_type),                      intent(in)  :: Time
    real,                                 intent(in)  :: dt
    real, dimension(i_max, j_max, q_max), intent(out) :: du_dt, dv_dt

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: dt_rays

    ! --------------------------------------------------------------------------

    if (is_first_step) then
        is_first_step = .false.
        dt_rays = dt
    else
        dt_rays = dt / 2.
    end if

    call track_ray(rays, 1)
    call mpp_clock_begin(clocks(1))

    call update_mean_fields(z_full, p_full, temp, uuu, vvv, &
        z_centers, z_faces, u_bar, v_bar, rho, N2, G2)

    call mpp_clock_end(clocks(1))

    if (steady_state) then

        call get_steady_fluxes(z_centers, z_faces, u_bar, v_bar, rho, N2, &
            G2, dt, rays, ghosts, last_meta, flux_x, flux_y)

    else

        call track_ray(rays, 2)
        call mpp_clock_begin(clocks(2))

        call take_RK4_step(z_centers, u_bar, v_bar, N2, G2, dt_rays, rays)

        call mpp_clock_end(clocks(2))
        call track_ray(rays, 3)
        call mpp_clock_begin(clocks(3))

        call apply_dissipation(z_centers, z_faces, rho, dt_rays, rays)
        call apply_breaking(z_faces, rho, rays)
        call check_boundaries(z_centers, rays)

        call mpp_clock_end(clocks(3))
        call track_ray(rays, 4)
        call mpp_clock_begin(clocks(4))

        call check_source(z_centers, u_bar, v_bar, N2, G2, dt_rays, &
            rays, ghosts, last_meta)

        call mpp_clock_end(clocks(4))
        call track_ray(rays, 5)
        call mpp_clock_begin(clocks(5))

        call project_fluxes(z_centers, rays, flux_x, flux_y)
        call track_ray(rays, 6)
        call check_rays(rays)

        call mpp_clock_end(clocks(5))

    end if

    call get_accelerations(z_faces, rho, flux_x, flux_y, du_dt, dv_dt)
    call send_nc_output(i_start, j_start, Time, rho, N2, G2, flux_x, flux_y, &
        du_dt, dv_dt)

end subroutine msgwam_calc

subroutine msgwam_end

    call save_ray_state(rays, ghosts, last_meta)
    is_initialized = .false.

end subroutine msgwam_end

end module msgwam_mod