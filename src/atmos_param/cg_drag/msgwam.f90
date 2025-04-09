module msgwam_mod

! ==============================================================================
! This module implements MS-GWaM as described in Bölöni et al. (2021) and as
! implemented in dsconnelly/python-msgwam on GitHub.
! ==============================================================================

use constants_mod,    only: constants_init, PI, RDGAS
use diag_manager_mod, only: diag_manager_init, register_diag_field, send_data
use fms_mod,          only: check_nml_error, CLOCK_ROUTINE, close_file, &
                            error_mesg, FATAL, file_exist, fms_init, &
                            mpp_clock_begin, mpp_clock_end, mpp_clock_id, &
                            MPP_CLOCK_SYNC, mpp_pe, mpp_root_pe, &
                            open_namelist_file, stdlog, write_version_number
use time_manager_mod, only: time_manager_init, time_type

implicit none
private

character(len=128) :: version = "msgwam.f90, 2025/03/17"
character(len=128) :: tagname = "cayuga"

! ==============================================================================
! public interfaces
! ==============================================================================

public msgwam_calc, msgwam_end, msgwam_init

! ==============================================================================
! namelist
! ==============================================================================

real :: boundary_flux = 0.01
logical :: break_waves = .true.
real :: cp_center = 0.
real :: cp_max = 50.
real :: cp_width = 10.
real :: dr_source = 1000.
real :: epsilon = 0.
real :: lat_extrinsic = 15.
integer :: max_age = 10 * 86400
real :: min_flux = 1.e-8
real :: mu = 1.e-3
integer :: n_max = 2500
integer :: n_source = 48
real :: source_pressure = 300.e+2
real :: T_hat_source = 10. * 3600
logical :: use_shapiro_filter = .true.

! The scheme should eventually be rewritten to calculate these values from the
! mean flow, but for now the ray tracer treats them as constants. 
real :: H_rho = 8.e+3
real :: N0 = 0.015

namelist / msgwam_nml / &
    boundary_flux, break_waves, cp_center, cp_max, cp_width, dr_source, & 
    epsilon, lat_extrinsic, max_age, min_flux, mu, n_max, n_source, &
    source_pressure, T_hat_source, use_shapiro_filter, H_rho, N0

! ==============================================================================
! derived type definitions
! ==============================================================================

! Note that as written, the scheme omits dk and dl from the definition of the
! ray volume type. If this scheme were to include lateral propagation, these
! would need to be added. For now, the action is simply dens * dm.

type :: t_ray
    real :: r, dr, k, l, m, dm, dens
    integer :: age, meta, q_lo, q_mid, q_hi
    logical :: is_ghost
end type t_ray

! The t_inc type is used for increments to the ray volume properties that are
! evolved in time. Again, for now this includes only r and m.

type :: t_inc
    real :: r, m
end type t_inc

! ==============================================================================
! module status and timing variables
! ==============================================================================

logical :: is_initialized = .false.
integer, dimension(5) :: clocks

! ==============================================================================
! mean state variables
! ==============================================================================

real :: hgamma_sq
integer :: i_max, j_max, q_max, q_source
real, dimension(:), allocatable :: coriolis_sq
logical, dimension(:), allocatable :: extrinsic
real, dimension(:, :, :), allocatable :: z_faces, u_bar, v_bar, rho
real, dimension(:, :, :), allocatable :: z_padded

! ==============================================================================
! ray volume state variables
! ==============================================================================

type(t_ray), dimension(:, :, :), allocatable :: rays
type(t_inc), dimension(:, :, :), allocatable :: increments
real, dimension(:, :, :), allocatable :: flux_x, flux_y
integer, dimension(:, :, :), allocatable :: ghosts

! ==============================================================================
! ray volume source variables
! ==============================================================================

integer :: n_per_dir
logical :: is_stochastic
real :: dc_source, omega_hat_source
integer, dimension(:, :), allocatable :: last_meta
type(t_ray), dimension(:, :, :), allocatable :: launches
real, dimension(:), allocatable :: cp_source

real, dimension(4) :: cos_phi = (/ 1., 0., -1., 0. /)
real, dimension(4) :: sin_phi = (/ 0., 1., 0., -1. /)

! ==============================================================================
! RK3 coefficients
! ==============================================================================

real, dimension(3) :: As = (/ 0., -5. / 9., -153. / 128. /)
real, dimension(3) :: Bs = (/ 1. / 3., 15. / 16., 8. / 15. /)

! ==============================================================================
! netCDF input/ouput variables
! ==============================================================================

character(len=7) :: mod_name = "msgwam"
integer :: id_flux_x, id_flux_y, id_accel_x, id_accel_y
real, parameter :: missing_value = -999.

! ==============================================================================
! dispersion relation interfaces
! ==============================================================================

interface get_cg_r
    module procedure cg_r_from_ray
    module procedure cg_r_from_reals
end interface

interface get_omega_hat
    module procedure omega_hat_from_ray
    module procedure omega_hat_from_reals
end interface

contains

! ==============================================================================
! public subroutines
! ==============================================================================

subroutine msgwam_init(lon_bounds, lat_bounds, p_ref, Time, axes)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:),    intent(in) :: lon_bounds, lat_bounds, p_ref
    type(time_type),       intent(in) :: Time
    integer, dimension(4), intent(in) :: axes

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i_err, io, j, log_unit, nml_unit, q
    real :: lat

    ! --------------------------------------------------------------------------

    if (is_initialized) then
        return
    end if

    call fms_init
    call time_manager_init
    call constants_init

    if (file_exist("input.nml")) then
        nml_unit = open_namelist_file()
        i_err = 1

        do while (i_err /= 0)
            read(nml_unit, nml=msgwam_nml, iostat=io)
            i_err = check_nml_error(io, "msgwam_nml")
        end do

        call close_file(nml_unit)
    end if

    call write_version_number(version, tagname)
    log_unit = stdlog()

    if (mpp_pe() == mpp_root_pe()) then
        write (log_unit, nml=msgwam_nml)
    end if
    
    i_max = size(lon_bounds) - 1
    j_max = size(lat_bounds) - 1
    q_max = size(p_ref) - 1

    allocate(z_faces(q_max + 1, i_max, j_max))
    allocate(z_padded(q_max + 2, i_max, j_max))

    allocate(u_bar(q_max, i_max, j_max))
    allocate(v_bar(q_max, i_max, j_max))
    allocate(rho(q_max, i_max, j_max))

    allocate(flux_x(q_max + 1, i_max, j_max))
    allocate(flux_y(q_max + 1, i_max, j_max))

    allocate(rays(n_max, i_max, j_max))
    allocate(ghosts(n_source, i_max, j_max))
    allocate(last_meta(i_max, j_max))

    call init_ray_state(rays, ghosts, last_meta)

    allocate(increments(n_max, i_max, j_max))
    allocate(launches(n_source, i_max, j_max))

    do q = 1, q_max
        if (p_ref(q) > source_pressure) then
            q_source = q
            exit
        end if
    end do

    allocate(coriolis_sq(j_max))
    allocate(extrinsic(j_max))

    do j = 1, j_max
        lat = 0.5 * (lat_bounds(j) + lat_bounds(j + 1))
        coriolis_sq(j) = (2 * PI * sin(lat) / 86400.) ** 2
        extrinsic = abs(180 * lat / PI) > lat_extrinsic
    end do

    hgamma_sq = ((1. / 2. - 2. / 7.) / H_rho) ** 2

    call init_source
    call init_nc_output(axes, Time)
    is_stochastic = epsilon > 0.

    clocks(1) = mpp_clock_id("      MS-GWaM mean state", grain=CLOCK_ROUTINE, flags=MPP_CLOCK_SYNC)
    clocks(2) = mpp_clock_id("      MS-GWaM RK3", grain=CLOCK_ROUTINE, flags=MPP_CLOCK_SYNC)
    clocks(3) = mpp_clock_id("      MS-GWaM sinks", grain=CLOCK_ROUTINE, flags=MPP_CLOCK_SYNC)
    clocks(4) = mpp_clock_id("      MS-GWaM source", grain=CLOCK_ROUTINE, flags=MPP_CLOCK_SYNC)
    clocks(5) = mpp_clock_id("      MS-GWaM fluxes", grain=CLOCK_ROUTINE, flags=MPP_CLOCK_SYNC)

end subroutine msgwam_init

subroutine msgwam_calc(i_start, j_start, lat, &
    p_full, z_full, temp, uuu, vvv, &
    Time, dt, du_dt, dv_dt)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                  intent(in)  :: i_start, j_start
    real, dimension(:, :),    intent(in)  :: lat
    real, dimension(:, :, :), intent(in)  :: p_full, z_full, temp, uuu, vvv
    type(time_type),          intent(in)  :: Time
    real,                     intent(in)  :: dt
    real, dimension(:, :, :), intent(out) :: du_dt, dv_dt

    ! --------------------------------------------------------------------------

    call mpp_clock_begin(clocks(1))

    call update_mean_state(z_full, p_full, temp, uuu, vvv, &
        z_padded, z_faces, u_bar, v_bar, rho)

    call mpp_clock_end(clocks(1))
    call mpp_clock_begin(clocks(2))

    call take_RK3_step(z_padded, u_bar, v_bar, dt / 2., rays, increments)

    call mpp_clock_end(clocks(2))
    call mpp_clock_begin(clocks(3))

    call apply_dissipation(z_padded, rho, dt / 2., rays)
    call apply_breaking(z_faces, rho, rays)

    call mpp_clock_end(clocks(3))
    call mpp_clock_begin(clocks(4))

    call check_boundaries(z_padded, rays)
    call check_source(z_padded, u_bar, v_bar, dt / 2., rays, ghosts, &
        last_meta, launches)

    call mpp_clock_end(clocks(4))
    call mpp_clock_begin(clocks(5))

    call update_fluxes(z_padded, rays, flux_x, flux_y)
    call get_accelerations(z_faces, rho, flux_x, flux_y, du_dt, dv_dt)
    call send_nc_output(i_start, j_start, Time, flux_x, flux_y, du_dt, dv_dt)

    call mpp_clock_end(clocks(5))

end subroutine msgwam_calc

subroutine msgwam_end

    call save_ray_state(rays, ghosts, last_meta)
    is_initialized = .false.

end subroutine msgwam_end

! ==============================================================================
! dispersion relation subroutines
! ==============================================================================

pure function cg_r_from_ray(ray, f2) result(cg_r)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    type(t_ray), intent(in) :: ray
    real,        intent(in) :: f2
    real                    :: cg_r

    ! --------------------------------------------------------------------------

    cg_r = cg_r_from_reals(ray%k, ray%l, ray%m, f2)

end function cg_r_from_ray

pure function cg_r_from_reals(k, l, m, f2) result(cg_r)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: k, l, m, f2
    real             :: cg_r

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: omega_hat, wvn_sq

    ! --------------------------------------------------------------------------

    wvn_sq = k ** 2 + l ** 2 + m ** 2 + hgamma_sq
    omega_hat = get_omega_hat(k, l, m, f2)

    cg_r = -m * (omega_hat ** 2 - f2) / (omega_hat * wvn_sq)

end function cg_r_from_reals

pure function get_dm(m) result(dm)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: m
    real             :: dm

    ! --------------------------------------------------------------------------

    dm = dc_source * (m ** 2) / N0

end function get_dm

pure function get_m(k, l, f2) result(m)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: k, l, f2
    real             :: m

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: omega_hat_sq

    ! --------------------------------------------------------------------------

    omega_hat_sq = omega_hat_source ** 2
    m = -sqrt(&
        (k ** 2 + l ** 2) * (N0 ** 2 - omega_hat_sq) / &
        (omega_hat_sq - f2) &
    )

end function get_m

pure function omega_hat_from_ray(ray, f2) result(omega_hat)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    type(t_ray), intent(in) :: ray
    real,        intent(in) :: f2
    real                    :: omega_hat

    ! --------------------------------------------------------------------------

    omega_hat = omega_hat_from_reals(ray%k, ray%l, ray%m, f2)

end function omega_hat_from_ray

pure function omega_hat_from_reals(k, l, m, f2) result(omega_hat)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: k, l, m, f2
    real             :: omega_hat

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: wvn_hor_sq, m2

    ! --------------------------------------------------------------------------

    wvn_hor_sq = k ** 2 + l ** 2
    m2 = m ** 2 + hgamma_sq

    omega_hat = sqrt( &
        (N0 ** 2 * wvn_hor_sq + f2 * m2) / &
        (wvn_hor_sq + m2) &
    )

end function omega_hat_from_reals

! ==============================================================================
! mean state helpers
! ==============================================================================

pure subroutine shapiro_filter(profile)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :), intent(inout) :: profile

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: q_max
    real, dimension(4, size(profile, 2), size(profile, 3)) :: tmp

    ! --------------------------------------------------------------------------

    q_max = size(profile, 1)
    tmp(1:2, :, :) = profile(1:2, :, :)
    tmp(3:4, :, :) = profile((q_max - 1):q_max, :, :)

    profile(2:(q_max - 1), :, :) = ( &
        profile(:(q_max - 2), :, :) + &
        2 * profile(2:(q_max - 1), :, :) + &
        profile(3:, :, :) &
    ) / 4.

    profile(1, :, :) = (3 * tmp(1, :, :) + tmp(2, :, :)) / 4.
    profile(q_max, :, :) = (tmp(3, :, :) + 3 * tmp(4, :, :)) / 4.

end subroutine shapiro_filter

pure subroutine get_accelerations(z_faces, rho, flux_x, flux_y, du_dt, dv_dt)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :), intent(in)  :: z_faces, rho, flux_x, flux_y
    real, dimension(:, :, :), intent(out) :: du_dt, dv_dt

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, q
    real :: dz, dFx_dz, dFy_dz

    ! --------------------------------------------------------------------------

    do j = 1, j_max
        do i = 1, i_max
            do q = 1, q_max
                dz = z_faces(q, i, j) - z_faces(q + 1, i, j)
                dFx_dz = (flux_x(q, i, j) - flux_x(q + 1, i, j)) / dz
                dFy_dz = (flux_y(q, i, j) - flux_y(q + 1, i, j)) / dz

                du_dt(i, j, q) = -dFx_dz / rho(q, i, j)
                dv_dt(i, j, q) = -dFy_dz / rho(q, i, j)
            end do
        end do
    end do

end subroutine

pure subroutine update_fluxes(z_padded, rays, flux_x, flux_y)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:, :, :),       intent(in)  :: z_padded
    type(t_ray), dimension(:, :, :), intent(in)  :: rays
    real, dimension(:, :, :),        intent(out) :: flux_x, flux_y

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n, q
    real :: cg, dz, frac, f_x, f_y, r_hi, r_lo, z_hi, z_lo

    ! --------------------------------------------------------------------------

    flux_x = 0.
    flux_y = 0.

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_max

                if (rays(n, i, j)%meta == -1) then
                    cycle
                end if

                associate (ray => rays(n, i, j))
                    r_lo = ray%r - ray%dr / 2.
                    r_hi = ray%r + ray%dr / 2.

                    cg = get_cg_r(ray, coriolis_sq(j))
                    f_x = ray%k * ray%dens * ray%dm * cg
                    f_y = ray%l * ray%dens * ray%dm * cg

                    do q = ray%q_hi, ray%q_lo
                        z_hi = z_padded(q, i, j)
                        z_lo = z_padded(q + 1, i, j)
                        dz = z_hi - z_lo
                    
                        frac = (min(r_hi, z_hi) - max(r_lo, z_lo)) / dz
                        flux_x(q, i, j) = flux_x(q, i, j) + frac * f_x
                        flux_y(q, i, j) = flux_y(q, i, j) + frac * f_y
                    end do

                end associate

            end do
        end do
    end do

    if (use_shapiro_filter) then
        call shapiro_filter(flux_x)
        call shapiro_filter(flux_y)
    end if

end subroutine update_fluxes

pure subroutine update_mean_state(z_full, p_full, temp, uuu, vvv, &
    z_padded, z_faces, u_bar, v_bar, rho)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :),  intent(in)  :: z_full, p_full, temp, uuu, vvv
    real, dimension(0:, :, :), intent(out) :: z_padded
    real, dimension(:, :, :),  intent(out) :: z_faces, u_bar, v_bar, rho

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, q

    ! --------------------------------------------------------------------------

    do q = 1, q_max
        do j = 1, j_max
            do i = 1, i_max
                u_bar(q, i, j) = uuu(i, j, q)
                v_bar(q, i, j) = vvv(i, j, q)
                rho(q, i, j) = p_full(i, j, q) / RDGAS / temp(i, j, q)
                z_padded(q, i, j) = z_full(i, j, q)
            end do
        end do
    end do

    z_padded(0, :, :) = 2 * z_padded(1, :, :) - z_padded(2, :, :)
    z_padded(q_max + 1, :, :) = 2 * z_padded(q_max, :, :) - &
        z_padded(q_max - 1, :, :)

    z_faces = (z_padded(1:, :, :) + z_padded(:q_max, :, :)) / 2.

end subroutine update_mean_state

! ==============================================================================
! ray volume state and sort helpers
! ==============================================================================

subroutine init_ray_state(rays, ghosts, last_meta)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    type(t_ray), dimension(:, :, :), intent(out) :: rays
    integer, dimension(:, :, :),     intent(out) :: ghosts
    integer, dimension(:, :),        intent(out) :: last_meta

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: iostat, unit
    character(len=32) :: fname, pe_str
    logical :: from_restart

    ! --------------------------------------------------------------------------

    write(pe_str, "(I2.2)") mpp_pe()
    fname = "INPUT/rays-" // trim(pe_str) // ".dat"
    inquire(file=trim(fname), exist=from_restart)

    if (.not. from_restart) then
        rays(:, :, :)%meta = -1
        ghosts(:, :, :) = -1
        last_meta(:, :) = 1
    
    else
        open(newunit=unit, file=trim(fname), form="unformatted", &
            iostat=iostat, action="read")

        if (iostat /= 0) then
            call error_mesg("msgwam_mod", "error loading ray state", FATAL)
        end if

        read(unit) rays
        read(unit) ghosts
        read(unit) last_meta

        close(unit)
    end if

end subroutine init_ray_state

pure subroutine delete_ray(n, i, j, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                         intent(in)    :: n, i, j
    type(t_ray), dimension(:, :, :), intent(inout) :: rays

    ! --------------------------------------------------------------------------

    associate (ray => rays(n, i, j))
        ray%r = 0
        ray%dr = 0
        
        ray%k = 0
        ray%l = 0
        ray%m = 0

        ray%dm = 0
        ray%dens = 0

        ray%age = 0
        ray%meta = -1

        ray%q_lo = -1
        ray%q_mid = -1
        ray%q_hi = -1

        ray%is_ghost = .false.
    end associate

end subroutine delete_ray

subroutine save_ray_state(rays, ghosts, last_meta)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    type(t_ray), dimension(:, :, :), intent(in) :: rays
    integer, dimension(:, :, :),     intent(in) :: ghosts
    integer, dimension(:, :),        intent(in) :: last_meta

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: iostat, unit
    character(len=32) :: fname, pe_str

    ! --------------------------------------------------------------------------

    write(pe_str, "(I2.2)") mpp_pe()
    fname = "RESTART/rays-" // trim(pe_str) // ".dat"

    open(newunit=unit, file=fname, form="unformatted", &
        iostat=iostat, action="write")

    if (iostat /= 0) then
        call error_mesg("msgwam_mod", "error saving ray state", FATAL)
    end if

    write(unit) rays
    write(unit) ghosts
    write(unit) last_meta

    close(unit)

end subroutine save_ray_state

pure function locate(r, z, q_guess) result(q)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real,               intent(in) :: r
    real, dimension(:), intent(in) :: z
    integer,            intent(in) :: q_guess
    integer                        :: q

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: q_max

    ! --------------------------------------------------------------------------

    q_max = size(z)
    if (r > z(1)) then
        q = 1
        return
    end if

    if (r < z(q_max)) then
        q = q_max - 1
        return
    end if
    
    q = q_guess
    do while (.true.)
        if (z(q) < r) then
            q = q - 1
        else if (z(q + 1) .ge. r) then
            q = q + 1
        else
            exit
        end if
    end do

end function locate

! ==============================================================================
! ray volume source subroutines
! ==============================================================================

subroutine check_source(z_padded, u_bar, v_bar, dt, &
    rays, ghosts, last_meta, launches)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:, :, :),       intent(in)    :: z_padded
    real, dimension(:, :, :),        intent(in)    :: u_bar, v_bar
    real,                            intent(in)    :: dt
    type(t_ray), dimension(:, :, :), intent(inout) :: rays
    integer, dimension(:, :, :),     intent(inout) :: ghosts
    integer, dimension(:, :),        intent(inout) :: last_meta
    type(t_ray), dimension(:, :, :), intent(out)   :: launches

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: add_at, i, j, n, s
    integer, dimension(size(rays, 2), size(rays, 3)) :: n_excess

    ! --------------------------------------------------------------------------

    do j = 1, j_max
        do i = 1, i_max
            n_excess(i, j) = 0

            do n = 1, n_max
                if (rays(n, i, j)%meta /= -1) then
                    n_excess(i, j) = n_excess(i, j) + 1
                end if
            end do
        end do
    end do

    call get_launches(z_padded, u_bar, v_bar, dt, rays, ghosts, &
        last_meta, launches, n_excess)

    n_excess = max(n_excess - n_max, 0)
    call prune(n_excess, rays)

    do j = 1, j_max
        do i = 1, i_max
            add_at = 1

            do n = 1, n_source
                if (launches(n, i, j)%meta == -1) then
                    cycle
                end if

                do s = add_at, n_max
                    if (rays(s, i, j)%meta == -1) then
                        add_at = s
                        exit
                    end if
                end do

                if (rays(add_at, i, j)%meta /= -1) then
                    call error_mesg("msgwam_mod", "too many rays", FATAL)
                end if

                rays(add_at, i, j) = launches(n, i, j)

                if (.not. is_stochastic) then
                    if (ghosts(n, i, j) > 0) then
                        rays(ghosts(n, i, j), i, j)%is_ghost = .false.
                    end if

                    ghosts(n, i, j) = add_at
                end if

            end do
        end do
    end do

end subroutine check_source

subroutine init_source

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: n

    ! --------------------------------------------------------------------------

    n_per_dir = n_source / 4
    dc_source = cp_max / n_per_dir
    allocate(cp_source(n_per_dir))

    do n = 1, n_per_dir
        cp_source(n) = (n - 0.5) * dc_source
    end do

    omega_hat_source = 2 * PI / T_hat_source

end subroutine init_source

pure subroutine find_lowest_energies(n_find, rays, idx)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer, dimension(:, :),        intent(in)  :: n_find
    type(t_ray), dimension(:, :, :), intent(in)  :: rays
    integer, dimension(:, :, :),     intent(out) :: idx

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: omega_hat, energy
    integer :: add_at, i, j, n, s
    real, dimension(size(idx, 1), size(rays, 2), size(rays, 3)) :: lowest

    ! --------------------------------------------------------------------------

    lowest = huge(lowest)

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_max
                associate (ray => rays(n, i, j))
                    if ((ray%meta == -1) .or. ray%is_ghost) then
                        cycle
                    end if

                    omega_hat = get_omega_hat(ray, coriolis_sq(j))
                    energy = ray%dens * ray%dm * omega_hat
                end associate

                add_at = -1
                do s = 1, n_find(i, j)
                    if (energy < lowest(s, i, j)) then
                        add_at = s
                    else
                        exit
                    end if
                end do

                if (add_at == -1) then
                    cycle
                end if

                do s = 1, add_at - 1
                    lowest(s, i, j) = lowest(s + 1, i, j)
                    idx(s, i, j) = idx(s + 1, i, j)
                end do

                lowest(add_at, i, j) = energy
                idx(add_at, i, j) = n
            end do
        end do
    end do

end subroutine find_lowest_energies

pure subroutine prune(n_excess, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer, dimension(:, :),        intent(in)    :: n_excess
    type(t_ray), dimension(:, :, :), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n
    integer, dimension(:, :, :), allocatable :: idx

    ! --------------------------------------------------------------------------

    allocate(idx(maxval(n_excess), i_max, j_max))
    call find_lowest_energies(n_excess, rays, idx)

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_excess(i, j)
                call delete_ray(idx(n, i, j), i, j, rays)
            end do
        end do
    end do

end subroutine prune

subroutine get_launches(z_padded, u_bar, v_bar, dt, rays, ghosts, &
    last_meta, launches, n_added)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:, :, :),       intent(in)    :: z_padded
    real, dimension(:, :, :),        intent(in)    :: u_bar, v_bar
    real,                            intent(in)    :: dt
    type(t_ray), dimension(:, :, :), intent(in)    :: rays
    integer, dimension(:, :, :),     intent(in)    :: ghosts
    integer, dimension(:, :),        intent(inout) :: last_meta
    type(t_ray), dimension(:, :, :), intent(out)   :: launches
    integer, dimension(:, :),        intent(out)   :: n_added

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    logical :: proceed
    integer :: dir, i, j, n, s, q_hi, q_mid, q_lo
    real :: cg, cp, flux, k, l, m, mag_cp_hat, mag_wvn_hor, r, total, u, v
    real, dimension(size(launches, 1), size(launches, 2), size(launches, 3)) &
        :: rand

    ! --------------------------------------------------------------------------

    call random_number(rand)

    do j = 1, j_max
        do i = 1, i_max
            r = z_padded(q_source, i, j)
            q_hi = locate(r, z_padded(1:, i, j), q_source)
            q_lo = locate(r - dr_source, z_padded(1:, i, j), q_source)
            q_mid = locate(r - dr_source / 2., z_padded(1:, i, j), q_source)

            total = 0.
            do dir = 1, 4
                do n = 1, n_per_dir
                    s = (dir - 1) * n_per_dir + n
                    proceed = .true.

                    if (.not. is_stochastic) then
                        if (ghosts(s, i, j) > 0) then

                            associate(ray => rays(ghosts(s, i, j), i, j))
                                if ((ray%r - ray%dr / 2.) < r) then
                                    launches(s, i, j)%meta = -1
                                    proceed = .false.
                                end if
                            end associate

                        end if
                    end if

                    cp = cp_source(n)
                    u = u_bar(q_source, i, j)
                    v = v_bar(q_source, i, j)
                    mag_cp_hat = cp - cos_phi(dir) * u - sin_phi(dir) * v

                    if (extrinsic(j)) then
                        flux = exp(-0.5 * ((cp / cp_width) ** 2))
                    else
                        flux = exp(-0.5 * ((mag_cp_hat / cp_width) ** 2))
                    end if

                    total = total + flux

                    if (.not. proceed) then
                        cycle
                    end if

                    mag_wvn_hor = omega_hat_source / mag_cp_hat
                    k = mag_wvn_hor * cos_phi(dir)
                    l = mag_wvn_hor * sin_phi(dir)

                    m = get_m(k, l, coriolis_sq(j))
                    cg = get_cg_r(k, l, m, coriolis_sq(j))

                    if (is_stochastic) then
                        proceed = rand(s, i, j) < epsilon * cg * dt / dr_source
                    end if
                        
                    if (proceed) then
                        launches(s, i, j)%r = r - dr_source / 2.
                        launches(s, i, j)%dr = dr_source

                        launches(s, i, j)%k = k
                        launches(s, i, j)%l = l
                        launches(s, i, j)%m = m

                        launches(s, i, j)%dm  = get_dm(m)
                        launches(s, i, j)%dens = flux / abs( &
                            mag_wvn_hor * launches(s, i, j)%dm * cg &
                        )

                        launches(s, i, j)%age = 0
                        launches(s, i, j)%meta = last_meta(i, j)
                        last_meta(i, j) = last_meta(i, j) + 1
                        n_added(i, j) = n_added(i, j) + 1

                        launches(s, i, j)%q_mid = q_mid
                        launches(s, i, j)%q_hi = q_hi
                        launches(s, i, j)%q_lo = q_lo

                        launches(s, i, j)%is_ghost = .not. is_stochastic
                    else
                        launches(s, i, j)%meta = -1
                    end if
                end do
            end do

            launches(:, i, j)%dens = launches(:, i, j)%dens &
                * 2 * boundary_flux / total

        end do
    end do

end subroutine get_launches

! ==============================================================================
! ray volume sink subroutines
! ==============================================================================

pure subroutine apply_breaking(z_faces, rho, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :),        intent(in)    :: z_faces, rho
    type(t_ray), dimension(:, :, :), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n, q
    real :: dz, frac, max_kappa, omega_hat, r_hi, r_lo, S, z_hi, z_lo
    real, dimension(size(rays, 1), size(rays, 2), size(rays, 3)) :: wvn_sq
    real, dimension(size(rho, 1), size(rho, 2), size(rho, 3)) :: num, den, kappa

    ! --------------------------------------------------------------------------

    if (.not. break_waves) then
        return
    end if

    num = -0.5 * rho * N0 ** 2
    den = 0.

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_max

                if (rays(n, i, j)%meta == -1) then
                    cycle
                end if

                associate (ray => rays(n, i, j))
                    omega_hat = get_omega_hat(ray, coriolis_sq(j))
                    wvn_sq(n, i, j) = ray%k ** 2 + ray%l ** 2 + ray%m ** 2
                    S = ray%m ** 2 * omega_hat * ray%dens * ray%dm

                    r_lo = ray%r - ray%dr / 2.
                    r_hi = ray%r + ray%dr / 2.

                    do q = max(1, ray%q_hi - 1), min(q_max, ray%q_lo + 1)
                        z_hi = z_faces(q, i, j)
                        z_lo = z_faces(q + 1, i, j)
                        dz = z_hi - z_lo

                        if (r_hi < z_lo) then
                            cycle
                        end if

                        if (r_lo > z_hi) then
                            exit
                        end if
                    
                        frac = (min(r_hi, z_hi) - max(r_lo, z_lo)) / dz
                        num(q, i, j) = num(q, i, j) + frac * S
                        den(q, i, j) = den(q, i, j) + frac * S * wvn_sq(n, i, j)
                    end do
                end associate

            end do
        end do
    end do

    kappa = merge(num / den, 0., den /= 0.)

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_max

                if (rays(n, i, j)%meta == -1) then
                    cycle
                end if

                associate(ray => rays(n, i, j))
                    r_lo = ray%r - ray%dr / 2.
                    r_hi = ray%r + ray%dr / 2.
                
                    max_kappa = 0.
                    do q = max(1, ray%q_hi - 1), min(q_max, ray%q_lo + 1)
                        z_hi = z_faces(q, i, j)
                        z_lo = z_faces(q + 1, i, j)

                        if (r_hi < z_lo) then
                            cycle
                        end if

                        if (r_lo > z_hi) then
                            exit
                        end if

                        if (kappa(q, i, j) > max_kappa) then
                            max_kappa = kappa(q, i, j)
                        end if
                    end do

                    ray%dens = ray%dens * max(0., &
                        1 - wvn_sq(n, i, j) * max_kappa &
                    )
                end associate

            end do
        end do
    end do

end subroutine apply_breaking

pure subroutine apply_dissipation(z_padded, rho, dt, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:, :, :),       intent(in)    :: z_padded
    real, dimension(:, :, :),        intent(in)    :: rho
    real,                            intent(in)    :: dt
    type(t_ray), dimension(:, :, :), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n
    real :: damping, omega_hat, nu, nu_hi, nu_lo, wvn_sq, z_hi, z_lo

    ! --------------------------------------------------------------------------

    if (mu == 0.) then
        return
    end if

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_max

                if (rays(n, i, j)%meta == -1) then
                    cycle
                end if

                associate (ray => rays(n, i, j))
                    z_hi = z_padded(ray%q_mid, i, j)
                    z_lo = z_padded(ray%q_mid + 1, i, j)

                    nu_hi = mu / rho(ray%q_mid, i, j)
                    nu_lo = mu / rho(ray%q_mid + 1, i, j)

                    nu = ( &
                        (nu_lo * (z_hi - ray%r) + nu_hi * (ray%r - z_lo)) / &
                        (z_hi - z_lo) &
                    )

                    wvn_sq = ray%k ** 2 + ray%l ** 2 + ray%m ** 2
                    omega_hat = get_omega_hat(ray, coriolis_sq(j))

                    damping = nu * wvn_sq * ( &
                        1 + coriolis_sq(j) / omega_hat ** 2 &
                    )

                    ray%dens = ray%dens * exp(-dt * damping)
                end associate

            end do
        end do
    end do

end subroutine apply_dissipation

 subroutine check_boundaries(z_padded, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:, :, :),       intent(in)    :: z_padded
    type(t_ray), dimension(:, :, :), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    logical :: delete
    integer :: i, j, n
    real :: cg, flux, wvn, z_lo, z_hi

    ! --------------------------------------------------------------------------

    do j = 1, j_max
        do i = 1, i_max
            z_lo = z_padded(q_max + 1, i, j)
            z_hi = z_padded(0, i, j)

            do n = 1, n_max
                associate (ray => rays(n, i, j))
                    if (ray%meta == -1 .or. ray%is_ghost) then
                        cycle
                    end if

                    cg = get_cg_r(ray, coriolis_sq(j))
                    wvn = sqrt(ray%k ** 2 + ray%l ** 2)
                    flux = wvn * ray%dens * ray%dm * cg

                    delete = ray%r + ray%dr / 2. < z_lo
                    delete = delete .or. (ray%r - ray%dr / 2.) > z_hi
                    delete = delete .or. (abs(flux) < min_flux)
                    delete = delete .or. (ray%age > max_age)
                end associate

                if (delete) then
                    call delete_ray(n, i, j, rays)
                end if
            
            end do
        end do
    end do

end subroutine check_boundaries

! ==============================================================================
! time stepping subroutines
! ==============================================================================

pure subroutine take_RK3_step(z_padded, u_bar, v_bar, dt, rays, increments)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:, :, :),       intent(in)    :: z_padded
    real, dimension(:, :, :),        intent(in)    :: u_bar, v_bar
    real,                            intent(in)    :: dt
    type(t_ray), dimension(:, :, :), intent(inout) :: rays
    type(t_inc), dimension(:, :, :), intent(out)   :: increments

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: dm_dt, dr_dt, dz
    integer i, j, n, q, stage
    real, dimension(size(u_bar, 1) - 1, size(u_bar, 2), size(u_bar, 3)) :: &
        du_dr, dv_dr

    ! --------------------------------------------------------------------------

    do j = 1, j_max
        do i = 1, i_max
            do q = 1, q_max - 1
                dz = z_padded(q, i, j) - z_padded(q + 1, i, j)
                du_dr(q, i, j) = (u_bar(q, i, j) - u_bar(q + 1, i, j)) / dz
                dv_dr(q, i, j) = (v_bar(q, i, j) - v_bar(q + 1, i, j)) / dz
            end do
        end do
    end do

    do stage = 1, 3

        do j = 1, j_max
            do i = 1, i_max
                do n = 1, n_max

                    if (rays(n, i, j)%meta == -1) then
                        cycle
                    end if

                    associate( &
                        ray => rays(n, i, j), &
                        inc => increments(n, i, j), &
                        z => z_padded(1:, i, j) &
                    )

                        dr_dt = get_cg_r(ray, coriolis_sq(j))
                        dm_dt = -( &
                            du_dr(ray%q_mid, i, j) * ray%k + &
                            dv_dr(ray%q_mid, i, j) * ray%l &
                        )

                        inc%r = As(stage) * inc%r + dt * dr_dt
                        ray%r = ray%r + Bs(stage) * inc%r

                        if (.not. ray%is_ghost) then
                            inc%m = As(stage) * inc%m + dt * dm_dt
                            ray%m = ray%m + Bs(stage) * inc%m
                        end if

                        ray%q_mid = locate(ray%r, z, ray%q_mid)

                        if (stage == 3) then
                            ray%q_hi = locate(ray%r + ray%dr / 2., z, ray%q_hi)
                            ray%q_lo = locate(ray%r - ray%dr / 2., z, ray%q_lo)

                            ray%age = ray%age + dt
                        end if

                    end associate

                end do
            end do
        end do

    end do

end subroutine take_RK3_step

! ==============================================================================
! netCDF output subroutines
! ==============================================================================

subroutine init_nc_output(axes, Time)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer, dimension(4), intent(in) :: axes
    type(time_type),       intent(in) :: Time

    ! --------------------------------------------------------------------------

    id_flux_x = register_diag_field(mod_name, "gw_flux_x", axes(1:3), Time, &
        "zonal GW flux", "Pa", missing_value=missing_value)
    id_flux_y = register_diag_field(mod_name, "gw_flux_y", axes(1:3), Time, &
        "meridional GW flux", "Pa", missing_value=missing_value)

    id_accel_x = register_diag_field(mod_name, "gw_accel_x", axes(1:3), Time, &
        "zonal GW acceleration", "m/s^2", missing_value=missing_value)
    id_accel_y = register_diag_field(mod_name, "gw_accel_y", axes(1:3), Time, &
        "meridional GW acceleration", "m/s^2", missing_value=missing_value)

end subroutine init_nc_output

pure subroutine reorder_axes(a, b)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :), intent(in)  :: a
    real, dimension(:, :, :), intent(out) :: b

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, q

    ! --------------------------------------------------------------------------

    do j = 1, j_max
        do i = 1, i_max
            do q = 1, q_max
                b(i, j, q) = a(q, i, j)
            end do
        end do
    end do

end subroutine reorder_axes

subroutine send_nc_output(i_start, j_start, Time, flux_x, flux_y, du_dt, dv_dt)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                  intent(in) :: i_start, j_start
    type(time_type),          intent(in) :: Time
    real, dimension(:, :, :), intent(in) :: flux_x, flux_y, du_dt, dv_dt

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i_err
    real, dimension(size(flux_x, 2), size(flux_x, 3), size(flux_x, 1) - 1) :: &
        tmp

    ! --------------------------------------------------------------------------

    if (id_flux_x > 0) then
        call reorder_axes(flux_x(2:, :, :), tmp)
        i_err = send_data(id_flux_x, tmp, Time, i_start, j_start)
    end if

    if (id_flux_y > 0) then
        call reorder_axes(flux_y(2:, :, :), tmp)
        i_err = send_data(id_flux_y, tmp, Time, i_start, j_start)
    end if

    if (id_accel_x > 0) then
        i_err = send_data(id_accel_x, du_dt, Time, i_start, j_start)
    end if

    if (id_accel_y > 0) then
        i_err = send_data(id_accel_y, dv_dt, Time, i_start, j_start)
    end if

end subroutine send_nc_output

end module msgwam_mod
