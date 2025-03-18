module cg_drag_mod

! ==============================================================================
! This version of cg_drag_mod implements MS-GWaM as described in Bölöni et al.
! (2021) and implemented in dsconnelly/python-msgwam on GitHub. A future version
! of this code should allow switching gravity wave schemes in the namelist.
! ==============================================================================

use constants_mod,    only: constants_init, PI, RDGAS
use diag_manager_mod, only: diag_manager_init, register_diag_field, send_data
use fms_mod,          only: check_nml_error, close_file, error_mesg, FATAL, &
                            file_exist, fms_init, mpp_pe, mpp_root_pe, &
                            open_namelist_file, stdlog, write_version_number
use time_manager_mod, only: time_manager_init, time_type

implicit none
private

character(len=128) :: version = "msgwam.f90, 2025/03/17"
character(len=128) :: tagname = "cayuga"

! ==============================================================================
! interfaces
! ==============================================================================

public cg_drag_calc, cg_drag_end, cg_drag_init

! ==============================================================================
! namelist
! ==============================================================================

real :: boundary_flux = 0.01
real :: dr_source = 1000.
real :: cp_center = 15.
real :: cp_max = 50.
real :: cp_width = 10.
real :: dk_source = 0.0001
real :: dl_source = 0.0001
real :: epsilon = 1.
logical :: extrinsic = .true.
integer :: max_age = 10 * 86400
real :: min_flux = 1.e-8
integer :: n_max = 2500
integer :: n_source = 48
real :: padding_z = 500.
real :: source_pressure = 300.e+2
real :: T_hat_source = 10. * 3600

! The scheme should eventually be rewritten to calculate these values from the
! mean flow, but for now the ray tracer treats them as constants. 
real :: H_rho = 8.e+3
real :: N0 = 0.015

namelist / cg_drag_nml / &
    boundary_flux, dr_source, cp_center, cp_max, cp_width, dk_source, &
    dl_source, epsilon, extrinsic, max_age, min_flux, n_max, n_source, &
    padding_z, source_pressure, T_hat_source, H_rho, N0

! ==============================================================================
! module-level private variables
! ==============================================================================

type :: t_inc
    real :: r, m
end type t_inc

type :: t_ray
    real :: r, dr, k, l, m, dk, dl, dm, dens
    integer :: age, meta
end type t_ray

type :: t_tend
    real :: r, m
end type t_tend

interface operator (+)
    module procedure add_inc_inc
    module procedure add_ray_inc
end interface

interface operator (*)
    module procedure mult_scalar_inc
    module procedure mult_scalar_tend
end interface

interface get_cg_r
    module procedure cg_r_from_ray
    module procedure cg_r_from_reals
end interface

interface get_omega_hat
    module procedure omega_hat_from_ray
    module procedure omega_hat_from_reals
end interface

type(t_ray), dimension(:, :, :), allocatable :: rays
type(t_tend), dimension(:, :, :), allocatable :: drays_dt
type(t_inc), dimension(:, :, :), allocatable :: increments
type(t_ray), dimension(:, :, :), allocatable :: launches

integer :: n_per_dir
real :: dc_source, omega_hat_source
real, dimension(:), allocatable :: cp_source, flux_source
integer, dimension(:, :), allocatable :: last_meta

real, dimension(4) :: cos_phi = (/ 1., 0., -1., 0. /)
real, dimension(4) :: sin_phi = (/ 0., 1., 0., -1. /)

real, dimension(:, :, :), allocatable :: flux_x, flux_y, z_faces, z_padded

real :: hgamma
integer :: q_source
real, dimension(:), allocatable :: coriolis
integer :: i_max, j_max, q_max

real, dimension(3) :: As = (/ 0., -5. / 9., -153. / 128. /)
real, dimension(3) :: Bs = (/ 1. / 3., 15. / 16., 8. / 15. /)

character(len=7) :: mod_name = "cg_drag"
integer :: id_flux_x, id_flux_y, id_accel_x, id_accel_y
real, parameter :: missing_value = -999.

logical :: is_initialized = .false.

contains

! ==============================================================================
! public subroutines
! ==============================================================================

subroutine cg_drag_init(lon_bounds, lat_bounds, p_ref, Time, axes)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:),    intent(in) :: lon_bounds, lat_bounds, p_ref
    type(time_type),       intent(in) :: Time
    integer, dimension(4), intent(in) :: axes

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i_err, io, log_unit, nml_unit
    integer :: j, q
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
            read(nml_unit, nml=cg_drag_nml, iostat=io)
            i_err = check_nml_error(io, "cg_drag_nml")
        end do

        call close_file(nml_unit)
    end if

    call write_version_number(version, tagname)
    log_unit = stdlog()

    if (mpp_pe() == mpp_root_pe()) then
        write (log_unit, nml=cg_drag_nml)
    end if

    i_max = size(lon_bounds) - 1
    j_max = size(lat_bounds) - 1
    q_max = size(p_ref(:)) - 1

    allocate(rays(i_max, j_max, n_max))
    allocate(drays_dt(i_max, j_max, n_max))
    allocate(increments(i_max, j_max, n_max))
    allocate(launches(i_max, j_max, n_source))

    allocate(flux_x(i_max, j_max, q_max + 1))
    allocate(flux_y(i_max, j_max, q_max + 1))
    allocate(z_faces(i_max, j_max, q_max + 1))
    allocate(z_padded(i_max, j_max, q_max + 2))

    rays(:, :, :)%meta = -1
    allocate(last_meta(i_max, j_max))
    last_meta = 0

    do q = 1, q_max
        if (p_ref(q) > source_pressure) then
            q_source = q
            exit
        end if
    end do

    allocate(coriolis(j_max))

    do j = 1, j_max
        lat = 0.5 * (lat_bounds(j) + lat_bounds(j + 1))
        coriolis(j) = 2 * PI * sin(lat) / 86400.
    end do

    ! When H_rho is computed correctly, this will have to be rewritten.
    hgamma = (1. / 2. - 2. / 7.) / H_rho

    call init_source
    call init_nc_output(axes, Time)

end subroutine cg_drag_init

subroutine cg_drag_calc(i_start, j_start, lat, &
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

    flux_x = 0.
    flux_y = 0.

    call take_RK3_step(z_full, uuu, vvv, dt, rays, drays_dt, increments)
    call check_boundaries(z_full, rays)
    call check_source(z_full, uuu, vvv, dt, rays)

    call pad_grid(z_full, z_faces)
    call pad_grid(z_faces, z_padded)

    call update_fluxes(z_padded, rays, flux_x, flux_y)
    call calc_accelerations(z_faces, p_full, temp, flux_x, flux_y, du_dt, dv_dt)
    call send_nc_output(i_start, j_start, Time, flux_x, flux_y, du_dt, dv_dt)

end subroutine cg_drag_calc

subroutine cg_drag_end

    is_initialized = .false.

end subroutine cg_drag_end

! ==============================================================================
! netCDF output
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

    ! --------------------------------------------------------------------------

    if (id_flux_x > 0) then
        i_err = send_data(id_flux_x, flux_x(:, :, 2:), Time, i_start, j_start)
    end if

    if (id_flux_y > 0) then
        i_err = send_data(id_flux_x, flux_y(:, :, 2:), Time, i_start, j_start)
    end if

    if (id_accel_x > 0) then
        i_err = send_data(id_accel_x, du_dt, Time, i_start, j_start)
    end if

    if (id_accel_y > 0) then
        i_err = send_data(id_accel_x, dv_dt, Time, i_start, j_start)
    end if

end subroutine send_nc_output

! ==============================================================================
! type(t_ray) and type(t_tend) helpers
! ==============================================================================

elemental function add_inc_inc(a, b) result(out)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    type(t_inc), intent(in) :: a, b
    type(t_inc)             :: out

    ! --------------------------------------------------------------------------

    out%r = a%r + b%r
    out%m = a%m + b%m

end function add_inc_inc

elemental function add_ray_inc(ray, inc) result(out)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    type(t_ray), intent(in) :: ray
    type(t_inc), intent(in) :: inc
    type(t_ray)             :: out

    ! --------------------------------------------------------------------------

    out = ray
    out%r = ray%r + inc%r
    out%m = ray%m + inc%m

end function add_ray_inc

subroutine delete_at(i, j, n, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                         intent(in)    :: i, j, n
    type(t_ray), dimension(:, :, :), intent(inout) :: rays

    ! --------------------------------------------------------------------------

    rays(i, j, n)%r = 0
    rays(i, j, n)%dr = 0

    rays(i, j, n)%k = 0
    rays(i, j, n)%l = 0
    rays(i, j, n)%m = 0

    rays(i, j, n)%dk = 0
    rays(i, j, n)%dl = 0
    rays(i, j, n)%dm = 0
    rays(i, j, n)%dens = 0

    rays(i, j, n)%age = 0
    rays(i, j, n)%meta = -1

end subroutine delete_at

elemental function mult_scalar_inc(c, inc) result(out)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real,        intent(in) :: c
    type(t_inc), intent(in) :: inc
    type(t_inc)             :: out

    ! --------------------------------------------------------------------------

    out%r = c * inc%r
    out%m = c * inc%m

end function mult_scalar_inc

elemental function mult_scalar_tend(c, tend) result(inc)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real,         intent(in) :: c
    type(t_tend), intent(in) :: tend
    type(t_inc)              :: inc

    ! --------------------------------------------------------------------------

    inc%r = c * tend%r
    inc%m = c * tend%m

end function mult_scalar_tend

! ==============================================================================
! dispersion relation
! ==============================================================================

function cg_r_from_ray(ray, f) result(cg_r)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    type(t_ray), intent(in) :: ray
    real,        intent(in) :: f
    real                    :: cg_r

    ! --------------------------------------------------------------------------
    
    cg_r = cg_r_from_reals(ray%k, ray%l, ray%m, f)

end function cg_r_from_ray

function cg_r_from_reals(k, l, m, f) result(cg_r)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: k, l, m, f
    real             :: cg_r

    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: omega_hat, wvn_sq

    ! --------------------------------------------------------------------------

    wvn_sq = k ** 2 + l ** 2 + m ** 2 + hgamma ** 2
    omega_hat = get_omega_hat(k, l, m, f)

    cg_r = -m * (omega_hat ** 2 - f ** 2) / omega_hat / wvn_sq

end function cg_r_from_reals

function get_dm(m) result(dm)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: m
    real             :: dm

    ! --------------------------------------------------------------------------

    dm = dc_source * (m ** 2) / N0

end function get_dm

function get_m(k, l, f) result(m)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: k, l, f
    real             :: m

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: omega_hat_sq

    ! --------------------------------------------------------------------------

    omega_hat_sq = omega_hat_source ** 2
    m = -sqrt(&
        (k ** 2 + l ** 2) * (N0 ** 2 - omega_hat_sq) / &
        (omega_hat_sq - f ** 2) &
    )

end function get_m

function omega_hat_from_ray(ray, f) result(omega_hat)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    type(t_ray), intent(in) :: ray
    real,        intent(in) :: f
    real                    :: omega_hat

    ! --------------------------------------------------------------------------

    omega_hat = omega_hat_from_reals(ray%k, ray%l, ray%m, f)

end function omega_hat_from_ray

function omega_hat_from_reals(k, l, m, f) result(omega_hat)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: k, l, m, f
    real             :: omega_hat

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: m2

    ! --------------------------------------------------------------------------

    m2 = m ** 2 + hgamma ** 2
    omega_hat = sqrt( &
        (N0 ** 2 * (k ** 2 + l ** 2) + f ** 2 * m2) / &
        (k ** 2 + l ** 2 + m2) &
    )

end function omega_hat_from_reals

! ==============================================================================
! wave source
! ==============================================================================

subroutine check_source(z_full, uuu, vvv, dt, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :),        intent(in)    :: z_full, uuu, vvv
    real,                            intent(in)    :: dt
    type(t_ray), dimension(:, :, :), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer, dimension(size(uuu, 1), size(uuu, 2)) :: n_active, n_excess
    integer :: add_at, i, j, n, s
    ! --------------------------------------------------------------------------

    call update_launches(z_full, uuu, vvv, dt, last_meta, launches)

    n_active = 0
    n_excess = 0

    do n = 1, n_max
        do j = 1, j_max
            do i = 1, i_max
                if (rays(i, j, n)%meta /= -1) then
                    n_active(i, j) = n_active(i, j) + 1
                end if
            end do
        end do
    end do

    do n = 1, n_source
        do j = 1, j_max
            do i = 1, i_max
                if (launches(i, j, n)%meta /= -1) then
                    n_excess(i, j) = n_excess(i, j) + 1
                end if
            end do
        end do
    end do

    n_excess = max(n_excess + n_active - n_max, 0)
    call prune(n_excess, rays)

    do j = 1, j_max
        do i = 1, i_max
            add_at = 1
            do n = 1, n_source
                if (launches(i, j, n)%meta == -1) then
                    cycle
                end if

                do s = add_at, n_max
                    if (rays(i, j, s)%meta == -1) then
                        add_at = s
                        exit
                    end if
                end do

                if (rays(i, j, add_at)%meta /= -1) then
                    call error_mesg("cg_drag_mod", "too many rays", FATAL)
                end if

                rays(i, j, add_at) = launches(i, j, n)
            end do
        end do
    end do
end subroutine check_source

subroutine find_lowest(values, n_find, idx)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :),    intent(in)  :: values
    integer, dimension(:, :),    intent(in)  :: n_find
    integer, dimension(:, :, :), intent(out) :: idx

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n
    integer :: add_at, s
    real, dimension(size(values, 1), size(values, 2), size(idx, 3)) :: lowest

    ! --------------------------------------------------------------------------

    lowest = maxval(values)

    do n = 1, n_max
        do j = 1, j_max
            do i = 1, i_max

                add_at = -1
                do s = 1, n_find(i, j)
                    if (values(i, j, n) < lowest(i, j, s)) then
                        add_at = s
                    else
                        exit
                    end if
                end do

                if (add_at == -1) then
                    cycle
                end if

                do s = 1, add_at - 1
                    lowest(i, j, s) = lowest(i, j, s + 1)
                end do

                lowest(i, j, add_at) = values(i, j, n)
                idx(i, j, add_at) = n
            end do
        end do
    end do

end subroutine find_lowest

subroutine init_source

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: n
    real :: arg, phi, total

    ! --------------------------------------------------------------------------

    n_per_dir = n_source / 4
    dc_source = cp_max / n_per_dir

    allocate(cp_source(n_per_dir))
    allocate(flux_source(n_per_dir))

    do n = 1, n_per_dir
        cp_source(n) = (n - 0.5) * dc_source
        arg = (cp_source(n) - cp_center) / cp_width
        flux_source(n) = exp(-0.5 * (arg ** 2))
    end do

    total = sum(flux_source)
    flux_source = flux_source * (boundary_flux / 4) / total
    omega_hat_source = 2 * PI / T_hat_source

end subroutine init_source

subroutine prune(n_excess, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer, dimension(:, :),        intent(in)    :: n_excess
    type(t_ray), dimension(:, :, :), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: omega_hat, volume
    real, dimension(size(rays, 1), size(rays, 2), size(rays, 3)) :: energy
    integer, dimension(:, :, :), allocatable :: idx
    integer :: i, j, n

    ! --------------------------------------------------------------------------

    do n = 1, n_max
        do j = 1, j_max
            do i = 1, i_max
                omega_hat = get_omega_hat(rays(i, j, n), coriolis(j))
                volume = rays(i, j, n)%dk * rays(i, j, n)%dl * rays(i, j, n)%dm
                energy(i, j, n) = rays(i, j, n)%dens * volume * omega_hat
            end do
        end do
    end do

    allocate(idx(i_max, j_max, maxval(n_excess)))
    call find_lowest(energy, n_excess, idx)

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_excess(i, j)
                call delete_at(i, j, idx(i, j, n), rays)
            end do
        end do
    end do

end subroutine prune

subroutine update_launches(z_full, uuu, vvv, dt, last_meta, launches)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :),        intent(in)    :: z_full, uuu, vvv
    real,                            intent(in)    :: dt
    integer, dimension(:, :),        intent(inout) :: last_meta
    type(t_ray), dimension(:, :, :), intent(out)   :: launches
    
    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n, dir, s
    real :: k, l, m, dm
    real :: cg, cp, prob, u, v, volume, wvn_hor, z_source
    real, dimension(size(launches, 1), size(launches, 2), size(launches, 3)) &
        :: rand
    ! --------------------------------------------------------------------------
    
    call random_number(rand)

    do n = 1, n_per_dir
        do dir = 1, 4
            do j = 1, j_max
                do i = 1, i_max

                    cp  = cp_source(n)
                    if (extrinsic) then
                        u = uuu(i, j, q_max)
                        v = vvv(i, j, q_max)
                        cp = cp - cos_phi(dir) * u - sin_phi(dir) * v
                    end if

                    wvn_hor = omega_hat_source / cp
                    k = wvn_hor * cos_phi(dir)
                    l = wvn_hor * sin_phi(dir)

                    m = get_m(k, l, coriolis(j))
                    cg = get_cg_r(k, l, m, coriolis(j))

                    z_source = z_full(i, j, q_source) 
                    prob = epsilon * cg * dt / dr_source
                    s = (dir - 1) * n_per_dir + n

                    if (rand(i, j, s) < prob) then
                        launches(i, j, s)%r = z_source - 0.5 * dr_source
                        launches(i, j, s)%dr = dr_source

                        launches(i, j, s)%k = k
                        launches(i, j, s)%l = l
                        launches(i, j, s)%m = m

                        launches(i, j, s)%dk = dk_source
                        launches(i, j, s)%dl = dl_source
                        launches(i, j, s)%dm = get_dm(m)

                        volume = dk_source * dl_source * launches(i, j, s)%dm
                        launches(i, j, s)%dens = flux_source(n) / abs( &
                            wvn_hor * volume * cg &
                        )

                        launches(i, j, s)%age = 0
                        launches(i, j, s)%meta = last_meta(i, j)
                        last_meta(i, j) = last_meta(i, j) + 1

                    else
                        launches(i, j, s)%meta = -1
                    end if
                end do
            end do
        end do
    end do


end subroutine update_launches

! ==============================================================================
! time stepping
! ==============================================================================

subroutine check_boundaries(z_full, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :),        intent(in)    :: z_full
    type(t_ray), dimension(:, :, :), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n
    type(t_ray) :: ray
    logical :: delete
    real :: cg, flux, volume, wvn

    ! --------------------------------------------------------------------------

    do n = 1, n_max
        do j = 1, j_max
            do i = 1, i_max
                if (rays(i, j, n)%meta == -1) then
                    cycle
                end if

                ray = rays(i, j, n)
                cg = get_cg_r(ray, coriolis(j))

                volume = ray%dk * ray%dl * ray%dm
                wvn = sqrt(ray%k ** 2 + ray%l ** 2)
                flux = wvn * ray%dens * volume * cg

                delete = ray%r - 0.5 * ray%dr > z_full(i, j, 1)
                delete = delete .or. (abs(flux) < min_flux)
                delete = delete .or. ray%age > max_age

                if (delete) then
                    call delete_at(i, j, n, rays)
                end if
            end do
        end do
    end do

end subroutine check_boundaries

subroutine take_RK3_step(z_full, uuu, vvv, dt, rays, drays_dt, increments)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :),         intent(in)    :: z_full, uuu, vvv
    real,                             intent(in)    :: dt
    type(t_ray), dimension(:, :, :),  intent(inout) :: rays
    type(t_tend), dimension(:, :, :), intent(out)   :: drays_dt
    type(t_inc), dimension(:, :, :),  intent(out)   :: increments

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: stage

    ! --------------------------------------------------------------------------

    call zero_increments(increments)

    do stage = 1, 3
        call update_tendencies(z_full, uuu, vvv, rays, drays_dt)
        increments = As(stage) * increments + dt * drays_dt
        rays = rays + Bs(stage) * increments
    end do

    rays(:, :, :)%age = rays(:, :, :)%age + dt

end subroutine take_RK3_step

subroutine update_tendencies(z_full, uuu, vvv, rays, drays_dt)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :),         intent(in)  :: z_full, uuu, vvv
    type(t_ray), dimension(:, :, :),  intent(in)  :: rays
    type(t_tend), dimension(:, :, :), intent(out) :: drays_dt

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    type(t_ray) :: ray
    real :: du_dr, dv_dr
    real :: z_hi, z_lo, dz
    integer :: i, j, n, q

    ! --------------------------------------------------------------------------

    do n = 1, n_max
        do j = 1, j_max
            do i = 1, i_max
                if (rays(i, j, n)%meta == -1) then
                    drays_dt(i, j, n)%r = 0
                    drays_dt(i, j, n)%m = 0
                    cycle
                end if

                ray = rays(i, j, n)
                do q = 1, q_max - 1
                    z_hi = z_full(i, j, q)
                    z_lo = z_full(i, j, q + 1)
                    dz = z_hi - z_lo

                    if ((z_lo < ray%r) .and. (ray%r < z_hi)) then
                        du_dr = (uuu(i, j, q) - uuu(i, j, q + 1)) / dz
                        dv_dr = (vvv(i, j, q) - vvv(i, j, q + 1)) / dz
                        exit
                    end if
                end do

                drays_dt(i, j, n)%r = get_cg_r(ray, coriolis(j))
                drays_dt(i, j, n)%m = -(ray%k * du_dr + ray%l * dv_dr)
            end do
        end do
    end do

end subroutine update_tendencies

subroutine zero_increments(increments)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    type(t_inc), dimension(:, :, :), intent(out) :: increments

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n

    ! --------------------------------------------------------------------------

    do n = 1, n_max
        do j = 1, j_max
            do i = 1, i_max
                increments(i, j, n)%r = 0
                increments(i, j, n)%m = 0
            end do
        end do
    end do

end subroutine zero_increments

! ==============================================================================
! flux calculations
! ==============================================================================

subroutine calc_accelerations(z_faces, p_full, temp, flux_x, flux_y, &
    du_dt, dv_dt)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :), intent(in)  :: z_faces, temp, p_full, &
        flux_x, flux_y
    real, dimension(:, :, :), intent(out) :: du_dt, dv_dt

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, q
    real :: dz, dFx_dz, dFy_dz, rho
    ! --------------------------------------------------------------------------

    do q = 1, q_max
        do j = 1, j_max
            do i = 1, i_max
                dz = z_faces(i, j, q) - z_faces(i, j, q + 1)
                dFx_dz = (flux_x(i, j, q) - flux_x(i, j, q + 1)) / dz
                dFy_dz = (flux_y(i, j, q) - flux_y(i, j, q + 1)) / dz

                rho = p_full(i, j, q) / temp(i, j, q) / RDGAS
                du_dt(i, j, q) = -dFx_dz / rho
                dv_dt(i, j, q) = -dFy_dz / rho
            end do
        end do
    end do

end subroutine calc_accelerations

subroutine project(z_faces, values, rays, output)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :),        intent(in)  :: z_faces
    real, dimension(:, :, :),        intent(in)  :: values
    type(t_ray), dimension(:, :, :), intent(in)  :: rays
    real, dimension(:, :, :),        intent(out) :: output

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n, q
    real :: r_lo, r_hi, frac
    real :: z_lo, z_hi, dz

    ! --------------------------------------------------------------------------

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_max
                r_lo = rays(i, j, n)%r - 0.5 * rays(i, j, n)%dr
                r_hi = rays(i, j, n)%r + 0.5 * rays(i, j, n)%dr

                do q = 1, q_max
                    z_hi = z_faces(i, j, q)
                    z_lo = z_faces(i, j, q + 1)
                    dz = z_hi - z_lo

                    if (r_hi < z_lo) then
                        cycle
                    end if

                    if (r_lo > z_hi) then
                        exit
                    end if

                    frac = (min(r_hi, z_hi) - max(r_lo, z_lo)) / dz
                    output(i, j, q) = output(i, j, q) + frac * values(i, j, n)
                end do
            end do
        end do
    end do

end subroutine project

subroutine pad_grid(z_in, z_out)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :), intent(in)  :: z_in
    real, dimension(:, :, :), intent(out) :: z_out

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, q, n_q

    ! --------------------------------------------------------------------------

    n_q = size(z_in, 3)

    do q = 2, n_q
        do j = 1, j_max
            do i = 1, i_max
                z_out(i, j, q) = 0.5 * (z_in(i, j, q - 1) + z_in(i, j, q))
            end do
        end do
    end do

    z_out(:, :, 1) = z_in(:, :, 1) + padding_z
    z_out(:, :, n_q + 1) = z_in(:, :, n_q) - padding_z

end subroutine pad_grid

subroutine update_fluxes(z_padded, rays, flux_x, flux_y)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :),        intent(in)  :: z_padded
    type(t_ray), dimension(:, :, :), intent(in)  :: rays
    real, dimension(:, :, :),        intent(out) :: flux_x, flux_y

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real, dimension(size(rays, 1), size(rays, 2), size(rays, 3)) :: action_flux
    integer :: i, j, n, q
    real :: cg, volume

    ! --------------------------------------------------------------------------

    do n = 1, n_max
        do j = 1, j_max
            do i = 1, i_max
                if (rays(i, j, n)%meta == -1) then
                    cycle
                end if

                cg = get_cg_r(rays(i, j, n), coriolis(j))
                volume = rays(i, j, n)%dk * rays(i, j, n)%dl * rays(i, j, n)%dm
                action_flux(i, j, n) = rays(i, j, n)%dens * volume * cg
            end do
        end do
    end do

    call project(z_padded, rays(:, :, :)%k * action_flux, rays, flux_x)
    call project(z_padded, rays(:, :, :)%l * action_flux, rays, flux_y)
    
end subroutine update_fluxes

end module cg_drag_mod