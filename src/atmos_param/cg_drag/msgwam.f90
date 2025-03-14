module cg_drag_mod

! ------------------------------------------------------------------------------
! This version of cg_drag_mod implements MS-GWaM as described in Boloni et al.
! (2021) and implemented in dsconnelly/python-msgwam on Github.
! ------------------------------------------------------------------------------

use constants_mod,    only: PI, constants_init
use fms_mod,          only: check_nml_error, close_file, file_exist, fms_init, &
                            mpp_pe, mpp_root_pe, open_namelist_file, stdlog
use time_manager_mod, only: time_manager_init, time_type

implicit none
private

character(len=128) :: version = 'cg_drag_msgwam.f90, 2025/03/14'
character(len=128) :: tagname = 'cayuga'

! ------------------------------------------------------------------------------
! interfaces
! ------------------------------------------------------------------------------

public cg_drag_init, cg_drag_calc, cg_drag_end

! ------------------------------------------------------------------------------
! namelist
! ------------------------------------------------------------------------------

real :: bc_flux = 0.01 ! (Pa)
real :: damp_level_pressure = 0.8e+02 ! (Pa)
real :: dr_source = 1000 ! (m)
real :: cp_center = 15 ! (m / s)
real :: cp_max = 50 ! (m / s)
real :: cp_width = 10 ! (m / s)
real :: dk_source = 0.0001 ! (1 / m)
real :: dl_source = 0.0001 ! (1 / m)
real :: epsilon = 1
logical :: extrinsic = .true.
real :: H_rho = 8.e3 ! (m)
real :: N0 = 0.015 ! (1 / s)
integer :: n_max = 2500
integer :: n_source = 48
real :: source_level_pressure = 300.e+02 ! (Pa)
real :: T_hat = 10 ! (hours)

namelist / cg_drag_nml / &
    bc_flux, damp_level_pressure, dr_source, cp_center, cp_max, cp_width, &
    dk_source, dl_source, epsilon, extrinsic, H_rho, N0, n_max, n_source, &
    source_level_pressure, T_hat

! ------------------------------------------------------------------------------
! private variables
! ------------------------------------------------------------------------------

integer, parameter :: N_PROPS = 11
logical :: module_is_initialized = .false.

integer :: k_source, k_damp
real, dimension(:), allocatable :: f_Cor

real :: dc_source
real, dimension(:), allocatable :: cp_source
real, dimension(:), allocatable :: flux_source
real, dimension(:, :), allocatable :: last_meta
integer :: n_per_dir
real :: omega_hat_source
real, dimension(4) :: phi_source

real, dimension(:, :, :, :), allocatable :: rays
integer, dimension(:, :, :), allocatable :: is_active

contains

! ------------------------------------------------------------------------------
! public subroutines
! ------------------------------------------------------------------------------

subroutine cg_drag_init(lonb, latb, pref, Time, axes)

    ! --------------------
    ! intent(in) variables
    ! --------------------

    real,            dimension(:), intent(in) :: lonb, latb, pref
    integer,         dimension(4), intent(in) :: axes
    type(time_type),               intent(in) :: Time

    ! lonb, latb : model longitudes and latitudes at cell corners
    ! pref : reference pressures at full levels, plus surface value
    ! Time : current time
    ! axes : unused

    ! ---------------
    ! local variables
    ! ---------------
    integer :: ierr, io, log_unit, nml_unit
    integer :: i_max, j_max, k_max
    integer :: i, j, k, n

    real :: arg, total
    real :: lat

    if (module_is_initialized) return

    call fms_init
    call time_manager_init
    call constants_init

    if (file_exist('input.nml')) then
        nml_unit = open_namelist_file()
        ierr = 1

        do while (ierr /= 0)
            read(nml_unit, nml=cg_drag_nml, iostat=io)
            ierr = check_nml_error(io, 'cg_drag_nml')
            ! maybe : add 'end" argument to read statement
        end do

        call close_file(nml_unit)
    end if

    call write_version_number(version, tagname)
    log_unit = stdlog()
    if (mpp_pe() == mpp_root_pe()) write (log_unit, nml=cg_drag_nml)

    k_max = size(pref(:)) - 1

    do k = 1, k_max
        if (pref(k) < damp_level_pressure) then
            k_damp = k
        end if

        if (pref(k) > source_level_pressure) then
            k_source = k
            exit
        end if
    end do

    do n = 1, 4
        phi_source(n) = (n - 1) * PI / 2
    end do

    n_per_dir = n_source / 4
    allocate(cp_source(n_per_dir))
    allocate(flux_source(n_per_dir))
    dc_source = cp_max / n_per_dir
    
    do n = 1, n_per_dir
        cp_source(n) = (n - 0.5) * dc_source
        arg = (cp_source(n) - cp_center) / cp_width
        flux_source(n) = exp(-0.5 * (arg ** 2))
    end do

    total = sum(flux_source)
    flux_source = flux_source * (bc_flux / 4) / total
    omega_hat_source = 2 * PI / (T_hat * 3600)

    i_max = size(lonb) - 1
    j_max = size(latb) - 1

    allocate(rays(i_max, j_max, N_PROPS, n_max))
    allocate(is_active(i_max, j_max, n_max))
    allocate(last_meta(i_max, j_max))

    allocate(f_Cor(j_max))
    do n = 1, j_max
        lat = 0.5 * (latb(n) + latb(n + 1))
        f_Cor(n) = 2 * PI * sin(lat) / 86400
    end do

end subroutine cg_drag_init

subroutine cg_drag_calc(is, js, lat, &
    p_full, z_full, temp, uuu, vvv, &
    Time, dt, du_dt, dv_dt)

    ! ---------
    ! arguments
    ! ---------
    integer,                   intent(in) :: is, js
    real, dimension(:, :),     intent(in) :: lat
    real, dimension(:, :, :),  intent(in) :: p_full, z_full, temp, uuu, vvv
    type(time_type),           intent(in) :: Time
    real,                      intent(in) :: dt
    real, dimension(:, :, :), intent(out) :: du_dt, dv_dt

    ! ---------------
    ! local variables
    ! ---------------

    call take_RK3_step(z_full, uuu, vvv, dt)

    du_dt = 1. / 86400.
    dv_dt = 1. / 86400.

end subroutine cg_drag_calc

subroutine cg_drag_end

    module_is_initialized = .false.

end subroutine cg_drag_end

! ------------------------------------------------------------------------------
! private subroutines
! ------------------------------------------------------------------------------

subroutine take_RK3_step(z_full, uuu, vvv, dt)

    ! ---------
    ! arguments
    ! ---------
    real, dimension(:, :, :) :: z_full, uuu, vvv
    real :: dt

    ! ---------------
    ! local variables
    ! ---------------
    integer :: step
    real, dimension(3) :: As, Bs
    real, dimension(:, :, :, :), allocatable :: drays_dt, increment

    As = (/ 0., -5. / 9., -153. / 128. /)
    Bs = (/ 1. / 3., 15. / 16., 8. / 15. /)

    allocate(drays_dt(size(rays, 1), size(rays, 2), N_PROPS - 2, n_max))
    allocate(increment(size(rays, 1), size(rays, 2), N_PROPS - 2, n_max))
    increment = 0

    do step = 1, 3
        call get_drays_dt(z_full, uuu, vvv, drays_dt)
        increment = drays_dt * dt + As(step) * increment
        rays(:, :, 1:9, :) = rays(:, :, 1:9, :) + Bs(step) * increment
    end do

    rays(:, :, 10, :) = rays(:, :, 10, :) + dt

end subroutine

subroutine get_drays_dt(z_full, uuu, vvv, drays_dt)

    ! ---------
    ! arguments
    ! ---------
    real, dimension(:, :, :), intent(in) :: z_full, uuu, vvv
    real, dimension(:, :, :, :), intent(out) :: drays_dt
    ! should have shape (n_i, n_j, N_PROPS - 2, n_max)

    ! ---------------
    ! local variables
    ! ---------------
    integer :: i, j, n, q
    integer :: i_max, j_max, q_max
    real :: r, k, l, m
    real :: du_dr, dv_dr
    real :: z_hi, z_lo, dz

    i_max = size(uuu, 1)
    j_max = size(uuu, 2)
    q_max = size(uuu, 3)

    do n = 1, n_max
        do j = 1, j_max
            do i = 1, i_max
                if (.not. is_active(i, j, n)) then
                    drays_dt(i, j, :, n) = 0
                    cycle
                end if

                r = rays(i, j, 1, n)
                k = rays(i, j, 3, n)
                l = rays(i, j, 4, n)
                m = rays(i, j, 5, n)

                drays_dt(i, j, 1, n) = get_cg_r(k, l, m, f_Cor(j))
                drays_dt(i, j, 2, n) = 0

                drays_dt(i, j, 3, n) = 0
                drays_dt(i, j, 4, n) = 0

                do q = 1, q_max - 1
                    z_hi = z_full(i, j, q)
                    z_lo = z_full(i, j, q + 1)
                    dz = z_hi - z_lo

                    if ((z_lo < r) .and. (r < z_hi)) then
                        du_dr = (uuu(i, j, q) - uuu(i, j, q + 1)) / dz
                        dv_dr = (vvv(i, j, q) - vvv(i, j, q + 1)) / dz

                        drays_dt(i, j, 5, q) = -(k * du_dr + l * dv_dr)
                        exit
                    end if
                end do

                drays_dt(i, j, 6, n) = 0
                drays_dt(i, j, 7, n) = 0
                drays_dt(i, j, 8, n) = 0
                drays_dt(i, j, 9, n) = 0
            end do
        end do
    end do

end subroutine get_drays_dt

function get_launches(lat, z_full, uuu, vvv, dt) result(launches)
    
    ! ---------
    ! arguments
    ! ---------
    real, dimension(:) :: lat
    real, dimension(:, :, :) :: z_full, uuu, vvv
    real :: dt
    
    real, dimension(:, :, :, :), allocatable :: launches

    ! ---------------
    ! local variables
    ! ---------------
    integer :: i, j, n, p, q
    integer :: bottom, i_max, j_max
    real :: cos_phi, sin_phi
    real :: prob
    real, dimension(:, :, :), allocatable :: rand
    real :: cg, cp, volume, wvn_hor
    real :: k, l, m, dm, dens
    real :: u, v
    
    i_max = size(uuu, 1)
    j_max = size(uuu, 2)
    bottom = size(uuu, 3)

    allocate(launches(i_max, j_max, N_PROPS - 1, n_source))
    allocate(rand(i_max, j_max, n_source))
    call random_number(rand)

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_per_dir
                do p = 1, 4

                    cos_phi = cos(phi_source(p))
                    sin_phi = sin(phi_source(p))
                
                    cp = cp_source(n)
                    if (extrinsic) then
                        u = uuu(i, j, bottom)
                        v = vvv(i, j, bottom)
                        cp = cp - cos_phi * u - sin_phi * v
                    end if

                    wvn_hor = omega_hat_source / cp
                    k = wvn_hor * cos_phi
                    l = wvn_hor * sin_phi

                    m = get_m(k, l, f_Cor(j))
                    cg = get_cg_r(k, l, m, f_Cor(j))
                    dm = get_dm(m)

                    volume = dk_source * dl_source * dm
                    dens = flux_source(n) / abs(wvn_hor * volume * cg)

                    prob = epsilon * cg * dt / dr_source
                    q = (p - 1) * n_per_dir + n

                    if (rand(i, j, q) < prob) then
                        launches(i, j, 1, q) = z_full(k_source) - 0.5 * dr_source
                        launches(i, j, 2, q) = dr_source

                        launches(i, j, 3, q) = k
                        launches(i, j, 4, q) = l
                        launches(i, j, 5, q) = m

                        launches(i, j, 6, q) = dk_source
                        launches(i, j, 7, q) = dl_source
                        launches(i, j, 8, q) = dm

                        launches(i, j, 9, q) = dens
                        launches(i, j, 10, q) = 0
                    else
                        launches(i, j, :, q) = -1
                    end if

                end do                
            end do
        end do
    end do

end function get_launches

function get_cg_r(k, l, m, f) result(cg_r)

    ! ---------
    ! arguments
    ! ---------
    real :: k, l, m, f
    real :: cg_r

    ! ---------------
    ! local variables
    ! ---------------
    real :: omega_hat, wvn_sq

    wvn_sq = k ** 2 + l ** 2 + m ** 2 + hgamma() ** 2
    omega_hat = get_omega_hat(k, l, m, f)

    cg_r = -m * (omega_hat ** 2 - f ** 2) / omega_hat / wvn_sq

end function get_cg_r

function get_m(k, l, f) result(m)

    ! ---------
    ! arguments
    ! ---------
    real :: k, l, f, m

    ! ---------------
    ! local variables
    ! ---------------
    real :: omega_hat_sq

    omega_hat_sq = omega_hat_source ** 2
    m = -sqrt( &
        (k ** 2 + l ** 2) * (N0 ** 2 - omega_hat_sq) / &
        (omega_hat_sq - f ** 2) &
    )

end function get_m

function get_dm(m) result(dm)

    ! ---------
    ! arguments
    ! ---------
    real :: m, dm

    dm = dc_source * (m ** 2) / N0

end function get_dm

real function hgamma()

    hgamma = (1. / 2. - 2. / 7.) / H_rho

end function hgamma

function get_omega_hat(k, l, m, f) result(omega_hat)

    ! ---------
    ! arguments
    ! ---------
    real :: k, l, m, f
    real :: omega_hat

    ! ---------------
    ! local variables
    ! ---------------
    real :: m2

    m2 = (m ** 2) + (hgamma() ** 2)

    omega_hat = sqrt( &
        (N0 ** 2 * (k ** 2 + l ** 2) + f ** 2 * m2) / &
        (k ** 2 + l ** 2 + m2) &
    )

end function get_omega_hat

end module cg_drag_mod