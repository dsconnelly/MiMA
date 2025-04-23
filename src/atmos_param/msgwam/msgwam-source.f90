module msgwam_source_mod

! ==============================================================================
! Implements subroutines to initialize and update ray volume sources, and add
! new ray volumes to the system.
! ==============================================================================

use constants_mod,        only: PI
use fms_mod,              only: error_mesg, FATAL

use msgwam_constants_mod, only: boundary_flux, cp_max, cp_width, dr_source, &
                                epsilon, f2, i_max, j_max, is_extrinsic, &
                                n_max, n_source, q_max, source_pressure, &
                                T_hat_source
use msgwam_rays_mod,      only: get_cg_r, get_dm, get_m, delete_ray, t_ray
use msgwam_utils_mod,     only: get_interp_coeffs, locate

implicit none
private

public check_source, init_source

logical :: is_stochastic
integer :: n_per_dir, q_source
real :: dc_source, omega_hat_source
real, dimension(:), allocatable :: cp_source
type(t_ray), dimension(:, :, :), allocatable :: launches

real, dimension(4) :: COS_PHI = (/ 1., 0., -1., 0. /)
real, dimension(4) :: SIN_PHI = (/ 0., 1., 0., -1. /)

contains

subroutine check_source(z_centers, u_bar, v_bar, N2, G2, dt, &
    rays, ghosts, last_meta)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:q_max + 1, i_max, j_max),  intent(in)    :: z_centers
    real, dimension(q_max, i_max, j_max),        intent(in)    :: u_bar, &
                                                                  v_bar, N2, G2
    real,                                        intent(in)    :: dt
    type(t_ray), dimension(n_max, i_max, j_max), intent(inout) :: rays
    integer, dimension(n_source, i_max, j_max),  intent(inout) :: ghosts
    integer, dimension(i_max, j_max),            intent(inout) :: last_meta

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: add_at, i, j, n, s
    integer, dimension(i_max, j_max) :: n_excess

    ! --------------------------------------------------------------------------

    n_excess = count(rays(:, :, :)%meta /= -1, dim=1)
    call update_launches(z_centers, u_bar, v_bar, N2, G2, dt, rays, ghosts, &
        last_meta, n_excess, launches)

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

pure subroutine find_lowest_energies(n_find, rays, idx)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer, dimension(i_max, j_max),            intent(in)  :: n_find
    type(t_ray), dimension(n_max, i_max, j_max), intent(in)  :: rays
    integer, dimension(:, :, :),                 intent(out) :: idx

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: energy
    integer :: add_at, i, j, n, s
    real, dimension(size(idx, 1), i_max, j_max) :: lowest

    ! --------------------------------------------------------------------------

    lowest = huge(lowest)

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_max

                associate(ray => rays(n, i, j))
                    if ((ray%meta == -1) .or. ray%is_ghost) then
                        cycle
                    end if

                    energy = ray%dens * ray%dm * ray%omega_hat * &
                        (ray%r_hi - ray%r_lo)
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

subroutine init_source(p_ref)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max) :: p_ref

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: n, q

    ! --------------------------------------------------------------------------

    n_per_dir = n_source / 4
    dc_source = cp_max / n_per_dir
    omega_hat_source = 2 * PI / T_hat_source
    is_stochastic = epsilon > 0.

    allocate(cp_source(n_per_dir))
    allocate(launches(n_source, i_max, j_max))

    do n = 1, n_per_dir
        cp_source(n) = (n - 0.5) * dc_source
    end do

    do q = 1, q_max
        if (p_ref(q) > source_pressure) then
            q_source = q
            exit
        end if
    end do

end subroutine init_source

pure subroutine prune(n_excess, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer, dimension(i_max, j_max),            intent(in)    :: n_excess
    type(t_ray), dimension(n_max, i_max, j_max), intent(inout) :: rays

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
                call delete_ray(rays(idx(n, i, j), i, j))
            end do
        end do
    end do

end subroutine prune

subroutine update_launches(z_centers, u_bar, v_bar, N2, G2, dt, rays, ghosts, &
    last_meta, n_added, launches)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:q_max + 1, i_max, j_max),     intent(in)    :: z_centers
    real, dimension(q_max, i_max, j_max),           intent(in)    :: u_bar, &
                                                                     v_bar, &
                                                                     N2, G2
    real,                                           intent(in)    :: dt
    type(t_ray), dimension(n_max, i_max, j_max),    intent(in)    :: rays
    integer, dimension(n_source, i_max, j_max),     intent(in)    :: ghosts
    integer, dimension(i_max, j_max),               intent(inout) :: last_meta
    integer, dimension(i_max, j_max),               intent(inout) :: n_added
    type(t_ray), dimension(n_source, i_max, j_max), intent(out)   :: launches
    
    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    logical :: cleared
    integer :: dir, i, j, n, q_hi, q_lo, q_mid, s
    real :: a, b, cg, cp, flux, G2_source, k, l, m, mag_cp_hat, mag_wvn_hor, &
            N2_source, prob, r, r_lo, r_hi, total, u, v
    real, dimension(n_source, i_max, j_max) :: rand

    ! --------------------------------------------------------------------------

    call random_number(rand)

    do j = 1, j_max
        do i = 1, i_max
            r_hi = z_centers(q_source, i, j)
            r_lo = r_hi - dr_source
            r = (r_lo + r_hi) / 2.

            associate( &
                z => z_centers(1:q_max, i, j), &
                N2_col => N2(:, i, j), &
                G2_col => G2(:, i, j) &
            )
                q_hi = locate(r_hi, z, q_source)
                q_lo = locate(r_lo, z, q_source)
                q_mid = locate(r, z, q_source)

                call get_interp_coeffs(z, r, q_mid, a, b)
                N2_source = a * N2_col(q_mid) + b * N2_col(q_mid + 1)
                G2_source = a * G2_col(q_mid) + b * G2_col(q_mid + 1)
            end associate

            total = 0.
            do dir = 1, 4
                do n = 1, n_per_dir

                    s = (dir - 1) * n_per_dir + n
                    cleared = .true.

                    if (.not. is_stochastic) then
                        if (ghosts(s, i, j) > 0) then
                            if (rays(ghosts(s, i, j), i, j)%r_lo < r_hi) then
                                launches(s, i, j)%meta = -1
                                cleared = .false.
                            end if
                        end if
                    end if

                    cp = cp_source(n)
                    u = u_bar(q_source, i, j)
                    v = v_bar(q_source, i, j)
                    mag_cp_hat = cp - COS_PHI(dir) * u - SIN_PHI(dir) * v

                    if (is_extrinsic(j)) then
                        flux = exp(-0.5 * ((cp / cp_width) ** 2))
                    else
                        flux = exp(-0.5 * ((mag_cp_hat / cp_width) ** 2))
                    end if

                    total = total + flux
                    if (.not. cleared) then
                        cycle
                    end if

                    mag_wvn_hor = omega_hat_source / mag_cp_hat
                    k = mag_wvn_hor * COS_PHI(dir)
                    l = mag_wvn_hor * SIN_PHI(dir)

                    m = get_m(k, l, omega_hat_source ** 2, N2_source, f2(j))
                    cg = get_cg_r(k, l, m, N2_source, f2(j), G2_source)

                    associate (ray => launches(s, i, j))
                        if (is_stochastic) then
                            prob = epsilon * cg * dt / dr_source
                            if (rand(s, i, j) .ge. prob) then
                                launches(s, i, j)%meta = -1
                                cycle
                            end if
                        end if

                        ray%r_hi = r_hi
                        ray%r_lo = r_lo

                        ray%k = k
                        ray%l = l
                        ray%m = m

                        ray%dm = get_dm(m, dc_source, N2_source)
                        ray%dens = flux / abs(mag_wvn_hor * ray%dm * cg)

                        ray%cg_r = cg
                        ray%omega_hat = omega_hat_source
                        ray%G2 = G2_source

                        ray%age = 0
                        ray%meta = last_meta(i, j)
                        ray%is_ghost = .not. is_stochastic

                        ray%q_hi = q_hi
                        ray%q_lo = q_lo
                        ray%q_mid = q_mid

                        last_meta(i, j) = last_meta(i, j) + 1
                        n_added(i, j) = n_added(i, j) + 1
                    end associate

                end do
            end do

            launches(:, i, j)%dens = launches(:, i, j)%dens &
                * 2 * boundary_flux / total

        end do
    end do

end subroutine update_launches

end module msgwam_source_mod