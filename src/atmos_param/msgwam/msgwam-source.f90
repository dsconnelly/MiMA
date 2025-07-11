module msgwam_source_mod

! ==============================================================================
! Implements subroutines to initialize and update ray volume sources, and add
! new ray volumes to the system.
! ==============================================================================

use constants_mod,        only: PI
use fms_mod,              only: error_mesg, FATAL, mpp_pe

use msgwam_constants_mod, only: boundary_flux_ex, boundary_flux_tr, cp_max, &
                                cp_width_ex, cp_width_tr, dr_ghost, dr_source, &
                                epsilon, f2, j_max, lat_tropics, n_max, &
                                n_source, print_prune_diag, q_max, &
                                r_source_ex, r_source_tr, source_dlat, &
                                steady_state, T_hat_source
use msgwam_rays_mod,      only: delete_ray, get_cg_r, get_dm, get_m, t_ray
use msgwam_utils_mod,     only: get_interp_coeffs, locate

implicit none
private

public check_source, init_source, r_source

logical :: is_stochastic
integer :: n_launches, n_per_dir
real :: dc_source, omega_hat_source

logical, dimension(:), allocatable :: is_extrinsic
real, dimension(:), allocatable :: cp_arg, cp_width, r_source
real, dimension(:, :), allocatable :: flux

real, dimension(4) :: COS_PHI = (/ 1., 0., -1., 0. /)
real, dimension(4) :: SIN_PHI = (/ 0., 1., 0., -1. /)

contains

subroutine check_source(j, z, u, v, N2, G2, dt, rays, ghosts, last_meta)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                       intent(in)    :: j
    real, dimension(q_max),        intent(in)    :: z, u, v, N2, G2
    real,                          intent(in)    :: dt
    type(t_ray), dimension(n_max), intent(inout) :: rays
    integer, dimension(n_source),  intent(inout) :: ghosts
    integer,                       intent(inout) :: last_meta

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: add_at, gid, n, n_excess, s
    type(t_ray), dimension(n_launches) :: launches

    ! --------------------------------------------------------------------------

    n_excess = count(rays(:)%meta /= -1)
    call get_launches(j, z, u, v, N2, G2, dt, rays, ghosts, last_meta, &
        n_excess, launches)

    n_excess = max(n_excess - n_max, 0)
    call prune(n_excess, rays)

    add_at = 1
    do n = 1, n_launches
        if (launches(n)%meta == -1) then
            cycle
        end if

        do s = add_at, n_max
            if (rays(s)%meta == -1) then
                add_at = s
                exit
            end if
        end do

        if (rays(add_at)%meta /= -1) then
            call error_mesg("msgwam", "too many rays", FATAL)
        end if

        rays(add_at) = launches(n)
        if (.not. is_stochastic) then
            gid = launches(n)%ghost_id

            if (gid /= -1) then
                if (ghosts(gid) > 0) then
                    rays(ghosts(gid))%ghost_id = -1
                end if

                ghosts(gid) = add_at
            end if
        end if
    end do

end subroutine check_source

subroutine init_source(lat_bounds)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(j_max + 1), intent(in) :: lat_bounds

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: j, n
    real :: arg, lat, total

    ! --------------------------------------------------------------------------

    n_per_dir = n_source / 4
    dc_source = cp_max / n_per_dir
    omega_hat_source = 2 * PI / T_hat_source
    is_stochastic = epsilon > 0.

    allocate(cp_arg(n_per_dir))
    allocate(is_extrinsic(j_max))
    allocate(flux(n_per_dir, j_max))
    allocate(cp_width(j_max))
    allocate(r_source(j_max))

    n_launches = n_source * max(1, int(dr_ghost / dr_source) + 1)

    do n = 1, n_per_dir
        cp_arg(n) = (n - 0.5) * dc_source
    end do

    do j = 1, j_max
        lat = abs(90 * (lat_bounds(j) + lat_bounds(j + 1)) / PI)
        arg = (lat - (lat_tropics - source_dlat)) / (2 * source_dlat)
        arg = min(max(arg, 0.), 1.) 

        is_extrinsic(j) = lat > lat_tropics
        cp_width(j) = cp_width_tr * (1 - arg) + cp_width_ex * arg
        r_source(j) = r_source_tr * (1 - arg) + r_source_ex * arg

        do n = 1, n_per_dir
            flux(n, j) = exp(-0.5 * ((cp_arg(n) / cp_width(j)) ** 2))
        end do

        total = boundary_flux_tr * (1 - arg) + boundary_flux_ex * arg
        flux(:, j) = 0.5 * flux(:, j) * total / sum(flux(:, j))
    end do

end subroutine init_source

pure subroutine find_lowest_energies(n_find, rays, idx)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                       intent(in)  :: n_find
    type(t_ray), dimension(n_max), intent(in)  :: rays
    integer, dimension(n_find),    intent(out) :: idx

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: dr, energy
    integer :: add_at, n, s
    real, dimension(n_find) :: lowest

    ! --------------------------------------------------------------------------

    lowest = huge(lowest)

    do n = 1, n_max
        if ((rays(n)%meta == -1) .or. (rays(n)%ghost_id /= -1)) then
            cycle
        end if

        dr = rays(n)%r_hi - rays(n)%r_lo
        energy = rays(n)%dens * rays(n)%dm * rays(n)%omega_hat * dr

        add_at = -1
        do s = 1, n_find
            if (energy < lowest(s)) then
                add_at = s
            else
                exit
            end if
        end do

        if (add_at == -1) then
            cycle
        end if

        lowest(:add_at - 1) = lowest(2:add_at)
        idx(:add_at - 1) = idx(2:add_at)

        lowest(add_at) = energy
        idx(add_at) = n
    end do

end subroutine find_lowest_energies

subroutine get_launches(j, z, u, v, N2, G2, dt, rays, ghosts, &
    last_meta, n_added, launches)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                            intent(in)    :: j
    real, dimension(q_max),             intent(in)    :: z, u, v, N2, G2
    real,                               intent(in)    :: dt
    type(t_ray), dimension(n_max),      intent(in)    :: rays
    integer, dimension(n_source),       intent(in)    :: ghosts
    integer,                            intent(inout) :: last_meta, n_added
    type(t_ray), dimension(n_launches), intent(out)   :: launches

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: add_at, dir, n, q_source, s
    real :: a, b, cg_r, cp_hat, dens, dm, G2_source, k, l, m, omega_hat, &
        N2_source, prob, r, r_ghost, r_hi, r_lo, u_source, v_source, wvn_hor

    real, dimension(n_source) :: rand
    real, dimension(q_max - 1) :: dz_inv

    ! --------------------------------------------------------------------------

    if (is_stochastic) then
        call random_number(rand)
    end if

    launches(:)%meta = -1
    launches(:)%ghost_id = -1

    r_ghost = r_source(j) - dr_ghost
    q_source = locate(z, r_source(j), 30)
    dz_inv = 1. / (z(:q_max - 1) - z(2:))

    call get_interp_coeffs(z, dz_inv, r_source(j), q_source, a, b)    
    N2_source = a * N2(q_source) + b * N2(q_source + 1)
    G2_source = a * G2(q_source) + b * G2(q_source + 1)

    if (is_extrinsic(j)) then
        u_source = a * u(q_source) + b * u(q_source + 1)
        v_source = a * v(q_source) + b * v(q_source + 1)
    end if

    add_at = 1
    do dir = 1, 4
        do n = 1, n_per_dir
            s = (dir - 1) * n_per_dir + n

            if (.not. (is_stochastic .or. steady_state)) then
                if (ghosts(s) > 0) then
                    if (rays(ghosts(s))%r_lo < r_ghost) then
                        cycle
                    end if
                end if
            end if

            cp_hat = cp_arg(n) * (COS_PHI(dir) + SIN_PHI(dir))
            if (is_extrinsic(j)) then
                cp_hat = cp_hat &
                    - u_source * abs(COS_PHI(dir)) &
                    - v_source * abs(SIN_PHI(dir))

            end if

            wvn_hor = omega_hat_source / cp_hat
            k = wvn_hor * abs(COS_PHI(dir))
            l = wvn_hor * abs(SIN_PHI(dir))

            if (is_stochastic .or. steady_state) then
                r_hi = r_source(j)
            else if (ghosts(s) > 0) then
                r_hi = rays(ghosts(s))%r_lo
                r_hi = min(r_hi, r_source(j))
            else
                r_hi = r_ghost
            end if

            m = get_m(k, l, omega_hat_source ** 2, N2_source, f2(j))
            call get_cg_r(m, k ** 2 + l ** 2, m ** 2 + G2_source, N2_source, &
                f2(j), omega_hat, cg_r)

            if (is_stochastic) then
                prob = epsilon * cg_r * dt / dr_source
                if (rand(s) .ge. prob) then
                    cycle
                end if
            end if

            dm = get_dm(m, dc_source, N2_source)
            dens = flux(n, j) / abs(wvn_hor * dm * cg_r)

            do
                r_lo = r_hi - dr_source
                r = 0.5 * (r_lo + r_hi)

                launches(add_at)%r_hi = r_hi
                launches(add_at)%r_lo = r_lo

                launches(add_at)%k = k
                launches(add_at)%l = l
                launches(add_at)%m = m

                launches(add_at)%dm = dm
                launches(add_at)%dens = dens

                launches(add_at)%cg_r = cg_r
                launches(add_at)%omega_hat = omega_hat
                launches(add_at)%G2 = G2_source

                launches(add_at)%age = 0
                launches(add_at)%meta = last_meta
                launches(add_at)%ghost_id = merge(s, -1, r_lo <= r_ghost)

                launches(add_at)%q_hi = locate(z, r_hi, 30)
                launches(add_at)%q_lo = locate(z, r_lo, 30)
                launches(add_at)%q_mid = locate(z, r, 30)

                last_meta = last_meta + 1
                n_added = n_added + 1
                add_at = add_at + 1

                if (is_stochastic .or. steady_state .or. (r_lo <= r_ghost)) then
                    exit
                end if

                r_hi = r_lo
            end do

        end do
    end do

end subroutine get_launches

pure subroutine prune(n_excess, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                       intent(in)    :: n_excess
    type(t_ray), dimension(n_max), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: n
    integer, dimension(n_excess) :: idx

    ! --------------------------------------------------------------------------

    if (n_excess > 0) then
        call find_lowest_energies(n_excess, rays, idx)
        do n = 1, n_excess
            call delete_ray(rays(idx(n)))
        end do
    end if

end subroutine prune

end module msgwam_source_mod