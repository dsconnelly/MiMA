module msgwam_source_mod

! ==============================================================================
! Implements subroutines to initialize and update ray volume sources, and add
! new ray volumes to the system.
! ==============================================================================

use constants_mod,        only: PI
use fms_mod,              only: error_mesg, FATAL, mpp_pe

use msgwam_constants_mod, only: boundary_flux_ex, boundary_flux_tr, cp_max, &
                                cp_width_ex, cp_width_tr, dr_ghost, dr_source, &
                                epsilon, equal_flux, f2, i_max, j_max, &
                                lat_tropics, n_max, n_source, print_prune_diag, &
                                prune_mode, q_max, source_dlat, T_hat_source
use msgwam_rays_mod,      only: delete_ray, get_cg_r, get_dm, get_m, t_ray
use msgwam_utils_mod,     only: get_interp_coeffs, locate

implicit none
private

public LONG_KIND, check_source, init_source, r_source, &
    t_prune_diag

logical :: is_stochastic
integer :: n_launches, n_per_dir
real :: omega_hat_source

real, dimension(:), allocatable :: r_source
real, dimension(:, :), allocatable :: cp_arg, dc_source, flux
logical, dimension(:), allocatable :: is_extrinsic

real, dimension(4) :: COS_PHI = (/ 1., 0., -1., 0. /)
real, dimension(4) :: SIN_PHI = (/ 0., 1., 0., -1. /)

integer, parameter :: LONG_KIND = selected_int_kind(12)
integer, parameter :: MAX_LAUNCHES = 4

type :: t_prune_diag
    integer :: n_pruned
    integer (kind=LONG_KIND) :: total_age
    real :: total_flux, total_r
end type t_prune_diag

contains

subroutine check_source(j, z, u, v, N2, G2, dt, rays, ghosts, last_meta, diag)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                       intent(in)    :: j
    real, dimension(q_max),        intent(in)    :: z, u, v, N2, G2
    real,                          intent(in)    :: dt
    type(t_ray), dimension(n_max), intent(inout) :: rays
    integer, dimension(n_source),  intent(inout) :: ghosts
    integer,                       intent(inout) :: last_meta
    type(t_prune_diag),            intent(out)   :: diag

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
    call prune(j, sum(flux(:, j)), n_excess, rays, diag)

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
    integer :: j, n, p
    real :: arg, closest, cp_width, dist, lat, total
    real, dimension(64) :: lat_in, r_source_in

    ! --------------------------------------------------------------------------

    n_per_dir = n_source / 4
    omega_hat_source = 2 * PI / T_hat_source
    is_stochastic = epsilon < 1.

    allocate(is_extrinsic(j_max))
    allocate(r_source(j_max))

    allocate(cp_arg(n_per_dir, j_max))
    allocate(dc_source(n_per_dir, j_max))
    allocate(flux(n_per_dir, j_max))

    n_launches = n_source * max(1, int(dr_ghost / dr_source) + 1)
    call read_source_level_file(lat_in, r_source_in)

    do j = 1, j_max
        lat = 90 * (lat_bounds(j) + lat_bounds(j + 1)) / PI
        arg = (abs(lat) - (lat_tropics - source_dlat)) / (2 * source_dlat)
        arg = min(max(arg, 0.), 1.) 

        is_extrinsic(j) = abs(lat) > lat_tropics
        cp_width = cp_width_tr * (1 - arg) + cp_width_ex * arg
        total = boundary_flux_tr * (1 - arg) + boundary_flux_ex * arg

        call get_cp_and_flux(cp_width, 0.5 * total, &
            cp_arg(:, j), dc_source(:, j), flux(:, j))

        closest = huge(closest)
        do p = 1, 64
            dist = abs(lat_in(p) - lat)
            if (dist < closest) then
                closest = dist
                r_source(j) = r_source_in(p)
            end if
        end do
    end do

end subroutine init_source

subroutine get_cp_and_flux(cp_width, total, cp_arg, dc_source, flux)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real,                       intent(in)  :: cp_width, total
    real, dimension(n_per_dir), intent(out) :: cp_arg, dc_source, flux

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------    
    integer :: b, n
    real :: accum, dc, flux_per
    real, dimension(n_per_dir + 1) :: edges
    real, dimension(1000001) :: grid, flux_fine

    ! --------------------------------------------------------------------------

    if (.not. equal_flux) then
        dc = cp_max / n_per_dir
        dc_source = dc

        do n = 1, n_per_dir
            cp_arg(n) = (n - 0.5) * dc
            flux(n) = exp(-0.5 * ((cp_arg(n) / cp_width) ** 2))
        end do

    else
        dc = cp_max / (size(grid) - 1)
        flux = 1.

        do n = 1, size(grid)
            grid(n) = (n - 1) * dc
            flux_fine(n) = exp(-0.5 * ((grid(n) / cp_width) ** 2))
        end do

        flux_per = sum(flux_fine) / n_per_dir
        accum = 0.
        
        edges(1) = 0.
        b = 2

        do n = 1, size(grid)
            accum = accum + flux_fine(n)

            if (accum >= flux_per) then
                edges(b) = grid(n)
                accum = 0.
                b = b + 1
            end if

            if (b > n_per_dir + 1) then
                exit
            end if
        end do

        if (b <= n_per_dir + 1) then
            edges(b:) = cp_max
        end if

        dc_source = edges(2:) - edges(:n_per_dir)
        cp_arg = (edges(:n_per_dir) + edges(2:)) / 2.
    end if

    flux = flux * total / sum(flux)

end subroutine get_cp_and_flux

subroutine read_source_level_file(lat_in, r_source_in)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(64), intent(out) :: lat_in, r_source_in

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: iostat, unit

    ! --------------------------------------------------------------------------

    open(newunit=unit, file="INPUT/source-levels.in", &
        iostat=iostat, action="read")
    
    if (iostat /= 0) then
        call error_mesg("msgwam_source_mod", &
            "error loading launch levels", FATAL)
    end if

    read(unit, *) lat_in
    read(unit, *) r_source_in

end subroutine

pure subroutine find_weakest_rays(j, n_find, rays, idx)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                       intent(in)  :: j
    integer,                       intent(in)  :: n_find
    type(t_ray), dimension(n_max), intent(in)  :: rays
    integer, dimension(n_find),    intent(out) :: idx

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: add_at, n, s
    real :: action, criterion, dr, r, wvn
    real, dimension(n_find) :: lowest

    integer, dimension(n_max) :: wdx
    real, dimension(n_max) :: precalc
    real, dimension(4, q_max) :: norms

    ! --------------------------------------------------------------------------

    lowest = huge(lowest)
    precalc = 0.
    norms = 0.
    wdx = 0

    if (prune_mode > 1) then
        do n = 1, n_max
            if (rays(n)%meta == -1) then
                cycle
            end if

            wvn = rays(n)%k + rays(n)%l
            dr = rays(n)%r_hi - rays(n)%r_lo
            action = rays(n)%dens * rays(n)%dm
            precalc(n) = abs(wvn * action * rays(n)%cg_r * dr)

            if (prune_mode == 3) then

                if (rays(n)%k > 0) then
                    wdx(n) = 1
                else if (rays(n)%k < 0) then
                    wdx(n) = 2
                else if (rays(n)%l > 0) then
                    wdx(n) = 3
                else
                    wdx(n) = 4
                end if

                norms(wdx(n), rays(n)%q_mid) = norms(wdx(n), rays(n)%q_mid) &
                    + precalc(n)

            end if
        end do    
    end if

    do n = 1, n_max
        if (rays(n)%meta == -1) then
            cycle
        end if

        if (prune_mode == 1) then
            ! prune by energy
            dr = rays(n)%r_hi - rays(n)%r_lo
            criterion = rays(n)%dens * rays(n)%dm * dr
            criterion = criterion * rays(n)%omega_hat

        else if (prune_mode == 2) then
            ! prune by momentum flux
            criterion = precalc(n)

        else if (prune_mode == 3) then
            ! prune by "importance"
            criterion = precalc(n) / norms(wdx(n), rays(n)%q_mid)

        end if

        r = 0.5 * (rays(n)%r_lo + rays(n)%r_hi)
        if ((r < r_source(j)) .and. (rays(n)%m < 0)) then
            ! exclude rays below the bottom domain unless all rays are such
            criterion = (0.5 * huge(criterion)) + criterion
        end if

        add_at = -1
        do s = 1, n_find
            if (criterion < lowest(s)) then
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

        lowest(add_at) = criterion
        idx(add_at) = n
    end do

end subroutine find_weakest_rays

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
    integer :: add_at, counter, dir, n, n_repeats, q_source, s
    real :: a, b, cg_r, cp_hat, dens, dm, G2_source, k, l, m, omega_hat, &
        N2_source, r, r_ghost, r_hi, r_lo, rate, u_source, v_source, wvn_hor

    real, dimension(MAX_LAUNCHES, n_source) :: rand
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

            if (.not. is_stochastic) then
                if (ghosts(s) > 0) then
                    if (rays(ghosts(s))%r_lo < r_ghost) then
                        cycle
                    end if
                end if
            end if

            cp_hat = cp_arg(n, j) * (COS_PHI(dir) + SIN_PHI(dir))
            if (is_extrinsic(j)) then
                cp_hat = cp_hat &
                    - u_source * abs(COS_PHI(dir)) &
                    - v_source * abs(SIN_PHI(dir))

            end if

            wvn_hor = omega_hat_source / cp_hat
            k = wvn_hor * abs(COS_PHI(dir))
            l = wvn_hor * abs(SIN_PHI(dir))

            r_hi = r_ghost
            if ((.not. is_stochastic)) then
                if (ghosts(s) > 0) then
                    r_hi = rays(ghosts(s))%r_lo
                    r_hi = min(r_hi, r_source(j))
                end if
            end if

            m = get_m(k, l, omega_hat_source ** 2, N2_source, f2(j))
            call get_cg_r(m, k ** 2 + l ** 2, m ** 2 + G2_source, N2_source, &
                f2(j), omega_hat, cg_r)

            if (is_stochastic) then
                rate = epsilon * cg_r * dt / dr_source
                n_repeats = sample_poisson(rand(:, s), rate)
                counter = 0

                if (n_repeats < 1) then
                    cycle
                end if
            end if

            dm = get_dm(m, dc_source(n, j), N2_source)
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

                if (is_stochastic) then
                    counter = counter + 1
                    if (counter == n_repeats) then
                        exit
                    end if
                else if (r_lo <= r_ghost) then
                    exit
                end if

                r_hi = r_lo
            end do

        end do
    end do

end subroutine get_launches

pure function sample_poisson(rand, rate) result(k)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, dimension(MAX_LAUNCHES), intent(in) :: rand
    real,                          intent(in) :: rate
    integer                                   :: k

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i
    real :: L, prod

    ! --------------------------------------------------------------------------

    L = exp(-rate)
    prod = 1.
    k = 0

    do i = 1, MAX_LAUNCHES
        prod = prod * rand(i)
        k = k + 1

        if (prod < L) then
            exit
        end if
    end do

    k = k - 1

end function sample_poisson

pure subroutine prune(j, norm, n_excess, rays, diag)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                       intent(in)    :: j
    real,                          intent(in)    :: norm
    integer,                       intent(in)    :: n_excess
    type(t_ray), dimension(n_max), intent(inout) :: rays
    type(t_prune_diag),            intent(out)   :: diag

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: n
    integer, dimension(n_excess) :: idx
    real :: flux, r

    ! --------------------------------------------------------------------------

    diag%n_pruned = 0
    diag%total_age = 0
    diag%total_flux = 0.
    diag%total_r = 0.

    if (n_excess > 0) then
        call find_weakest_rays(j, n_excess, rays, idx)

        if (print_prune_diag) then
            do n = 1, n_excess
                diag%n_pruned = diag%n_pruned + 1
                diag%total_age = diag%total_age + rays(idx(n))%age

                r = 0.5 * (rays(idx(n))%r_lo + rays(idx(n))%r_hi)
                flux = sqrt(rays(idx(n))%k ** 2 + rays(idx(n))%l ** 2)
                flux = flux * rays(idx(n))%dens * rays(idx(n))%dm
                flux = flux * rays(idx(n))%cg_r

                diag%total_flux = diag%total_flux + abs(flux / norm)
                diag%total_r = diag%total_r + r / 1000.
            end do
        end if

        do n = 1, n_excess
            !DIR$ NOINLINE
            call delete_ray(rays(idx(n)))
        end do
    end if

end subroutine prune

end module msgwam_source_mod