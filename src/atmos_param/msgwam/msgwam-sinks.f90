module msgwam_sinks_mod

! ==============================================================================
! This module implements ray volume dissipation and breaking.
! ==============================================================================

use msgwam_constants_mod, only: break_waves, f2, i_max, j_max, max_age, &
                                min_flux, mu, n_max, n_sponge, q_max
use msgwam_rays_mod,      only: delete_ray, t_ray
use msgwam_utils_mod,     only: get_interp_coeffs

implicit none
private

public apply_breaking, apply_dissipation, check_boundaries

contains

 subroutine apply_breaking(z_faces, rho, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max + 1, i_max, j_max),    intent(in)    :: z_faces
    real, dimension(q_max, i_max, j_max),        intent(in)    :: rho
    type(t_ray), dimension(n_max, i_max, j_max), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n, q
    real :: dz, frac, max_kappa, S, wvn_hor_sq, wvn_ver_sq, z_hi, z_lo
    real, dimension(q_max) :: num, den, kappa
    real, dimension(n_max) :: wvn_sq

    ! --------------------------------------------------------------------------

    if (.not. break_waves) then
        return
    end if

    do j = 1, j_max
        do i = 1, i_max
            num = -0.5 * rho(:, i, j)
            den = 0.

            do n = 1, n_max
                if (rays(n, i, j)%meta == -1) then
                    cycle
                end if

                associate(ray => rays(n, i, j))
                    wvn_hor_sq = ray%k ** 2 + ray%l ** 2
                    wvn_ver_sq = ray%m ** 2 + ray%G2

                    wvn_sq(n) = wvn_hor_sq + wvn_ver_sq
                    S = ray%dens * ray%dm * wvn_ver_sq * wvn_hor_sq / &
                        (ray%omega_hat * wvn_sq(n))

                    do q = max(1, ray%q_hi - 1), min(q_max, ray%q_lo + 1)

                        z_hi = z_faces(q, i, j)
                        z_lo = z_faces(q + 1, i, j)
                        dz = z_hi - z_lo

                        if (ray%r_hi < z_lo) then
                            cycle
                        else if (ray%r_lo > z_hi) then
                            exit
                        end if

                        frac = (min(ray%r_hi, z_hi) - max(ray%r_lo, z_lo)) / dz
                        num(q) = num(q) + frac * S
                        den(q) = den(q) + frac * S * wvn_sq(n)
                    end do
                end associate
            end do

            kappa = merge(num / den, 0., den /= 0.)

            do n = 1, n_max
                if (rays(n, i, j)%meta == -1) then
                    cycle
                end if

                max_kappa = 0.
                associate(ray => rays(n, i, j))
                    do q = max(1, ray%q_hi - 1), min(q_max, ray%q_lo + 1)
                        z_hi = z_faces(q, i, j)
                        z_lo = z_faces(q + 1, i, j)

                        if (ray%r_hi < z_lo) then
                            cycle
                        else if (ray%r_lo > z_hi) then
                            exit
                        end if

                        if (kappa(q) > max_kappa) then
                            max_kappa = kappa(q)
                        end if
                    end do

                    ray%dens = ray%dens * max(0., 1 - wvn_sq(n) * max_kappa)
                end associate
            end do

        end do
    end do

end subroutine apply_breaking

pure subroutine apply_dissipation(z_centers, z_faces, rho, dt, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:q_max + 1, i_max, j_max),  intent(in)    :: z_centers
    real, dimension(q_max + 1, i_max, j_max),    intent(in)    :: z_faces
    real, dimension(q_max, i_max, j_max),        intent(in)    :: rho
    real,                                        intent(in)    :: dt
    type(t_ray), dimension(n_max, i_max, j_max), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n
    real :: a, b, damping, dz_inv, nu, r, sponge, wvn_sq, z_sponge, z_zero

    ! --------------------------------------------------------------------------

    if (mu == 0.) then
        return
    end if

    do j = 1, j_max
        do i = 1, i_max

            if (n_sponge > 0) then
                z_zero = z_faces(1, i, j)
                z_sponge = z_faces(n_sponge + 1, i, j)
            end if

            associate( &
                z_col => z_centers(1:q_max, i, j), &
                nu_col => mu / rho(:, i, j) &
            )

                do n = 1, n_max
                    if (rays(n, i, j)%meta == -1) then
                        cycle
                    end if

                    associate (ray => rays(n, i, j))
                        r = (ray%r_lo + ray%r_hi) / 2.

                        if (n_sponge > 0 .and. r > z_sponge) then
                            damping = 1 - ray%cg_r * dt / (z_zero - r)
                            damping = min(max(damping, 0.), 1.)

                        else
                            call get_interp_coeffs(z_col, r, ray%q_mid, a, b)
                            nu = a * nu_col(ray%q_mid) + b * nu_col(ray%q_mid + 1)
                            wvn_sq = ray%k ** 2 + ray%l ** 2 + ray%m ** 2 + ray%G2

                            damping = exp(-dt * nu * wvn_sq * ( &
                                1 + f2(j) / ray%omega_hat ** 2))

                        end if

                        ray%dens = ray%dens * damping
                    end associate
                end do

            end associate
        end do
    end do

end subroutine apply_dissipation

pure subroutine check_boundaries(dt, z_centers, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real,                                        intent(in)    :: dt
    real, dimension(0:q_max + 1, i_max, j_max),  intent(in)    :: z_centers
    type(t_ray), dimension(n_max, i_max, j_max), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    logical :: delete
    integer :: i, j, n
    real :: flux, wvn, z_hi, z_lo

    ! --------------------------------------------------------------------------

    do j = 1, j_max
        do i = 1, i_max
            z_hi = z_centers(0, i, j)
            z_lo = z_centers(q_max + 1, i, j)

            do n = 1, n_max
                associate (ray => rays(n, i, j))
                    if ((ray%meta == -1) .or. ray%is_ghost) then
                        cycle
                    end if

                    wvn = sqrt(ray%k ** 2 + ray%l ** 2)
                    flux = wvn * ray%dens * ray%dm * ray%cg_r

                    delete = ray%r_lo > z_hi
                    delete = delete .or. (ray%r_hi < z_lo)
                    delete = delete .or. (abs(flux) < min_flux)
                    delete = delete .or. (ray%age > max_age)
                end associate

                if (delete) then
                    call delete_ray(rays(n, i, j))
                end if
            end do

        end do
    end do

end subroutine check_boundaries

end module msgwam_sinks_mod