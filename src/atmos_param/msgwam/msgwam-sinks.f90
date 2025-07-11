module msgwam_sinks_mod

! ==============================================================================
! This module implements ray volume dissipation and breaking.
! ==============================================================================

use msgwam_constants_mod, only: break_waves, f2, max_age, min_flux, mu, n_max, &
                                n_sponge, q_max
use msgwam_rays_mod,      only: delete_ray, t_ray
use msgwam_utils_mod,     only: get_interp_coeffs

implicit none
private

public apply_breaking, apply_dissipation, check_boundaries

contains

pure subroutine apply_breaking(z_faces, rho, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max + 1),    intent(in) :: z_faces
    real, dimension(q_max),        intent(in) :: rho
    type(t_ray), dimension(n_max), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: n, q, q_hi, q_lo
    real :: D, dz, frac, max_kappa, wvn_hor_sq, wvn_ver_sq, z_hi, z_lo
    
    real, dimension(q_max) :: num, den, kappa
    real, dimension(n_max) :: wvn_sq

    ! --------------------------------------------------------------------------

    if (.not. break_waves) then
        return
    end if

    num = -0.5 * rho
    den = 0.

    do n = 1, n_max
        if (rays(n)%meta == -1) then
            cycle
        end if

        wvn_hor_sq = rays(n)%k ** 2 + rays(n)%l ** 2
        wvn_ver_sq = rays(n)%m ** 2 + rays(n)%G2
        wvn_sq(n) = wvn_hor_sq + wvn_ver_sq

        D = rays(n)%dens * rays(n)%dm * wvn_hor_sq * wvn_ver_sq / &
            rays(n)%omega_hat

        q_hi = max(1, rays(n)%q_hi - 1)
        q_lo = min(q_max, rays(n)%q_lo + 1)

        do q = q_hi, q_lo
            z_hi = z_faces(q)
            z_lo = z_faces(q + 1)
            dz = z_hi - z_lo

            if (rays(n)%r_hi < z_lo) then
                cycle
            else if (rays(n)%r_lo > z_hi) then
                exit
            end if

            frac = (min(rays(n)%r_hi, z_hi) - max(rays(n)%r_lo, z_lo)) / dz
            num(q) = num(q) + frac * D / wvn_sq(n)
            den(q) = den(q) + frac * D
        end do
    end do

    kappa = merge(num / den, 0., den /= 0.)

    do n = 1, n_max
        if (rays(n)%meta == -1) then
            cycle
        end if

        max_kappa = 0.
        q_hi = max(1, rays(n)%q_hi - 1)
        q_lo = min(q_max, rays(n)%q_lo + 1)

        do q = q_hi, q_lo
            z_hi = z_faces(q)
            z_lo = z_faces(q + 1)

            if (rays(n)%r_hi < z_lo) then
                cycle
            else if (rays(n)%r_lo > z_hi) then
                exit
            end if

            if (kappa(q) > max_kappa) then
                max_kappa = kappa(q)
            end if
        end do

        rays(n)%dens = rays(n)%dens * max(0., 1 - wvn_sq(n) * max_kappa)
    end do

end subroutine apply_breaking

pure subroutine apply_dissipation(j, z_faces, z, rho, dt, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                       intent(in)    :: j
    real, dimension(q_max + 1),    intent(in)    :: z_faces
    real, dimension(q_max),        intent(in)    :: z, rho
    real,                          intent(in)    :: dt
    type(t_ray), dimension(n_max), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: n
    real :: a, b, damping, nu, r, wvn_sq, z_sponge, z_zero
    real, dimension(q_max - 1) :: dz_inv
    real, dimension(q_max) :: rho_inv

    ! --------------------------------------------------------------------------

    if (mu == 0.) then
        return
    end if

    rho_inv = 1. / rho
    dz_inv = 1. / (z(:q_max - 1) - z(2:))

    if (n_sponge > 0) then
        z_zero = z_faces(1)
        z_sponge = z_faces(n_sponge + 1)
    end if

    do n = 1, n_max
        if (rays(n)%meta == -1) then
            cycle
        end if

        r = 0.5 * (rays(n)%r_lo + rays(n)%r_hi)
        if (n_sponge > 0 .and. r > z_sponge) then
            damping = 1 - rays(n)%cg_r * dt / (z_zero - r)
            damping = min(1., max(0., damping))
        
        else
            call get_interp_coeffs(z, dz_inv, r, rays(n)%q_mid, a, b)
            nu = a * rho_inv(rays(n)%q_mid) + b * rho_inv(rays(n)%q_mid + 1)
            nu = mu * nu

            wvn_sq = rays(n)%k ** 2 + rays(n)%l ** 2 &
                + rays(n)%m ** 2 + rays(n)%G2

            damping = exp(-dt * nu * wvn_sq * ( &
                1 + f2(j) / rays(n)%omega_hat ** 2))
        end if

        rays(n)%dens = rays(n)%dens * damping
    end do

end subroutine apply_dissipation

pure subroutine check_boundaries(z_centers, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:q_max + 1),  intent(in)    :: z_centers
    type(t_ray), dimension(n_max), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: n
    logical :: delete
    real :: flux, wvn, z_hi, z_lo

    ! --------------------------------------------------------------------------

    z_hi = z_centers(0)
    z_lo = z_centers(q_max + 1)

    do n = 1, n_max
        if ((rays(n)%meta == -1) .or. (rays(n)%ghost_id /= -1)) then
            cycle
        end if

        wvn = sqrt(rays(n)%k ** 2 + rays(n)%l ** 2)
        flux = wvn * rays(n)%dens * rays(n)%dm * rays(n)%cg_r

        delete = rays(n)%r_lo > z_hi
        delete = delete .or. (rays(n)%r_hi < z_lo)
        delete = delete .or. (abs(flux) < min_flux)
        delete = delete .or. (rays(n)%age > max_age)

        if (delete) then
            call delete_ray(rays(n))
        end if
    end do

end subroutine check_boundaries

end module msgwam_sinks_mod