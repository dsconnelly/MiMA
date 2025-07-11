module msgwam_mean_mod

! ==============================================================================
! This module provides a subroutine that calculates the required mean state
! variables with axis order optimized for performance. Also contains routines
! for calculating fluxes and accelerations on the mean grid.
! ==============================================================================

use constants_mod,        only: CP_AIR, GRAV, RDGAS

use msgwam_constants_mod, only: min_N2, n_max, q_max, stoch_factor, &
                                use_shapiro_filter
use msgwam_rays_mod,      only: t_ray
use msgwam_utils_mod,     only: shapiro_filter

implicit none
private

public get_accelerations, project_fluxes, update_mean_fields, &
       update_mean_gradients

contains

pure subroutine get_accelerations(z_faces, rho, flux_x, flux_y, du_dt, dv_dt)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max + 1), intent(in)  :: z_faces, flux_x, flux_y
    real, dimension(q_max),     intent(in)  :: rho
    real, dimension(q_max),     intent(out) :: du_dt, dv_dt

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real, dimension(q_max) :: dz_inv, rho_inv

    ! --------------------------------------------------------------------------

    rho_inv = 1. / rho
    dz_inv = 1. / (z_faces(:q_max) - z_faces(2:))
    du_dt = -(flux_x(:q_max) - flux_x(2:)) * dz_inv * rho_inv
    dv_dt = -(flux_y(:q_max) - flux_y(2:)) * dz_inv * rho_inv

end subroutine get_accelerations

pure subroutine project_fluxes(z_centers, rays, flux_x, flux_y)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:q_max + 1),  intent(in)  :: z_centers
    type(t_ray), dimension(n_max), intent(in)  :: rays
    real, dimension(q_max + 1),    intent(out) :: flux_x, flux_y

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: n, q
    real :: aflux, dz, frac, kflux, lflux, z_hi, z_lo

    ! --------------------------------------------------------------------------

    flux_x = 0.
    flux_y = 0.

    do n = 1, n_max
        if (rays(n)%meta == -1) then
            cycle
        end if

        aflux = rays(n)%dens * rays(n)%dm * rays(n)%cg_r * stoch_factor
        kflux = rays(n)%k * aflux
        lflux = rays(n)%l * aflux

        do q = rays(n)%q_hi, rays(n)%q_lo + 2
            z_hi = z_centers(q - 1)
            z_lo = z_centers(q)

            if (rays(n)%r_hi < z_lo) then
                cycle
            else if (rays(n)%r_lo > z_hi) then
                exit
            end if

            dz = z_hi - z_lo
            frac = (min(rays(n)%r_hi, z_hi) - max(rays(n)%r_lo, z_lo)) / dz
            flux_x(q) = flux_x(q) + frac * kflux
            flux_y(q) = flux_y(q) + frac * lflux
        end do
    end do

    if (use_shapiro_filter) then
        call shapiro_filter(flux_x)
        call shapiro_filter(flux_y)
    end if

end subroutine project_fluxes

pure subroutine update_mean_fields(z_full, p_full, temp, uuu, vvv, &
    z_centers, z_faces, u, v, rho, N2, G2)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max),       intent(in)  :: z_full, p_full, temp, uuu, vvv
    real, dimension(0:q_max + 1), intent(out) :: z_centers
    real, dimension(q_max + 1),   intent(out) :: z_faces
    real, dimension(q_max),       intent(out) :: u, v, rho, N2, G2
    
    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: q, q_hi, q_lo
    real :: dT_dz, drho_dz, dz

    ! --------------------------------------------------------------------------

    u = uuu
    v = vvv
    z_centers(1:q_max) = z_full
    rho = p_full / RDGAS / temp

    z_centers(q_max + 1) = 0.
    z_centers(0) = 2. * z_centers(1) - z_centers(2)
    z_faces = 0.5 * (z_centers(1:) + z_centers(:q_max))

    do q = 1, q_max
        q_hi = max(1, q - 1)
        q_lo = min(q_max, q + 1)

        dz = z_centers(q_hi) - z_centers(q_lo)
        drho_dz = (rho(q_hi) - rho(q_lo)) / dz
        dT_dz = (temp(q_hi) - temp(q_lo)) / dz

        N2(q) = max(min_N2, (GRAV / temp(q)) * (dT_dz + GRAV / CP_AIR))
        G2(q) = ((1. / 2. - 2. / 7.) * drho_dz / (-2. * rho(q))) ** 2
    end do

end subroutine update_mean_fields

pure subroutine update_mean_gradients(z, u, v, N2, G2, &
    du_dz, dv_dz, dN2_dz, dG2_dz)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max), intent(in)  :: z, u, v, N2, G2
    real, dimension(q_max), intent(out) :: du_dz, dv_dz, dN2_dz, dG2_dz

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: q, q_hi, q_lo
    real :: dz_inv

    ! --------------------------------------------------------------------------

    do q = 1, q_max
        q_hi = max(1, q - 1)
        q_lo = min(q_max, q + 1)
        dz_inv = 1. / (z(q_hi) - z(q_lo))

        du_dz(q) = (u(q_hi) - u(q_lo)) * dz_inv
        dv_dz(q) = (v(q_hi) - v(q_lo)) * dz_inv
        dN2_dz(q) = (N2(q_hi) - N2(q_lo)) * dz_inv
        dG2_dz(q) = (G2(q_hi) - G2(q_lo)) * dz_inv
    end do

end subroutine update_mean_gradients

end module msgwam_mean_mod
