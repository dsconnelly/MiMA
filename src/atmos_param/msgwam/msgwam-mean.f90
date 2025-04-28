module msgwam_mean_mod

! ==============================================================================
! This module provides a subroutine that calculates the required mean state
! variables with axis order optimized for performance.
! ==============================================================================

use constants_mod,        only: CP_AIR, GRAV, RDGAS

use msgwam_constants_mod, only: i_max, j_max, min_N2, n_max, n_sponge, q_max, &
                                use_shapiro_filter
use msgwam_rays_mod,      only: t_ray
use msgwam_utils_mod,     only: shapiro_filter

implicit none
private

public get_accelerations, project_fluxes, update_mean_fields, &
       update_mean_gradients

contains

pure subroutine get_accelerations(z_faces, rho, flux_x, flux_y, &
    sponge_x, sponge_y, du_dt, dv_dt)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max + 1, i_max, j_max), intent(in)  :: z_faces
    real, dimension(q_max, i_max, j_max),     intent(in)  :: rho
    real, dimension(q_max + 1, i_max, j_max), intent(in)  :: flux_x, flux_y
    real, dimension(i_max, j_max),            intent(in)  :: sponge_x, sponge_y
    real, dimension(i_max, j_max, q_max),     intent(out) :: du_dt, dv_dt

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

            if (n_sponge > 0) then
                do q = 1, n_sponge
                    du_dt(i, j, q) = du_dt(i, j, q) + sponge_x(i, j)
                    dv_dt(i, j, q) = dv_dt(i, j, q) + sponge_y(i, j)
                end do
            end if

        end do
    end do

end subroutine get_accelerations

pure subroutine project_fluxes(z_centers, rays, flux_x, flux_y)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:q_max + 1, i_max, j_max),  intent(in)  :: z_centers
    type(t_ray), dimension(n_max, i_max, j_max), intent(in)  :: rays
    real, dimension(q_max + 1, i_max, j_max),    intent(out) :: flux_x, flux_y

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n, q
    real :: aflux, dz, frac, z_hi, z_lo

    ! --------------------------------------------------------------------------

    flux_x = 0.
    flux_y = 0.

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_max

                if (rays(n, i, j)%meta == -1) then
                    cycle
                end if

                associate(ray => rays(n, i, j))
                    aflux = ray%dens * ray%dm * ray%cg_r

                    do q = ray%q_hi, ray%q_lo + 2
                        z_hi = z_centers(q - 1, i, j)
                        z_lo = z_centers(q, i, j)
                        dz = z_hi - z_lo

                        if (ray%r_hi < z_lo) then
                            cycle
                        else if (ray%r_lo > z_hi) then
                            exit
                        end if

                        frac = (min(ray%r_hi, z_hi) - max(ray%r_lo, z_lo)) / dz
                        flux_x(q, i, j) = flux_x(q, i, j) + frac * ray%k * aflux
                        flux_y(q, i, j) = flux_y(q, i, j) + frac * ray%l * aflux
                    end do
                end associate

            end do
        end do
    end do

    if (use_shapiro_filter) then
        call shapiro_filter(flux_x)
        call shapiro_filter(flux_y)
    end if

end subroutine project_fluxes

pure subroutine update_mean_gradients(z, u_bar, v_bar, N2, G2, &
    du_dr, dv_dr, dN2_dr, dG2_dr)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max), intent(in)  :: z
    real, dimension(q_max), intent(in)  :: u_bar, v_bar, N2, G2
    real, dimension(q_max), intent(out) :: du_dr, dv_dr, dN2_dr, dG2_dr

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: q, q_hi, q_lo
    real :: dz

    ! --------------------------------------------------------------------------

    do q = 1, q_max
        q_hi = max(1, q - 1)
        q_lo = min(q_max, q + 1)
        dz = z(q_hi) - z(q_lo)

        du_dr(q) = (u_bar(q_hi) - u_bar(q_lo)) / dz
        dv_dr(q) = (v_bar(q_hi) - v_bar(q_lo)) / dz
        dN2_dr(q) = (N2(q_hi) - N2(q_lo)) / dz
        dG2_dr(q) = (G2(q_hi) - G2(q_lo)) / dz
    end do
    
end subroutine update_mean_gradients

pure subroutine update_mean_fields(z_full, p_full, temp, uuu, vvv, &
    z_centers, z_faces, u_bar, v_bar, rho, N2, G2)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(i_max, j_max, q_max),       intent(in)  :: z_full, p_full, &
                                                               temp, uuu, vvv
    real, dimension(0:q_max + 1, i_max, j_max), intent(out) :: z_centers
    real, dimension(q_max + 1, i_max, j_max),   intent(out) :: z_faces
    real, dimension(q_max, i_max, j_max),       intent(out) :: u_bar, v_bar, &
                                                               rho, N2, G2

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, q, q_hi, q_lo
    real :: dT_dz, drho_dz, dz

    ! --------------------------------------------------------------------------

    do q = 1, q_max
        do j = 1, j_max
            do i = 1, i_max
                u_bar(q, i, j) = uuu(i, j, q)
                v_bar(q, i, j) = vvv(i, j, q)
                z_centers(q, i, j) = z_full(i, j, q)
                rho(q, i, j) = p_full(i, j, q) / RDGAS / temp(i, j, q)
            end do
        end do
    end do

    z_centers(q_max + 1, :, :) = 0.
    z_centers(0, :, :) = 2 * z_centers(1, :, :) - z_centers(2, :, :)
    z_faces = (z_centers(1:, :, :) + z_centers(:q_max, :, :)) / 2.

    do j = 1, j_max
        do i = 1, i_max
            do q = 1, q_max
                q_hi = max(1, q - 1)
                q_lo = min(q_max, q + 1)
                
                dz = z_centers(q_hi, i, j) - z_centers(q_lo, i, j)
                drho_dz = (rho(q_hi, i, j) - rho(q_lo, i, j)) / dz
                dT_dz = (temp(i, j, q_hi) - temp(i, j, q_lo)) / dz

                N2(q, i, j) = max(min_N2, &
                    (GRAV / temp(i, j, q)) * (dT_dz + GRAV / CP_AIR))

                G2(q, i, j) = ((1. / 2. - 2. / 7.) * drho_dz / &
                    (-2. * rho(q, i, j))) ** 2

            end do
        end do
    end do

end subroutine update_mean_fields

end module msgwam_mean_mod