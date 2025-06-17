module msgwam_steady_mod

! ==============================================================================
! Implements the steady-state saturation scheme from Bölöni et al. (2021).
! ==============================================================================

use msgwam_constants_mod, only: f2, i_max, j_max, n_max, n_source, n_sponge, &
                                q_max, use_shapiro_filter
use msgwam_rays_mod,      only: get_cg_r, get_m, t_ray
use msgwam_source_mod,    only: q_source, update_launches
use msgwam_utils_mod,     only: get_interp_coeffs, shapiro_filter

implicit none
private

public get_steady_fluxes

contains

subroutine get_steady_fluxes(z_centers, z_faces, u_c, v_c, rho_c, N2_c, G2, &
    dt, rays, ghosts, last_meta, flux_x, flux_y)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:q_max + 1, i_max, j_max),  intent(in)    :: z_centers
    real, dimension(q_max + 1, i_max, j_max),    intent(in)    :: z_faces
    real, dimension(q_max, i_max, j_max),        intent(in)    :: u_c, v_c, &
                                                                  rho_c, N2_c, G2
    real,                                        intent(in)    :: dt
    type(t_ray), dimension(n_max, i_max, j_max), intent(in)    :: rays
    integer, dimension(n_source, i_max, j_max),  intent(in)    :: ghosts
    integer, dimension(i_max, j_max),            intent(inout) :: last_meta
    real, dimension(q_max + 1, i_max, j_max),    intent(out)   :: flux_x, flux_y

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n, q
    real :: action, aflux, dec_x, dec_y, f_abs, m, omega, omega_hat, &
        omega_hat_sq, threshold, wvn_hor_sq

    integer, dimension(i_max, j_max) :: n_added
    type(t_ray), dimension(n_source, i_max, j_max) :: launches
    real, dimension(q_max + 1) :: u_f, v_f, rho_f, N2_f

    ! --------------------------------------------------------------------------

    n_added = 0
    call update_launches(z_centers, u_c, v_c, N2_c, G2, dt, rays, ghosts, &
        last_meta, n_added, launches)

    flux_x = 0.
    flux_y = 0.

    do j = 1, j_max
        f_abs = sqrt(abs(f2(j)))

        do i = 1, i_max

            call interp_to_faces( &
               z_centers(1:q_max, i, j), z_faces(:, i, j), &
               u_c(:, i, j), v_c(:, i, j), rho_c(:, i, j), N2_c(:, i, j), &
               u_f, v_f, rho_f, N2_f &
            )

            do n = 1, n_source
            associate(ray => launches(n, i, j))
                wvn_hor_sq = ray%k ** 2 + ray%l ** 2
                action = ray%dens * ray%dm
                aflux = action * ray%cg_r

                omega = ray%omega_hat &
                    + ray%k * u_f(q_source + 1) &
                    + ray%l * v_f(q_source + 1)

                do q = q_source + 1, 1, -1
                    omega_hat = omega - ray%k * u_f(q) - ray%l * v_f(q)
                    omega_hat_sq = omega_hat ** 2

                    if (omega_hat <= f_abs) then
                        exit
                    end if

                    if (omega_hat_sq >= N2_f(q)) then
                        flux_x(q + 1:q_max + 1, i, j) = &
                            flux_x(q + 1:q_max + 1, i, j) - ray%k * aflux

                        flux_y(q + 1:q_max + 1, i, j) = &
                            flux_y(q + 1:q_max + 1, i, j) - ray%l * aflux

                        exit
                    end if

                    m = get_m(ray%k, ray%l, omega_hat_sq, N2_f(q), f2(j))
                    threshold = 0.5 * rho_f(q) * omega_hat * ( &
                        1 / (m ** 2) + 1 / wvn_hor_sq)

                    action = min(threshold, action)

                    aflux = -m * ( &
                        (omega_hat_sq - f2(j)) / &
                        omega_hat / (wvn_hor_sq + m ** 2) &
                    ) * action  

                    flux_x(q, i, j) = flux_x(q, i, j) + ray%k * aflux
                    flux_y(q, i, j) = flux_y(q, i, j) + ray%l * aflux                    
                end do

                if (n_sponge > 0) then
                    if ((q == 0) .and. (aflux > 0)) then
                        dec_x = ray%k * aflux / n_sponge
                        dec_y = ray%l * aflux / n_sponge

                        do q = 1, n_sponge
                            flux_x(q, i, j) = flux_x(q, i, j) &
                                - (n_sponge - q + 1) * dec_x

                            flux_y(q, i, j) = flux_y(q, i, j) &
                                - (n_sponge - q + 1) * dec_y

                        end do

                    end if
                end if

            end associate
            end do

        end do
    end do

    if (use_shapiro_filter) then
        call shapiro_filter(flux_x)
        call shapiro_filter(flux_y)
    end if

end subroutine get_steady_fluxes

subroutine interp_to_faces(z_centers, z_faces, u_c, v_c, rho_c, N2_c, &
    u_f, v_f, rho_f, N2_f)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max),     intent(in)  :: z_centers
    real, dimension(q_max + 1), intent(in)  :: z_faces
    real, dimension(q_max),     intent(in)  :: u_c, v_c, rho_c, N2_c
    real, dimension(q_max + 1), intent(out) :: u_f, v_f, rho_f, N2_f

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: q, q_guess
    real :: a, b

    real, dimension(q_max - 1) :: dz_inv

    ! --------------------------------------------------------------------------

    dz_inv = 1. / (z_centers(:q_max - 1) - z_centers(2:))

    do q = 1, q_max + 1
        q_guess = min(q_max - 1, max(1, q - 1))
        call get_interp_coeffs(z_centers, dz_inv, z_faces(q), q_guess, a, b)

        u_f(q) = a * u_c(q_guess) + b * u_c(q_guess + 1)
        v_f(q) = a * v_c(q_guess) + b * v_c(q_guess + 1)

        rho_f(q) = a * rho_c(q_guess) + b * rho_c(q_guess + 1)
        N2_f(q) = a * N2_c(q_guess) + b * N2_c(q_guess + 1)
    end do

end subroutine interp_to_faces

end module msgwam_steady_mod