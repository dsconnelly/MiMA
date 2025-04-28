module msgwam_RK4_mod

! ==============================================================================
! This module implements the RK4 time stepping scheme for the ray volumes.
! ==============================================================================

use msgwam_constants_mod, only: dr_min, f2, i_max, j_max, n_max, q_max
use msgwam_mean_mod,      only: update_mean_gradients
use msgwam_rays_mod,      only: get_cg_r, get_omega_hat, t_ray
use msgwam_utils_mod,     only: get_interp_coeffs, locate

implicit none
private

public take_RK4_step

type :: t_inc
    real :: r_hi, r_lo, m
end type t_inc

real, parameter :: ONE_SIXTH = 1. / 6.
real, parameter :: ONE_THIRD = 1. / 3.
real, dimension(4), parameter :: COEFFS = (/ 1. / 2, 1. / 2, 1., 1. / 6 /)

contains

pure subroutine take_RK4_step(z_centers, u_bar, v_bar, N2, G2, dt, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:q_max + 1, i_max, j_max),  intent(in)    :: z_centers
    real, dimension(q_max, i_max, j_max),        intent(in)    :: u_bar, &
                                                                  v_bar, N2, G2
    real,                                        intent(in)    :: dt
    type(t_ray), dimension(n_max, i_max, j_max), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n, stage, swap_q
    type(t_inc), dimension(0:4) :: incs

    real, dimension(q_max) :: du_dr, dv_dr, dN2_dr, dG2_dr
    real :: a, area, b, cg_hi, cg_lo, cg_mid, dG2_mid, dN2_mid, du_mid, dv_mid, &
            G2_hi, G2_lo, G2_mid, N2_hi, N2_lo, N2_mid, omega_hat, r, swap_r, &
            wvn_hor_sq, wvn_ver_sq

    ! --------------------------------------------------------------------------

    do j = 1, j_max
        do i = 1, i_max

            associate( &
                N2_col => N2(:, i, j), &
                G2_col => G2(:, i, j), &
                z => z_centers(1:q_max, i, j) &
            )

            call update_mean_gradients(z_centers(1:q_max, i, j), &
                u_bar(:, i, j), v_bar(:, i, j), N2_col, G2_col, &
                du_dr, dv_dr, dN2_dr, dG2_dr)
        
            do n = 1, n_max

                if (rays(n, i, j)%meta == -1) then
                    cycle
                end if

                associate(ray => rays(n, i, j))

                incs(0)%r_lo = ray%r_lo
                incs(0)%r_hi = ray%r_hi
                incs(0)%m = ray%m

                area = (ray%r_hi - ray%r_lo) * ray%dm
                r = 0.5 * (ray%r_hi + ray%r_lo)

                ray%q_hi = locate(ray%r_hi, z, ray%q_hi)
                ray%q_lo = locate(ray%r_lo, z, ray%q_lo)
                ray%q_mid = locate(r, z, ray%q_mid)

                do stage = 1, 4

                    call get_interp_coeffs(z, ray%r_hi, ray%q_hi, a, b)
                    N2_hi = a * N2_col(ray%q_hi) + b * N2_col(ray%q_hi + 1)
                    G2_hi = a * G2_col(ray%q_hi) + b * G2_col(ray%q_hi + 1)

                    call get_interp_coeffs(z, ray%r_lo, ray%q_lo, a, b)
                    N2_lo = a * N2_col(ray%q_lo) + b * N2_col(ray%q_lo + 1)
                    G2_lo = a * G2_col(ray%q_lo) + b * G2_col(ray%q_lo + 1)

                    call get_interp_coeffs(z, r, ray%q_mid, a, b)
                    N2_mid = a * N2_col(ray%q_mid) + b * N2_col(ray%q_mid + 1)
                    G2_mid = a * G2_col(ray%q_mid) + b * G2_col(ray%q_mid + 1)

                    cg_hi = get_cg_r(ray, N2_hi, f2(j), G2_hi)
                    cg_lo = get_cg_r(ray, N2_lo, f2(j), G2_lo)

                    if (.not. ray%is_ghost) then
                        incs(stage)%r_lo = dt * cg_lo
                        incs(stage)%r_hi = dt * cg_hi

                        omega_hat = get_omega_hat(ray, N2_mid, f2(j), G2_mid)
                        wvn_hor_sq = ray%k ** 2 + ray%l ** 2
                        wvn_ver_sq = ray%m ** 2 + G2_mid

                        dG2_mid = a * dG2_dr(ray%q_mid) + b * dG2_dr(ray%q_mid + 1)
                        dN2_mid = a * dN2_dr(ray%q_mid) + b * dN2_dr(ray%q_mid + 1)
                        du_mid = a * du_dr(ray%q_mid) + b * du_dr(ray%q_mid + 1)
                        dv_mid = a * dv_dr(ray%q_mid) + b * dv_dr(ray%q_mid + 1)

                        incs(stage)%m = -dt * ( &
                            ray%k * du_mid + ray%l * dv_mid + &
                            ( &
                                wvn_hor_sq * dN2_mid + &
                                (f2(j) - omega_hat ** 2) * dG2_mid &
                            ) / (2 * omega_hat * (wvn_hor_sq + wvn_ver_sq)) &
                        )

                    else
                        cg_mid = 0.5 * (cg_hi + cg_lo)
                        incs(stage)%r_hi = dt * cg_mid
                        incs(stage)%r_lo = dt * cg_mid
                        incs(stage)%m = 0.                        
                    end if

                    ray%r_lo = incs(0)%r_lo + COEFFS(stage) * incs(stage)%r_lo
                    ray%r_hi = incs(0)%r_hi + COEFFS(stage) * incs(stage)%r_hi
                    ray%m = incs(0)%m + COEFFS(stage) * incs(stage)%m

                    if (stage == 4) then
                        ray%r_lo = ray%r_lo + ONE_SIXTH * incs(1)%r_lo + &
                            ONE_THIRD * (incs(2)%r_lo + incs(3)%r_lo)

                        ray%r_hi = ray%r_hi + ONE_SIXTH * incs(1)%r_hi + &
                            ONE_THIRD * (incs(2)%r_hi + incs(3)%r_hi)

                        ray%m = ray%m + ONE_SIXTH * incs(1)%m + &
                            ONE_THIRD * (incs(2)%m + incs(3)%m)
                    end if

                    r = 0.5 * (ray%r_lo + ray%r_hi)
                    ray%q_hi = locate(ray%r_hi, z, ray%q_hi)
                    ray%q_lo = locate(ray%r_lo, z, ray%q_lo)
                    ray%q_mid = locate(r, z, ray%q_mid)
                end do

                if (ray%r_hi < ray%r_lo) then
                    swap_r = ray%r_hi
                    ray%r_hi = ray%r_lo
                    ray%r_lo = swap_r

                    swap_q = ray%q_hi
                    ray%q_hi = ray%q_lo
                    ray%q_lo = swap_q
                end if

                if (ray%r_hi - ray%r_lo < dr_min) then
                    ray%r_hi = r + 0.5 * dr_min
                    ray%r_lo = r - 0.5 * dr_min

                    ray%q_hi = locate(ray%r_hi, z, ray%q_hi)
                    ray%q_lo = locate(ray%r_lo, z, ray%q_lo)
                end if

                ray%age = ray%age + int(dt)
                ray%dm = area / (ray%r_hi - ray%r_lo)

                call get_interp_coeffs(z, r, ray%q_mid, a, b)
                N2_mid = a * N2_col(ray%q_mid) + b * N2_col(ray%q_mid + 1)
                G2_mid = a * G2_col(ray%q_mid) + b * G2_col(ray%q_mid + 1)

                ray%cg_r = get_cg_r(ray, N2_mid, f2(j), G2_mid)
                ray%omega_hat = get_omega_hat(ray, N2_mid, f2(j), G2_mid)
                ray%G2 = G2_mid

                end associate

            end do

            end associate

        end do
    end do

end subroutine take_RK4_step

end module msgwam_RK4_mod