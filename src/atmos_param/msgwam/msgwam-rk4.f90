module msgwam_RK4_mod

! ==============================================================================
! This module implements the RK4 time stepping scheme for the ray volumes.
! ==============================================================================

use msgwam_constants_mod, only: dr_min, f2, Gamma2, i_max, j_max, n_max, q_max
use msgwam_rays_mod,      only: get_cg_r, get_omega_hat, t_ray
use msgwam_utils_mod,     only: interp, locate

implicit none
private

public take_RK4_step

type :: t_inc
    real :: m, r_hi, r_lo
end type t_inc

real, parameter :: DECAY_RATE = 0.01
real, dimension(4), parameter :: COEFFS = (/ 1 / 2., 1 / 2., 1., 1 / 6. /)

contains

 subroutine take_RK4_step(z_centers, u_bar, v_bar, N2, dt, rays)
    use fms_mod, only: mpp_pe
    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:q_max + 1, i_max, j_max),  intent(in)    :: z_centers
    real, dimension(q_max, i_max, j_max),        intent(in)    :: u_bar, &
                                                                  v_bar, N2
    real,                                        intent(in)    :: dt
    type(t_ray), dimension(n_max, i_max, j_max), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n, q, stage, swap_q
    type(t_inc), dimension(0:4) :: incs

    real, dimension(q_max - 1, i_max, j_max) :: du_dr, dv_dr, dN2_dr
    real :: area, cg_hi, cg_lo, cg_mid, dz, N2_hi, N2_lo, N2_mid, &
        omega_hat, r, swap_r, weight, wvn_hor_sq, wvn_ver_sq

    ! --------------------------------------------------------------------------

    do j = 1, j_max
        do i = 1, i_max
            do q = 1, q_max - 1
                dz = z_centers(q, i, j) - z_centers(q + 1, i, j)
                du_dr(q, i, j) = (u_bar(q, i, j) - u_bar(q + 1, i, j)) / dz
                dv_dr(q, i, j) = (v_bar(q, i, j) - v_bar(q + 1, i, j)) / dz
                dN2_dr(q, i, j) = (N2(q, i, j) - N2(q + 1, i, j)) / dz
            end do
        end do
    end do

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_max

                if (rays(n, i, j)%meta == -1) then
                    cycle
                end if

                associate( &
                    ray => rays(n, i, j), &
                    N2_col => N2(:, i, j), &
                    z => z_centers(1:q_max, i, j) &
                )

                    incs(0)%r_lo = ray%r_lo
                    incs(0)%r_hi = ray%r_hi
                    incs(0)%m = ray%m

                    area = (ray%r_hi - ray%r_lo) * ray%dm
                    r = (ray%r_lo + ray%r_hi) / 2.

                    do stage = 1, 4
                        N2_mid = interp(z, N2_col, r, ray%q_mid)
                        N2_lo = interp(z, N2_col, ray%r_lo, min(q_max, ray%q_lo))
                        N2_hi = interp(z, N2_col, ray%r_hi, max(1, ray%q_hi))

                        cg_lo = get_cg_r(ray, N2_lo, f2(j), Gamma2)
                        cg_hi = get_cg_r(ray, N2_hi, f2(j), Gamma2)
                        cg_mid = (cg_hi + cg_lo) / 2.

                        if (.not. ray%is_ghost) then
                            weight = min(1., exp(-DECAY_RATE &
                                * (ray%r_hi - ray%r_lo - dr_min)))

                            incs(stage)%r_lo = dt * (cg_mid * weight + &
                                cg_lo * (1 - weight))

                            incs(stage)%r_hi = dt * (cg_mid * weight + &
                                cg_hi * (1 - weight))

                            wvn_hor_sq = ray%k ** 2 + ray%l ** 2
                            wvn_ver_sq = ray%m ** 2 + Gamma2

                            omega_hat = get_omega_hat( &
                                ray, N2_mid, f2(j), Gamma2)

                            incs(stage)%m = -dt * ( &
                                du_dr(ray%q_mid, i, j) * ray%k + &
                                dv_dr(ray%q_mid, i, j) * ray%l + &
                                wvn_hor_sq * dN2_dr(ray%q_mid, i, j) / &
                                (2 * omega_hat * (wvn_hor_sq + wvn_ver_sq)) &
                            )

                        else
                            incs(stage)%r_lo = dt * cg_mid
                            incs(stage)%r_hi = dt * cg_mid
                            incs(stage)%m = 0
                        end if

                        ray%r_lo = incs(0)%r_lo + &
                                COEFFS(stage) * incs(stage)%r_lo

                        ray%r_hi = incs(0)%r_hi + &
                            COEFFS(stage) * incs(stage)%r_hi

                        ray%m = incs(0)%m + COEFFS(stage) * incs(stage)%m

                        if (stage == 4) then
                            ray%r_lo = ray%r_lo + (incs(1)%r_lo + &
                                2 * (incs(2)%r_lo + incs(3)%r_lo)) / 6.

                            ray%r_hi = ray%r_hi + (incs(1)%r_hi + &
                                2 * (incs(2)%r_hi + incs(3)%r_hi)) / 6.

                            ray%m = ray%m + (incs(1)%m + &
                                2 * (incs(2)%m + incs(3)%m)) / 6.
                        end if

                        r = (ray%r_lo + ray%r_hi) / 2.
                        ray%q_lo = locate(ray%r_lo, z, ray%q_lo)
                        ray%q_hi = locate(ray%r_hi, z, ray%q_hi)
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

                    ray%age = ray%age + int(dt)
                    ray%dm = area / (ray%r_hi - ray%r_lo)

                    N2_mid = interp(z, N2_col, r, ray%q_mid)
                    ray%cg_r = get_cg_r(ray, N2_mid, f2(j), Gamma2)
                    ray%omega_hat = get_omega_hat(ray, N2_mid, f2(j), Gamma2)

                end associate

            end do
        end do
    end do

end subroutine take_RK4_step

end module msgwam_RK4_mod