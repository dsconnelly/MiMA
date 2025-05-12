module msgwam_RK4_mod

! ==============================================================================
! This module implements the RK4 time stepping scheme for the ray volumes.
! ==============================================================================

use msgwam_constants_mod, only: dr_max, dr_min, f2, i_max, j_max, n_max, q_max
use msgwam_mean_mod,      only: update_mean_gradients
use msgwam_rays_mod,      only: get_cg_r, get_omega_hat_sq, t_ray
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

pure subroutine take_RK4_step(z_centers, u_bar, v_bar, N2_all, G2_all, dt, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:q_max + 1, i_max, j_max),  intent(in)    :: z_centers
    real, dimension(q_max, i_max, j_max),        intent(in)    :: u_bar, &
                                                                  v_bar, &
                                                                  N2_all, &
                                                                  G2_all
    real,                                        intent(in)    :: dt
    type(t_ray), dimension(n_max, i_max, j_max), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n, stage

    real, dimension(q_max - 1) :: dz_inv
    real, dimension(q_max) :: dG2_dr, dN2_dr, du_dr, dv_dr
    real :: a, area, b, cg_hi, cg_lo, cg_mid, dG2_mid, dN2_mid, du_mid, dv_mid, &
        G2, ignore, K2pG2_inv, m2, N2, omega_hat_sq, r, swap_r, wvn_hor_sq

    type(t_inc), dimension(0:4) :: incs

    ! --------------------------------------------------------------------------

    do j = 1, j_max
        do i = 1, i_max
            associate( &
                z => z_centers(1:q_max, i, j), &
                N2_col => N2_all(:, i, j), &
                G2_col => G2_all(:, i, j) &
            )

            dz_inv = 1. / (z(:q_max - 1) - z(2:))
            call update_mean_gradients(z, u_bar(:, i, j), v_bar(:, i, j), &
                N2_col, G2_col, du_dr, dv_dr, dN2_dr, dG2_dr)

            do n = 1, n_max
                if (rays(n, i, j)%meta == -1) then
                    cycle
                end if

                associate(ray => rays(n, i, j))

                    incs(0)%r_lo = ray%r_lo
                    incs(0)%r_hi = ray%r_hi
                    incs(0)%m = ray%m

                    wvn_hor_sq = ray%k ** 2 + ray%l ** 2
                    area = ray%dm * (ray%r_hi - ray%r_lo)
                    r = 0.5 * (ray%r_lo + ray%r_hi)
                    call locate_all(z, r, ray)

                    do stage = 1, 4

                        m2 = ray%m ** 2

                        call get_interp_coeffs(z, dz_inv, ray%r_hi, ray%q_hi, a, b)
                        N2 = a * N2_col(ray%q_hi) + b * N2_col(ray%q_hi + 1)
                        G2 = a * G2_col(ray%q_hi) + b * G2_col(ray%q_hi + 1)
                        call get_cg_r(ray%m, wvn_hor_sq, m2, N2, f2(j), G2, &
                            ignore, cg_hi)

                        call get_interp_coeffs(z, dz_inv, ray%r_lo, ray%q_lo, a, b)
                        N2 = a * N2_col(ray%q_lo) + b * N2_col(ray%q_lo + 1)
                        G2 = a * G2_col(ray%q_lo) + b * G2_col(ray%q_lo + 1)
                        call get_cg_r(ray%m, wvn_hor_sq, m2, N2, f2(j), G2, &
                            ignore, cg_lo)

                        if (ray%is_ghost) then
                            cg_mid = 0.5 * (cg_lo + cg_hi)
                            incs(stage)%r_lo = dt * cg_mid
                            incs(stage)%r_hi = dt * cg_mid
                            incs(stage)%m = 0.
                        else
                            incs(stage)%r_lo = dt * cg_lo
                            incs(stage)%r_hi = dt * cg_hi

                            call get_interp_coeffs(z, dz_inv, r, ray%q_mid, a, b)
                            N2 = a * N2_col(ray%q_mid) + b * N2_col(ray%q_mid + 1)
                            G2 = a * G2_col(ray%q_mid) + b * G2_col(ray%q_mid + 1)
                            call get_omega_hat_sq(wvn_hor_sq, m2, N2, f2(j), &
                                G2, K2pG2_inv, omega_hat_sq)

                            dG2_mid = a * dG2_dr(ray%q_mid) + b * dG2_dr(ray%q_mid + 1)
                            dN2_mid = a * dN2_dr(ray%q_mid) + b * dN2_dr(ray%q_mid + 1)
                            du_mid = a * du_dr(ray%q_mid) + b * du_dr(ray%q_mid + 1)
                            dv_mid = a * dv_dr(ray%q_mid) + b * dv_dr(ray%q_mid + 1)

                            incs(stage)%m = -dt * ( &
                                ray%k * du_mid + ray%l * dv_mid + ( &
                                    wvn_hor_sq * dN2_mid + &
                                    (f2(j) - omega_hat_sq) * dG2_mid &
                                ) * 0.5 * K2pG2_inv / sqrt(omega_hat_sq) &
                            )
                        end if

                        call update_stage(stage, incs(:)%r_hi, ray%r_hi)
                        call update_stage(stage, incs(:)%r_lo, ray%r_lo)
                        call update_stage(stage, incs(:)%m, ray%m)

                        if (stage < 4) then
                            r = 0.5 * (ray%r_lo + ray%r_hi)
                            call locate_all(z, r, ray)
                        end if

                    end do

                    call update_stage_final(incs(:)%r_lo, ray%r_lo)
                    call update_stage_final(incs(:)%r_hi, ray%r_hi)
                    call update_stage_final(incs(:)%m, ray%m)
                    r = 0.5 * (ray%r_lo + ray%r_hi)

                    if (ray%r_hi < ray%r_lo) then
                        swap_r = ray%r_hi
                        ray%r_hi = ray%r_lo
                        ray%r_lo = swap_r
                    end if

                    if (ray%r_hi - ray%r_lo < dr_min) then
                        ray%r_hi = r + 0.5 * dr_min
                        ray%r_lo = r - 0.5 * dr_min
                    else if (ray%r_hi - ray%r_lo > dr_max) then
                        ray%r_hi = r + 0.5 * dr_max
                        ray%r_lo = r - 0.5 * dr_max
                    end if

                    ray%age = ray%age + int(dt)
                    ray%dm = area / (ray%r_hi - ray%r_lo)
                    call locate_all(z, r, ray)

                    call get_interp_coeffs(z, dz_inv, r, ray%q_mid, a, b)
                    N2 = a * N2_col(ray%q_mid) + b * N2_col(ray%q_mid + 1)
                    G2 = a * G2_col(ray%q_mid) + b * G2_col(ray%q_mid + 1)

                    ray%G2 = G2
                    call get_cg_r(ray%m, wvn_hor_sq, ray%m ** 2, N2, f2(j), &
                        G2, ray%omega_hat, ray%cg_r)

                end associate
            end do

            end associate
        end do
    end do

end subroutine take_RK4_step

pure subroutine locate_all(z, r, ray)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max), intent(in)    :: z
    real,                   intent(in)    :: r
    type(t_ray),            intent(inout) :: ray

    ! --------------------------------------------------------------------------

    ray%q_mid = locate(z, r, ray%q_mid)
    ray%q_lo = locate(z, ray%r_lo, ray%q_lo)
    ray%q_hi = locate(z, ray%r_hi, ray%q_hi)    

end subroutine locate_all

pure subroutine update_stage(stage, incs, out)
    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,              intent(in)  :: stage
    real, dimension(0:4), intent(in)  :: incs
    real,                 intent(out) :: out

    ! --------------------------------------------------------------------------

    out = incs(0) + COEFFS(stage) * incs(stage)

end subroutine update_stage

pure subroutine update_stage_final(incs, out)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(0:4), intent(in)    :: incs
    real,                 intent(inout) :: out

    ! --------------------------------------------------------------------------

    out = out + ONE_SIXTH * incs(1) + ONE_THIRD * (incs(2) + incs(3))

end subroutine update_stage_final

end module msgwam_RK4_mod