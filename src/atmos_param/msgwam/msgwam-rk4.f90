module msgwam_RK4_mod

! ==============================================================================
! This module implements the RK4 time stepping scheme for the ray volumes.
! ==============================================================================

use msgwam_constants_mod, only: dr_max, dr_min, f2, n_max, q_max
use msgwam_mean_mod,      only: update_mean_gradients
use msgwam_rays_mod,      only: get_cg_r, get_omega_hat_sq, t_ray
use msgwam_source_mod,    only: r_source
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

pure subroutine take_RK4_step(j, z, u, v, N2, G2, dt, rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                       intent(in)    :: j
    real, dimension(q_max),        intent(in)    :: z, u, v, N2, G2
    real,                          intent(in)    :: dt
    type(t_ray), dimension(n_max), intent(inout) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: n
    real, dimension(q_max - 1) :: dz_inv
    real, dimension(q_max) :: du_dz, dv_dz, dN2_dz, dG2_dz

    ! --------------------------------------------------------------------------

    dz_inv = 1. / (z(:q_max - 1) - z(2:))
    call update_mean_gradients(z, u, v, N2, G2, du_dz, dv_dz, dN2_dz, dG2_dz)

    do n = 1, n_max
        if (rays(n)%meta == -1) then
            cycle
        end if

        call take_RK4_step_ray(j, z, dz_inv, N2, G2, du_dz, dv_dz, &
            dN2_dz, dG2_dz, dt, rays(n))
    end do

end subroutine take_RK4_step

pure subroutine take_RK4_step_ray(j, z, dz_inv, N2_col, G2_col, &
    du_dz, dv_dz, dN2_dz, dG2_dz, dt, ray)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                    intent(in)    :: j
    real, dimension(q_max),     intent(in)    :: z, N2_col, G2_col, &
                                                 du_dz, dv_dz, dN2_dz, dG2_dz
    real, dimension(q_max - 1), intent(in)    :: dz_inv
    real,                       intent(in)    :: dt
    type(t_ray),                intent(inout) :: ray

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: stage
    type(t_inc), dimension(0:4) :: incs

    real :: a, area, b, cg_hi, cg_lo, dG2_dr, dN2_dr, du_dr, dv_dr, G2, ignore, &
        K2pG2_inv, m2, N2, omega_hat_sq, r, swap_r, wvn_hor_sq

    ! --------------------------------------------------------------------------

    wvn_hor_sq = ray%k ** 2 + ray%l ** 2
    area = ray%dm * (ray%r_hi - ray%r_lo)
    r = 0.5 * (ray%r_lo + ray%r_hi)
    call locate_all(z, r, ray)

    incs(0)%r_lo = ray%r_lo
    incs(0)%r_hi = ray%r_hi
    incs(0)%m = ray%m

    do stage = 1, 4
        m2 = ray%m ** 2

        call get_interp_coeffs(z, dz_inv, ray%r_hi, ray%q_hi, a, b)
        N2 = a * N2_col(ray%q_hi) + b * N2_col(ray%q_hi + 1)
        G2 = a * G2_col(ray%q_hi) + b * G2_col(ray%q_hi + 1)
        call get_cg_r(ray%m, wvn_hor_sq, m2 + G2, N2, f2(j), ignore, cg_hi)

        call get_interp_coeffs(z, dz_inv, ray%r_lo, ray%q_lo, a, b)
        N2 = a * N2_col(ray%q_lo) + b * N2_col(ray%q_lo + 1)
        G2 = a * G2_col(ray%q_lo) + b * G2_col(ray%q_lo + 1)
        call get_cg_r(ray%m, wvn_hor_sq, m2 + G2, N2, f2(j), ignore, cg_lo)

        if (r < r_source(j)) then
            incs(stage)%r_lo = dt * 0.5 * (cg_lo + cg_hi)
            incs(stage)%r_hi = incs(stage)%r_lo
            incs(stage)%m = 0.
        else

            incs(stage)%r_lo = dt * cg_lo
            incs(stage)%r_hi = dt * cg_hi

            call get_interp_coeffs(z, dz_inv, r, ray%q_mid, a, b)
            N2 = a * N2_col(ray%q_mid) + b * N2_col(ray%q_mid + 1)
            G2 = a * G2_col(ray%q_mid) + b * G2_col(ray%q_mid + 1)
            call get_omega_hat_sq(wvn_hor_sq, m2 + G2, N2, f2(j), &
                K2pG2_inv, omega_hat_sq)

            du_dr = a * du_dz(ray%q_mid) + b * du_dz(ray%q_mid + 1)
            dv_dr = a * dv_dz(ray%q_mid) + b * dv_dz(ray%q_mid + 1)
            dN2_dr = a * dN2_dz(ray%q_mid) + b * dN2_dz(ray%q_mid + 1)
            dG2_dr = a * dG2_dz(ray%q_mid) + b * dG2_dz(ray%q_mid + 1)

            incs(stage)%m = -dt * ( &
                ray%k * du_dr + ray%l * dv_dr + ( &
                    wvn_hor_sq * dN2_dr + &
                    (f2(j) - omega_hat_sq) * dG2_dr &
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
    call get_cg_r(ray%m, wvn_hor_sq, ray%m ** 2 + G2, N2, f2(j), &
        ray%omega_hat, ray%cg_r)

end subroutine take_RK4_step_ray

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