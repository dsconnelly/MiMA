module msgwam_rays_mod

! ==============================================================================
! This module implements the t_ray types which stores ray volume information. It 
! also implements a range of functions related to the dispersion relation.
! ==============================================================================

implicit none
private

public delete_ray, get_cg_r, get_dm, get_m, &
       get_omega_hat_sq, t_ray

type :: t_ray
    real :: r_hi, r_lo, k, l, m, dm, dens, cg_r, omega_hat, G2
    integer :: age, meta, q_hi, q_lo, q_mid, ghost_id
end type t_ray

contains

pure subroutine delete_ray(ray)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    type(t_ray), intent(inout) :: ray

    ! --------------------------------------------------------------------------

    ray%r_lo = 0.
    ray%r_hi = 0.

    ray%k = 0.
    ray%l = 0.
    ray%m = 0.

    ray%dm = 0.
    ray%dens = 0.

    ray%cg_r = 0.
    ray%omega_hat = 0.
    ray%G2 = 0.

    ray%age = 0
    ray%meta = -1
    ray%ghost_id = -1

    ray%q_lo = 0
    ray%q_mid = 0
    ray%q_hi = 0

end subroutine delete_ray

pure subroutine get_cg_r(m, wvn_hor_sq, m2, N2, f2, G2, &
    omega_hat, cg_r)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, intent(in)  :: m, wvn_hor_sq, m2, N2, f2, G2
    real, intent(out) :: omega_hat, cg_r

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: omega_hat_sq, K2pG2_inv

    ! --------------------------------------------------------------------------

    call get_omega_hat_sq(wvn_hor_sq, m2, N2, f2, G2, K2pG2_inv, omega_hat_sq)
    
    omega_hat = sqrt(omega_hat_sq)
    cg_r = m * (f2 - omega_hat_sq) * K2pG2_inv / omega_hat

end subroutine get_cg_r

pure function get_dm(m, dc, N2) result(dm)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: m, dc, N2
    real             :: dm

    ! --------------------------------------------------------------------------

    dm = dc * (m ** 2) / sqrt(N2)

end function get_dm

pure function get_m(k, l, omega_hat_sq, N2, f2) result(m)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: k, l, omega_hat_sq, N2, f2
    real             :: m

    ! --------------------------------------------------------------------------

    m = -sqrt( &
        (k ** 2 + l ** 2) * (N2 - omega_hat_sq) / &
        (omega_hat_sq - f2) &
    )

end function get_m

pure subroutine get_omega_hat_sq(wvn_hor_sq, m2, N2, f2, G2, &
    K2pG2_inv, omega_hat_sq)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, intent(in)  :: wvn_hor_sq, m2, N2, f2, G2
    real, intent(out) :: K2pG2_inv, omega_hat_sq

    ! --------------------------------------------------------------------------

    K2pG2_inv = 1. / (wvn_hor_sq + m2 + G2)
    omega_hat_sq = (N2 * wvn_hor_sq + f2 * (m2 + G2)) * K2pG2_inv

end subroutine get_omega_hat_sq

end module msgwam_rays_mod