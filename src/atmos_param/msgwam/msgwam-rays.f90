module msgwam_rays_mod

! ==============================================================================
! This module implements the t_ray types which stores ray volume information. It 
! also implements a range of functions related to the dispersion relation.
! ==============================================================================

implicit none
private

public delete_ray, get_cg_r, get_dm, get_m, get_omega_hat, t_ray

type :: t_ray
    real :: r_hi, r_lo, k, l, m, dm, dens, cg_r, omega_hat, G2
    integer :: age, meta, q_hi, q_lo, q_mid
    logical :: is_ghost
end type t_ray

interface get_cg_r
    module procedure get_cg_r_from_ray
    module procedure get_cg_r_from_reals
end interface get_cg_r

interface get_omega_hat
    module procedure get_omega_hat_from_ray
    module procedure get_omega_hat_from_reals
end interface get_omega_hat

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
    ray%is_ghost = .false.

    ray%q_lo = 0
    ray%q_mid = 0
    ray%q_hi = 0

end subroutine delete_ray

pure function get_cg_r_from_ray(ray, N2, f2, Gamma2) result(cg_r)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    type(t_ray), intent(in) :: ray
    real,        intent(in) :: N2, f2, Gamma2
    real                    :: cg_r

    ! --------------------------------------------------------------------------

    cg_r = get_cg_r_from_reals(ray%k, ray%l, ray%m, N2, f2, Gamma2)

end function get_cg_r_from_ray

pure function get_cg_r_from_reals(k, l, m, N2, f2, Gamma2) result(cg_r)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: k, l, m, N2, f2, Gamma2
    real             :: cg_r

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: omega_hat, wvn_sq

    ! --------------------------------------------------------------------------

    wvn_sq = k ** 2 + l ** 2 + m ** 2 + Gamma2
    omega_hat = get_omega_hat(k, l, m, N2, f2, Gamma2)

    cg_r = -m * (omega_hat ** 2 - f2) / (omega_hat * wvn_sq)

end function get_cg_r_from_reals

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

pure function get_omega_hat_from_ray(ray, N2, f2, Gamma2) result(omega_hat)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    type(t_ray), intent(in) :: ray
    real,        intent(in) :: N2, f2, Gamma2
    real                    :: omega_hat

    ! --------------------------------------------------------------------------
    
    omega_hat = get_omega_hat_from_reals(ray%k, ray%l, ray%m, N2, f2, Gamma2)

end function get_omega_hat_from_ray

pure function get_omega_hat_from_reals(k, l, m, N2, f2, Gamma2) result(omega_hat)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: k, l, m, N2, f2, Gamma2
    real             :: omega_hat

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: wvn_hor_sq, wvn_ver_sq

    ! --------------------------------------------------------------------------
    
    wvn_hor_sq = k ** 2 + l ** 2
    wvn_ver_sq = m ** 2 + Gamma2

    omega_hat = sqrt( &
        (N2 * wvn_hor_sq + f2 * wvn_ver_sq) / &
        (wvn_hor_sq + wvn_ver_sq) &
    )

end function get_omega_hat_from_reals

end module msgwam_rays_mod