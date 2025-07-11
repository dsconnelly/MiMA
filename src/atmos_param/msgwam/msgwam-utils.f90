module msgwam_utils_mod

! ==============================================================================
! Various utility functions for MS-GWaM.
! ==============================================================================

use msgwam_constants_mod, only: q_max

implicit none
private

public get_interp_coeffs, locate, shapiro_filter

contains

pure subroutine get_interp_coeffs(z, dz_inv, r, q, a, b)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max),     intent(in)  :: z
    real, dimension(q_max - 1), intent(in)  :: dz_inv
    real,                       intent(in)  :: r
    integer,                    intent(in)  :: q
    real,                       intent(out) :: a, b

    ! --------------------------------------------------------------------------

    a = (r - z(q + 1)) * dz_inv(q)

    a = merge(0., a, r < z(q_max))
    a = merge(1., a, r > z(1))

    ! a = min(1., max(0., a))
    b = 1. - a

end subroutine get_interp_coeffs

pure function locate(z, r, q_guess) result(q)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real,                   intent(in) :: r
    real, dimension(q_max), intent(in) :: z
    integer,                intent(in) :: q_guess
    integer                            :: q

    ! --------------------------------------------------------------------------

    q = q_guess

    do while (r <= z(q + 1) .and. q < q_max - 1)
        q = q + 1
    end do

    do while (r > z(q) .and. q > 1)
        q = q - 1
    end do

end function locate

pure subroutine shapiro_filter(profile)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max + 1), intent(inout) :: profile

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    real :: a, b

    ! --------------------------------------------------------------------------

    a = (3. * profile(1) + profile(2)) / 4.
    b = (profile(q_max) + 3. * profile(q_max + 1)) / 4.

    profile(2:q_max) = ( &
        profile(:(q_max - 1)) + &
        2. * profile(2:q_max) + &
        profile(3:(q_max + 1)) &
    ) / 4.

    profile(1) = a
    profile(q_max + 1) = b

end subroutine

end module msgwam_utils_mod