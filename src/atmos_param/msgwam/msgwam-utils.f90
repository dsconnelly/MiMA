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
    integer :: dir

    ! --------------------------------------------------------------------------

    q = q_guess
    if (z(q) >= r .and. r > z(q + 1)) return

    if (r > z(1)) then
        q = 1
        return
    end if

    if (r < z(q_max)) then
        q = q_max - 1
        return
    end if

    dir = merge(-1, 1, r > z(q))

    do while (.true.)
        q = q + dir
        if (z(q) >= r .and. r > z(q + 1)) return
    end do

end function locate

pure subroutine shapiro_filter(profile)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:, :, :), intent(inout) :: profile

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: q_max
    real, dimension(4, size(profile, 2), size(profile, 3)) :: temp

    ! --------------------------------------------------------------------------

    q_max = size(profile, 1)
    temp(1:2, :, :) = profile(1:2, :, :)
    temp(3:4, :, :) = profile((q_max - 1):q_max, :, :)

    profile(2:(q_max - 1), :, :) = ( &
        profile(:(q_max - 2), :, :) + &
        2 * profile(2:(q_max - 1), :, :) + &
        profile(3:, :, :) &
    ) / 4.

    profile(1, :, :) = (3 * temp(1, :, :) + temp(2, :, :)) / 4.
    profile(q_max, :, :) = (temp(3, :, :) + 3 * temp(4, :, :)) / 4.

end subroutine shapiro_filter

end module msgwam_utils_mod