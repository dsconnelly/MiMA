module msgwam_utils_mod

! ==============================================================================
! Various utility functions for MS-GWaM.
! ==============================================================================

implicit none
private

public get_interp_coeffs, interp, locate, shapiro_filter

contains

pure subroutine get_interp_coeffs(z, r, q, a, b)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:), intent(in) :: z
    real,               intent(in) :: r
    integer,            intent(in) :: q
    real,               intent(out) :: a, b

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: q_max
    real :: dz

    ! --------------------------------------------------------------------------

    q_max = size(z)

    if (r > z(1)) then
        a = 1.
        b = 0.

    else if (r < z(q_max)) then
        a = 0.
        b = 1.

    else
        dz = z(q) - z(q + 1)
        a = (r - z(q + 1)) / dz
        b = (z(q) - r) / dz

    end if

end subroutine get_interp_coeffs

pure function interp(z, profile, r, q) result(v)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:), intent(in) :: z, profile
    real,               intent(in) :: r
    integer,            intent(in) :: q
    real                           :: v

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: q_max
    real :: dz

    ! --------------------------------------------------------------------------

    q_max = size(z)

    if (r > z(1)) then
        v = profile(1)
    else if (r < z(q_max)) then
        v = profile(q_max)
    else
        dz = z(q) - z(q + 1)
        v = profile(q + 1) * (z(q) - r) / dz + profile(q) * (r - z(q + 1)) / dz
    end if

end function interp

pure function locate(r, z, q_guess) result(q)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real,               intent(in) :: r
    real, dimension(:), intent(in) :: z
    integer,            intent(in) :: q_guess
    integer                        :: q

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: q_max

    ! --------------------------------------------------------------------------

    q_max = size(z)

    if (r > z(1)) then
        q = 1
        return
    end if

    if (r < z(q_max)) then
        q = q_max - 1
        return
    end if

    q = q_guess
    do while (.true.)
        if (z(q) < r) then
            q = q - 1
        else if (z(q + 1) .ge. r) then
            q = q + 1
        else ! z(q) => r > z(q + 1)
            exit
        end if
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