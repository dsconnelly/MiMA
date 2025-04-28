module msgwam_debug_mod

! ==============================================================================
! Some useful debugging tools for the MS-GWaM code.
! ==============================================================================

use fms_mod,              only: error_mesg, FATAL, mpp_pe

use msgwam_constants_mod, only: debug_mode, i_max, j_max, n_max, track
use msgwam_rays_mod,      only: t_ray

implicit none
private

public check_rays, track_ray

contains

subroutine check_rays(rays)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    type(t_ray), dimension(n_max, i_max, j_max), intent(in) :: rays

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, n

    ! --------------------------------------------------------------------------

    if (debug_mode /= 1) return

    do j = 1, j_max
        do i = 1, i_max
            do n = 1, n_max

                if (rays(n, i, j)%meta == -1) cycle

                associate(ray => rays(n, i, j))

                    if (nan_or_inf(ray%m)) then
                        write(*, *) mpp_pe(), n, i, j, ray%meta
                        call error_mesg("msgwam", "error in check_rays", FATAL)
                    end if

                end associate

            end do
        end do
    end do

end subroutine check_rays

subroutine track_ray(rays, loc)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    type(t_ray), dimension(n_max, i_max, j_max), intent(in) :: rays
    integer,                                     intent(in) :: loc

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: pe, n, i, j, meta

    ! --------------------------------------------------------------------------

    if (debug_mode /= 2) return
    
    pe = track(1)
    n = track(2)
    i = track(3)
    j = track(4)
    meta = track(5)
    
    if (mpp_pe() /= pe) return
    if (rays(n, i, j)%meta /= meta) return

    associate (ray => rays(n, i, j))

        if (.not. nan_or_inf(ray%m)) return

        write (*, *) "at location", loc
        write (*, *) "r_lo =", ray%r_lo
        write (*, *) "r_hi =", ray%r_hi
        write (*, *) "m =", ray%m
        write (*, *) "age =", ray%age
        call error_mesg("msgwam", "tracked ray error", FATAL)

    end associate

end subroutine track_ray

pure function nan_or_inf(v) result(res)

    ! --------------------------------------------------------------------------
    ! arguments and result
    ! --------------------------------------------------------------------------
    real, intent(in) :: v
    logical          :: res

    ! --------------------------------------------------------------------------

    res = v /= v
    res = res .or. v > huge(v)
    res = res .or. v < -huge(v)

end function nan_or_inf

end module msgwam_debug_mod
