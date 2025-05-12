module msgwam_constants_mod

! ==============================================================================
! This module contains various constants that can be set at the start of the
! integration and then referenced by other modules. The MS-GWaM namelist is also
! read and exposed by this module. This module's members are public by default.
! ==============================================================================

use constants_mod, only: PI
use fms_mod,       only: check_nml_error, close_file, file_exist, mpp_pe, &
                         mpp_root_pe, open_namelist_file, stdlog

implicit none
public

! ==============================================================================
! namelist
! ==============================================================================

real    :: boundary_flux_ex   = 0.005
real    :: boundary_flux_tr   = 0.005
logical :: break_waves        = .true.
real    :: cp_max             = 50.
real    :: cp_width_ex        = 35.
real    :: cp_width_tr        = 35.
real    :: dr_max             = 10000.
real    :: dr_min             = 50.
real    :: dr_source          = 1000.
real    :: epsilon            = 0.
real    :: lat_tropics        = 15.
integer :: max_age            = 10 * 86400
real    :: min_flux           = 1.e-8
real    :: mu                 = 1.e-3
integer :: n_max              = 2500
integer :: n_source           = 48
integer :: n_sponge           = 3
real    :: source_dlat        = 5.
real    :: source_pressure    = 300.e2
real    :: T_hat_source       = 10. * 3600
logical :: use_shapiro_filter = .true.

! These namelist parameters are for debugging only.
integer               :: debug_mode = 0
logical               :: print_prune_diag = .false.
integer, dimension(5) :: track      = (/ 0, 1, 1, 1, 1/)

namelist / msgwam_nml / &
    boundary_flux_ex, boundary_flux_tr, break_waves, cp_max, cp_width_ex, &
    cp_width_tr, dr_max, dr_min, dr_source, epsilon, lat_tropics, max_age, &
    min_flux, mu, n_max, n_source, n_sponge, source_dlat, source_pressure, &
    T_hat_source, use_shapiro_filter, debug_mode, print_prune_diag, track

private :: msgwam_nml

! ==============================================================================
! other global constants, set at initialization
! ==============================================================================

integer :: i_max, j_max, q_max
real, dimension(:), allocatable :: f2
real :: min_N2, one_over_epsilon

contains

subroutine init_msgwam_constants(lon_bounds, lat_bounds, p_ref)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(:), intent(in) :: lon_bounds, lat_bounds, p_ref

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i_err, io, j, log_unit, nml_unit
    real :: lat

    ! --------------------------------------------------------------------------

    if (file_exist("input.nml")) then
        nml_unit = open_namelist_file()
        i_err = 1

        do while (i_err /= 0)
            read(nml_unit, nml=msgwam_nml, iostat=io)
            i_err = check_nml_error(io, "msgwam_nml")
        end do

        call close_file(nml_unit)
    end if

    log_unit = stdlog()
    if (mpp_pe() == mpp_root_pe()) then
        write (log_unit, nml=msgwam_nml)
    end if

    i_max = size(lon_bounds) - 1
    j_max = size(lat_bounds) - 1
    q_max = size(p_ref) - 1

    allocate(f2(j_max))

    do j = 1, j_max
        lat = 0.5 * (lat_bounds(j) + lat_bounds(j + 1))
        f2(j) = (2 * PI * sin(lat) / 86400.) ** 2
    end do

    min_N2 = (2 * PI / (2 * 3600)) ** 2

    if (epsilon > 0) then
        one_over_epsilon = 1. / epsilon
    else
        one_over_epsilon = 1.
    end if

end subroutine init_msgwam_constants

end module msgwam_constants_mod