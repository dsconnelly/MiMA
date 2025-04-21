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

real    :: boundary_flux      = 0.01
logical :: break_waves        = .true.
real    :: cp_max             = 50.
real    :: cp_width           = 35.
real    :: dr_min             = 50.
real    :: dr_source          = 1000.
real    :: epsilon            = 0.
real    :: lat_extrinsic      = 15.
integer :: max_age            = 10 * 86400
real    :: min_flux           = 1.e-8
real    :: min_N2             = 0.005 ** 2
real    :: mu                 = 1.e-3
integer :: n_max              = 2500
integer :: n_source           = 48
real    :: source_pressure    = 300.e2
real    :: T_hat_source       = 10. * 3600
logical :: use_shapiro_filter = .true.

! To be removed, eventually.
real    :: H_rho              = 8.e+3

namelist / msgwam_nml / &
    boundary_flux, break_waves, cp_max, cp_width, dr_min, dr_source, epsilon, &
    lat_extrinsic, max_age, min_flux, min_N2, mu, n_max, n_source, &
    source_pressure, T_hat_source, use_shapiro_filter, H_rho

private :: lat_extrinsic, msgwam_nml

! ==============================================================================
! other global constants, set at initialization
! ==============================================================================

integer :: i_max, j_max, q_max
real, dimension(:), allocatable :: f2
logical, dimension(:), allocatable :: is_extrinsic
real :: Gamma2

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
    allocate(is_extrinsic(j_max))

    do j = 1, j_max
        lat = 0.5 * (lat_bounds(j) + lat_bounds(j + 1))
        is_extrinsic(j) = abs(180 * lat / PI) > lat_extrinsic
        f2(j) = (2 * PI * sin(lat) / 86400.) ** 2
    end do

    Gamma2 = ((1. / 2. - 2. / 7.) / H_rho) ** 2
        
end subroutine init_msgwam_constants

end module msgwam_constants_mod