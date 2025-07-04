module msgwam_io_mod

! ==============================================================================
! This module implements subroutines for handling netCDF outputs and for reading
! and writing ray volume state for restarts.
! ==============================================================================

use diag_manager_mod,     only: register_diag_field, send_data
use fms_mod,              only: error_mesg, FATAL, mpp_pe
use time_manager_mod,     only: time_type

use msgwam_constants_mod, only: i_max, j_max, n_max, n_source, q_max, &
                                steady_state
use msgwam_rays_mod,      only: t_ray

implicit none
private

public init_nc_output, init_ray_state, save_ray_state, send_nc_output

! ==============================================================================
! netCDF input/ouput variables
! ==============================================================================

character(len=7) :: mod_name = "msgwam"
real, parameter :: missing_value = -999.
integer :: id_u, id_v, id_rho, id_N2, id_G2, id_flux_x, id_flux_y, &
    id_accel_x, id_accel_y

! ==============================================================================

contains

subroutine init_nc_output(axes, Time)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer, dimension(4), intent(in) :: axes
    type(time_type),       intent(in) :: Time

    ! --------------------------------------------------------------------------

    id_u = register_diag_field(mod_name, "gw_u", axes(1:3), Time, &
        "MS-GWaM zonal wind", "m / s", missing_value=missing_value)
    id_v = register_diag_field(mod_name, "gw_v", axes(1:3), Time, &
        "MS-GWaM meridional wind", "m / s", missing_value=missing_value)

    id_rho = register_diag_field(mod_name, "gw_rho", axes(1:3), Time, &
        "MS-GWaM density", "kg/m^3", missing_value=missing_value)
    id_N2 = register_diag_field(mod_name, "gw_N2", axes(1:3), Time, &
        "MS-GWaM squared BV frequency", "1/s^2", missing_value=missing_value)
    id_G2 = register_diag_field(mod_name, "gw_G2", axes(1:3), Time, &
        "MS-GWaM squared scale height", "1/m^2", missing_value=missing_value)

    id_flux_x = register_diag_field(mod_name, "gw_flux_x", axes(1:3), Time, &
        "zonal GW flux", "Pa", missing_value=missing_value)
    id_flux_y = register_diag_field(mod_name, "gw_flux_y", axes(1:3), Time, &
        "meridional GW flux", "Pa", missing_value=missing_value)

    id_accel_x = register_diag_field(mod_name, "gw_accel_x", axes(1:3), Time, &
        "zonal GW acceleration", "m/s^2", missing_value=missing_value)
    id_accel_y = register_diag_field(mod_name, "gw_accel_y", axes(1:3), Time, &
        "meridional GW acceleration", "m/s^2", missing_value=missing_value)

end subroutine init_nc_output

subroutine init_ray_state(rays, ghosts, last_meta)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    type(t_ray), dimension(n_max, i_max, j_max), intent(out) :: rays
    integer, dimension(n_source, i_max, j_max),  intent(out) :: ghosts
    integer, dimension(i_max, j_max),            intent(out) :: last_meta

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: iostat, unit
    character(len=32) :: fname, pe_str
    logical :: from_restart

    ! --------------------------------------------------------------------------

    write(pe_str, "(I2.2)") mpp_pe()
    fname = "INPUT/rays-" // trim(pe_str) // ".dat"
    inquire(file=trim(fname), exist=from_restart)

    if (steady_state .or. (.not. from_restart)) then
        rays(:, :, :)%meta = -1
        ghosts(:, :, :) = -1
        last_meta(:, :) = 1

        return
    end if

    open(newunit=unit, file=trim(fname), form="unformatted", &
        iostat=iostat, action="read")

    if (iostat /= 0) then
        call error_mesg("msgwam_io_mod", "error loading ray state", FATAL)
    end if

    read(unit) rays
    read(unit) ghosts
    read(unit) last_meta

    close(unit)    

end subroutine init_ray_state

pure subroutine reorder_axes(a, b)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    real, dimension(q_max, i_max, j_max), intent(in)  :: a
    real, dimension(i_max, j_max, q_max), intent(out) :: b

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i, j, q

    ! --------------------------------------------------------------------------

    do j = 1, j_max
        do i = 1, i_max
            do q = 1, q_max
                b(i, j, q) = a(q, i, j)
            end do
        end do
    end do

end subroutine reorder_axes

subroutine save_ray_state(rays, ghosts, last_meta)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    type(t_ray), dimension(n_max, i_max, j_max), intent(in) :: rays
    integer, dimension(n_source, i_max, j_max),  intent(in) :: ghosts
    integer, dimension(i_max, j_max),            intent(in) :: last_meta

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: iostat, unit
    character(len=32) :: fname, pe_str

    ! --------------------------------------------------------------------------

    if (steady_state) then
        return
    end if

    write(pe_str, "(I2.2)") mpp_pe()
    fname = "RESTART/rays-" // trim(pe_str) // ".dat"

    open(newunit=unit, file=fname, form="unformatted", &
        iostat=iostat, action="write")

    if (iostat /= 0) then
        call error_mesg("msgwam_io_mod", "error saving ray state", FATAL)
    end if

    write(unit) rays
    write(unit) ghosts
    write(unit) last_meta

    close(unit)

end subroutine save_ray_state

subroutine send_nc_output(i_start, j_start, Time, u_bar, v_bar, rho, N2, G2, &
    flux_x, flux_y, du_dt, dv_dt)

    ! --------------------------------------------------------------------------
    ! arguments
    ! --------------------------------------------------------------------------
    integer,                                  intent(in) :: i_start, j_start
    type(time_type),                          intent(in) :: Time
    real, dimension(q_max, i_max, j_max),     intent(in) :: u_bar, v_bar, N2, &
                                                            G2, rho
    real, dimension(q_max + 1, i_max, j_max), intent(in) :: flux_x, flux_y
    real, dimension(i_max, j_max, q_max),     intent(in) :: du_dt, dv_dt

    ! --------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------
    integer :: i_err
    real, dimension(i_max, j_max, q_max) :: temp

    ! --------------------------------------------------------------------------

    if (id_u > 0) then
        call reorder_axes(u_bar, temp)
        i_err = send_data(id_u, temp, Time, i_start, j_start)
    end if

    if (id_v > 0) then
        call reorder_axes(v_bar, temp)
        i_err = send_data(id_v, temp, Time, i_start, j_start)
    end if

    if (id_rho > 0) then
        call reorder_axes(rho, temp)
        i_err = send_data(id_rho, temp, Time, i_start, j_start)
    end if

    if (id_N2 > 0) then
        call reorder_axes(N2, temp)
        i_err = send_data(id_N2, temp, Time, i_start, j_start)
    end if

    if (id_G2 > 0) then
        call reorder_axes(G2, temp)
        i_err = send_data(id_G2, temp, Time, i_start, j_start)
    end if

    if (id_flux_x > 0) then
        call reorder_axes(flux_x(:q_max, :, :), temp)
        i_err = send_data(id_flux_x, temp, Time, i_start, j_start)
    end if

    if (id_flux_y > 0) then
        call reorder_axes(flux_y(:q_max, :, :), temp)
        i_err = send_data(id_flux_y, temp, Time, i_start, j_start)
    end if

    if (id_accel_x > 0) then
        i_err = send_data(id_accel_x, du_dt, Time, i_start, j_start)
    end if

    if (id_accel_y > 0) then
        i_err = send_data(id_accel_y, dv_dt, Time, i_start, j_start)
    end if

end subroutine send_nc_output

end module msgwam_io_mod