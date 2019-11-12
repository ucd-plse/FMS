!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!-----------------------------------------------------------------------
!                 ParallelIO (PIO) Module for FMS
!
! AUTHOR: A. Altuntas (altuntas@ucar.edu)
!         National Center for Atmospheric Research
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! For the full text of the GNU General Public License,
! write to: Free Software Foundation, Inc.,
!           675 Mass Ave, Cambridge, MA 02139, USA.
!-----------------------------------------------------------------------

! <CONTACT EMAIL="altuntas@ucar.edu">
!   A. Altuntas
! </CONTACT>

! <OVERVIEW>
!   <TT>mpp_pio_mod</TT>, is ... TODO
! </OVERVIEW>

! <DESCRIPTION>
!   TODO
! </DESCRIPTION>

module mpp_pio_mod

! PIO requires MPI:
#ifndef use_libMPI
#undef use_PIO
#endif

#ifdef use_PIO

#include <fms_platform.h>

  use mpp_mod,  only : input_nml_file
  use mpp_mod,  only : mpp_pe, mpp_root_pe, mpp_npes, get_mpp_comm
  use mpp_mod,  only : mpp_error, FATAL, WARNING, lowercase
  use pio,      only : PIO_init, pio_createfile
  use pio,      only : pio_createfile, pio_openfile, pio_file_is_open
  use pio,      only : File_desc_t, var_desc_t, iosystem_desc_t
  use pio,      only : pio_write, pio_clobber, pio_nowrite
  use pio,      only : PIO_put_att
  use pio,      only : pio_iotype_netcdf, pio_iotype_pnetcdf
  use mpp_parameter_mod,  only : MPP_WRONLY, MPP_RDONLY, MPP_APPEND, MPP_OVERWR

  implicit none
  private

  public :: mpp_pio_init
  public :: mpp_pio_openfile
  public :: write_attribute_pio

  character(len=64)   :: pio_netcdf_format, pio_typename
  integer             :: pio_numiotasks, pio_rearranger, pio_root, pio_stride
  namelist /mpp_pio_nml/ pio_netcdf_format, pio_numiotasks, pio_rearranger, &
                         pio_root, pio_stride, pio_typename

  type(iosystem_desc_t) :: pio_iosystem     ! The ParallelIO system set up by PIO_init
  integer               :: pio_iotype       ! PIO_IOTYPE_NETCDF or PNETCDF
  integer               :: pio_optbase = 1  ! Start index of I/O processors
  logical               :: pio_cf_compliance

  contains

  subroutine mpp_pio_init(cf_compliance)
    logical, intent(in) :: cf_compliance
    ! local
    integer :: unit_begin, unit_end, unit_nml, io_status, unit
    logical :: opened
    integer :: pe, npes, localcomm
    integer :: numAggregator = 0 !TODO

    pio_cf_compliance = cf_compliance

    ! read mpp_pio_nml namelist
#ifdef INTERNAL_FILE_NML
      read (input_nml_file, mpp_pio_nml, iostat=io_status)
#else
      unit_begin = 103
      unit_end   = 512
      do unit_nml = unit_begin, unit_end
         inquire(unit_nml, OPENED=opened )
         if( .NOT.opened )exit
      end do

      open(unit_nml,file='input.nml', iostat=io_status)
      read(unit_nml,mpp_pio_nml,iostat=io_status)
      close(unit_nml)
#endif

    if (trim(pio_typename) == "netcdf") then
      pio_iotype = pio_iotype_netcdf
    elseif (trim(pio_typename) == "pnetcdf") then
      pio_iotype = pio_iotype_pnetcdf
    else
      call mpp_error(FATAL,'Unknown PIO filetype')
    endif

    pe = mpp_pe()
    npes = mpp_npes()
    localcomm = get_mpp_comm()

    ! TODO
    ! add checks regarding the consistency of npes, numiotasksm stride, etc.

    !initialize PIO
    call PIO_init(    &
      pe,             & ! MPI rank
      localcomm,      & ! MPI communicator
      pio_numiotasks, & ! Number of iotasks
      numAggregator,  & ! number of aggregators to use
      pio_stride,     & ! stride
      pio_rearranger, & ! form of rearrangement
      pio_iosystem,   & ! The ParallelIO system set up by PIO_init
      pio_optbase )     ! Start index of I/O processors (optional)

    print *, "initialized PIO: ", localcomm

  end subroutine mpp_pio_init

  function mpp_pio_openfile(file_desc, file_name, action_flag)
    type(File_desc_t),  intent(inout)         :: file_desc
    character(len=*),   intent(in)            :: file_name
    integer,            intent(in)            :: action_flag
    ! local
    integer :: mpp_pio_openfile, stat
    integer :: nmode
    logical :: file_exists
    character(len=128) :: file_name_debug ! TODO remove this

    print *, "OPENING >>>>>>>>>>>>>>>> ", trim(file_name)

    !TODO Remove below if statement:
    if (trim(file_name) == "./ocean_geometry.nc") then
      file_name_debug = "./ocean_geometry_pio.nc"
    else if (trim(file_name) == "./Vertical_coordinate.nc") then
      file_name_debug = "./Vertical_coordinate_pio.nc"
    else if (trim(file_name) == "./MOM_IC.nc") then
      file_name_debug = "./MOM_IC_pio.nc"
    else if (trim(file_name) == "./ocean.stats.nc") then
      file_name_debug = "./ocean.stats_pio.nc"
    else
      print *, "debug ", __FILE__, __LINE__
      call mpp_error(FATAL,'TODO - NOT_IMPLEMENTED')
    endif

    inquire(file=file_name_debug, exist=file_exists)

    if (action_flag == MPP_WRONLY) then
      nmode = pio_write
    else if (action_flag == MPP_OVERWR) then
      if (.not. file_exists) then
        print *, "NOT_IMPLEMENTED ", __FILE__, __LINE__
        call mpp_error(FATAL,'TODO - NOT_IMPLEMENTED')
      endif
      nmode = pio_clobber
    else if (action_flag == MPP_RDONLY) then
      nmode = pio_nowrite
    else
      print *, "NOT_IMPLEMENTED ", __FILE__, __LINE__
      call mpp_error(FATAL,'TODO - NOT_IMPLEMENTED')
    endif

    if (pio_file_is_open(file_desc)) then
      mpp_pio_openfile = 0
      return
    endif

    if (action_flag == MPP_WRONLY .or. action_flag == MPP_OVERWR ) then
      stat = pio_createfile(pio_iosystem, file_desc, pio_iotype, trim(file_name_debug), nmode)
    else if (action_flag == MPP_RDONLY) then
      stat = pio_openfile(pio_iosystem, file_desc, pio_iotype, trim(file_name_debug), nmode)
    else
      print *, "NOT_IMPLEMENTED ", __FILE__, __LINE__
      call mpp_error(FATAL,'TODO - NOT_IMPLEMENTED')
    endif

    mpp_pio_openfile = stat
  end function mpp_pio_openfile

  subroutine write_attribute_pio(file_desc, file_name, attr_id, attr_name, rval, ival, cval, pack)
    type(File_desc_t),intent(inout) :: file_desc
    character(len=*), intent(in)    :: file_name
    type(var_desc_t), intent(in)    :: attr_id
    character(len=*), intent(in)    :: attr_name
    real,             intent(in), optional :: rval(:)
    integer,          intent(in), optional :: ival(:)
    character(len=*), intent(in), optional :: cval
    integer,          intent(in), optional :: pack
    ! local:
    integer, allocatable :: rval_i(:)
    integer :: error ! error code

    ! this subroutine body is copied from mpp_io_write.inc::write_attribute_netcdf
    ! and is modified for PIO

    if( PRESENT(rval) )then ! attributes of type real
      if( PRESENT(pack) )then
        if( pack== 0 ) then !! here be dragons, use ival branch!...
          call mpp_error( FATAL, &
                  'write_attribute_pio: attempting to write internal NF_INT, currently int32, as real.' )
        else if( pack.EQ.1 )then
          if( KIND(rval).EQ.DOUBLE_KIND )then
            error = PIO_put_att(file_desc, attr_id, attr_name, rval)
          else if( KIND(rval).EQ.FLOAT_KIND )then
            call mpp_error( WARNING, &
                 'write_attribute_pio: attempting to write internal 32-bit real as external 64-bit.' )
            error = PIO_put_att(file_desc, attr_id, attr_name, rval)
          end if
          call pio_err( error, file_name=file_name, attr_name=attr_name)
        else if( pack.EQ.2 )then
          error = PIO_put_att(file_desc, attr_id, attr_name, rval)
          call pio_err( error, file_name=file_name, attr_name=attr_name)
        else
          call mpp_error( FATAL, 'write_attribute_pio: only legal packing values are 1,2.' )
        end if
      else ! no pack provided
        error = PIO_put_att(file_desc, attr_id, attr_name, rval)
        call pio_err( error, file_name=file_name, attr_name=attr_name)
      end if
    else if( PRESENT(ival) )then ! attributes of type integer
      if( PRESENT(pack) ) then
        if (pack ==0) then
          if (KIND(ival).EQ.LONG_KIND ) then
             call mpp_error(FATAL,'only use NF_INTs with pack=0 for now')
          end if
          error = PIO_put_att(file_desc, attr_id, attr_name, ival)
          call pio_err( error, file_name=file_name, attr_name=attr_name)
        else
          call mpp_error( FATAL, 'write_attribute_pio: only implimented ints when pack=0, else use reals.' )
        endif
      else ! no pack provided
        error = PIO_put_att(file_desc, attr_id, attr_name, ival)
        call pio_err( error, file_name=file_name, attr_name=attr_name)
      end if
    else if( present(cval) )then
      if (.NOT.pio_cf_compliance .or. trim(attr_name).NE.'calendar') then
        error = PIO_put_att(file_desc, attr_id, attr_name, cval)
      else
        error = PIO_put_att(file_desc, attr_id, attr_name, lowercase(cval))
      endif
      call pio_err( error, file_name=file_name, attr_name=attr_name)
    else
      call mpp_error( FATAL, 'write_attribute_pio: one of rval, ival, cval must be present.' )
    end if

  end subroutine write_attribute_pio

  subroutine pio_err( err, file_name, axis_name, field_name, attr_name, string)
    integer,          intent(in)           :: err
    character(len=*), intent(in), optional :: file_name
    character(len=*), intent(in), optional :: axis_name
    character(len=*), intent(in), optional :: field_name
    character(len=*), intent(in), optional :: attr_name
    character(len=*),             optional :: string
    ! local
    character(len=256) :: errmsg

    errmsg = ""

    if (err/=0) then ! TODO confirm this
      if (present(file_name))   errmsg = trim(errmsg)//' File='//file_name
      if (present(axis_name))   errmsg = trim(errmsg)//' Axis='//axis_name
      if (present(field_name))  errmsg = trim(errmsg)//' Field='//field_name
      if (present(attr_name))   errmsg = trim(errmsg)//' Attribute='//attr_name
      if (present(string))      errmsg = trim(errmsg)//string

      !call mpp_io_exit('NOSYNC')        !make sure you close all open files !TODO
      call mpp_error( FATAL, 'PIO ERROR: '//trim(errmsg) )
    endif
  end subroutine pio_err

#endif
end module mpp_pio_mod
