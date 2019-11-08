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

! PIO only works with MPI:
#ifndef use_libMPI
#undef use_PIO
#endif

#ifdef use_PIO

#include <fms_platform.h>

  use mpp_mod,  only : input_nml_file
  use mpp_mod,  only : mpp_pe, mpp_root_pe, mpp_npes, get_mpp_comm
  use pio,      only : PIO_init
  use pio,      only : iosystem_desc_t

  implicit none
  private

  public :: mpp_pio_init

  character(len=64)   :: pio_netcdf_format, pio_typename
  integer             :: pio_numiotasks, pio_rearranger, pio_root, pio_stride
  namelist /mpp_pio_nml/ pio_netcdf_format, pio_numiotasks, pio_rearranger, &
                         pio_root, pio_stride, pio_typename

  type(iosystem_desc_t) :: pio_iosystem     ! The ParallelIO system set up by PIO_init
  integer               :: pio_optbase = 1  ! Start index of I/O processors

  contains

  subroutine mpp_pio_init()
    integer :: unit_begin, unit_end, unit_nml, io_status, unit
    logical :: opened
    integer :: pe, npes, localcomm
    integer :: numAggregator = 0 !TODO

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

#endif
end module mpp_pio_mod
