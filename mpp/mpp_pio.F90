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
  use pio,      only : PIO_init, pio_createfile, PIO_initdecomp
  use pio,      only : pio_createfile, pio_openfile, pio_file_is_open
  use pio,      only : File_desc_t, var_desc_t, iosystem_desc_t, IO_desc_t
  use pio,      only : pio_write, pio_clobber, pio_nowrite
  use pio,      only : PIO_put_att
  use pio,      only : pio_iotype_netcdf, pio_iotype_pnetcdf
  use pio,      only : PIO_DOUBLE, PIO_REAL, PIO_INT
  use mpp_parameter_mod,  only : MPP_WRONLY, MPP_RDONLY, MPP_APPEND, MPP_OVERWR
  use mpp_parameter_mod,  only : CENTER, EAST, NORTH, CORNER
  use mpp_domains_mod,    only : domain2d, mpp_get_compute_domain, mpp_get_global_domain
  use mpp_domains_mod,    only : mpp_get_data_domain, mpp_get_memory_domain

  implicit none
  private

#include <netcdf.inc>

  public :: mpp_pio_init
  public :: mpp_pio_stage_ioDesc
  public :: mpp_pio_openfile

  character(len=64)   :: pio_netcdf_format, pio_typename
  integer             :: pio_numiotasks, pio_rearranger, pio_root, pio_stride
  namelist /mpp_pio_nml/ pio_netcdf_format, pio_numiotasks, pio_rearranger, &
                         pio_root, pio_stride, pio_typename

  type(iosystem_desc_t) :: pio_iosystem     ! The ParallelIO system set up by PIO_init
  integer               :: pio_iotype       ! PIO_IOTYPE_NETCDF or PNETCDF
  integer               :: pio_optbase = 1  ! Start index of I/O processors

  integer, parameter :: npos = 8 ! EAST=3, NORTH=5, CENTER=7, CORNER=8
  integer, parameter :: ndimc = 5 ! 2: 2d in space
                                  ! 3: 3d in space
                                  ! 4: 3d in space, 1d in time
                                  ! 5: 2d in space, 1d in time

  ! IO descriptions
  type ioDesc_t
     type (IO_desc_t), pointer :: ptr => NULL()
  end type ioDesc_t
  type (ioDesc_t), dimension(npos,ndimc) :: ioDesc_i
  type (ioDesc_t), dimension(npos,ndimc) :: ioDesc_r
  type (ioDesc_t), dimension(npos,ndimc) :: ioDesc_d


  contains

  subroutine mpp_pio_init()
    ! local
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

  subroutine mpp_pio_stage_ioDesc(basetype_nf, domain, ioDesc,pos, ndim3, ndim4)
    integer,          intent(in) :: basetype_nf
    type(domain2D),   intent(in) :: domain
    type (IO_desc_t), pointer    :: ioDesc
    integer,          intent(in) :: pos ! EAST=3, NORTH=5, CENTER=7, CORNER=8
    integer,          intent(in) :: ndim3 ! size of third dimension
    integer,          intent(in), optional :: ndim4 ! size of fourth dimension
    ! local
    type (IO_desc_t), pointer :: ioDesc_wrk
    integer :: basetype_pio
    integer :: ni, nj
    integer :: local_size
    integer :: global_size
    integer :: nig, njg, nk
    integer :: is, ie, js, je
    integer :: isd, ied, jsd, jed
    integer :: ism, iem, jsm, jem
    integer :: i,j,k,n
    integer, dimension(:), allocatable :: dof3d
    logical :: decomp_init_needed
    integer :: dimc

    if (associated(ioDesc)) return

    ! check if cell position is valid
    select case (pos)
      case(EAST)
      case(NORTH)
      case(CENTER)
      case(CORNER)
      case default
        call mpp_error(FATAL,'mpp_pio_stage_ioDesc - Unknown cell position')
    end select

    ! now determine dimension configuration index:
    dimc = -1
    if (.not. present(ndim4)) then
      if (ndim3==1) then
        dimc = 2
      else
        dimc = 3
      endif
    else
      if (ndim3>1 .and. ndim4==1) then
        dimc = 3
      else if (ndim3>1 .and. ndim4>1) then
        dimc = 4
      else if (ndim3==1 .and. ndim4>1) then
        dimc = 5
      endif
    endif
    if (dimc == -1) then
      call mpp_error(FATAL,'mpp_pio_stage_ioDesc - Cannot determine dimension configuration')
    endif

    ! determine which ioDesc is the corresponding one
    decomp_init_needed = .false.
    select case (basetype_nf)
      case(NF_INT)
        basetype_pio = PIO_INT
        if (.not. associated(ioDesc_i(pos, dimc)%ptr)) then
          allocate(ioDesc_i(pos, dimc)%ptr)
          decomp_init_needed = .true.
        endif
        ioDesc_wrk => ioDesc_i(pos,dimc)%ptr
      case(NF_REAL)
        basetype_pio = PIO_REAL
        if (.not. associated(ioDesc_r(pos, dimc)%ptr)) then
          allocate(ioDesc_r(pos,dimc)%ptr)
          decomp_init_needed = .true.
        endif
        ioDesc_wrk => ioDesc_r(pos,dimc)%ptr
      case(NF_DOUBLE)
        basetype_pio = PIO_DOUBLE
        if (.not. associated(ioDesc_d(pos, dimc)%ptr)) then
          allocate(ioDesc_d(pos,dimc)%ptr)
          decomp_init_needed = .true.
        endif
        ioDesc_wrk => ioDesc_d(pos,dimc)%ptr
      case default
        call mpp_error(FATAL,'mpp_pio_stage_ioDesc - Unknown kind')
    end select

    ! If not already done, initialize the corresponding ioDesc
    if (decomp_init_needed) then

      if (present(ndim4)) then
        call mpp_error(FATAL,'mpp_pio_decomp_init - 4th dimension not supported yet.')

      else ! 3 dimensions only
        ! construct dof3d
        call mpp_get_compute_domain( domain, is, ie, js, je, xsize = ni, ysize = nj, position=pos)
        call mpp_get_global_domain ( domain, xsize = nig, ysize = njg, position=pos)
        call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, position=pos )
        call mpp_get_memory_domain ( domain, ism, iem, jsm, jem, position=pos )

        global_size = nig*njg*ndim3
        local_size = ni*nj*ndim3
        allocate(dof3d(local_size))

        n = 1
        do k=1,ndim3
          do j=js,je
            do i=is,ie
              dof3d(n) = i + (j-1)*nig + (k-1)*nig*njg
              n = n+1
            enddo
          enddo
        enddo

        ! sanity check:
        if (any(dof3d<1) .or. any(dof3d>global_size)) then
          call mpp_error(FATAL,'error in dof3d construction')
        endif

        call mpp_pio_decomp_init(basetype_pio, domain, pos, dof3d, ioDesc_wrk, nig, njg, ndim3)
      endif
    endif

    ioDesc => ioDesc_wrk

  end subroutine mpp_pio_stage_ioDesc

  subroutine mpp_pio_decomp_init(basetype_pio, domain, pos, compdof, ioDesc, nig, njg, ndim3, ndim4)
    integer,                intent(in) :: basetype_pio
    type(domain2D),         intent(in) :: domain
    integer,                intent(in) :: pos ! EAST=3, NORTH=5, CENTER=7, CORNER=8
    integer, dimension(:),  intent(in) :: compdof(:)
    type (IO_desc_t),       pointer    :: ioDesc
    integer,                intent(in) :: nig, njg, ndim3
    integer, intent(in), optional      :: ndim4
    ! local

    if (.not. associated(ioDesc)) then
        call mpp_error(FATAL,'ioDesc must be allocated before initialization.')
    endif

    if (present(ndim4)) then
      call PIO_initdecomp(pio_iosystem, basetype_pio, (/nig, njg, ndim3, ndim4/), compdof, ioDesc)
    else
      call PIO_initdecomp(pio_iosystem, basetype_pio, (/nig, njg, ndim3/), compdof, ioDesc)
    endif

  end subroutine mpp_pio_decomp_init

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
        nmode = pio_write
      else
        nmode = pio_clobber
      endif
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

#endif
end module mpp_pio_mod
