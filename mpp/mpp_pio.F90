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
  use mpp_mod,  only : mpp_pe, mpp_root_pe, mpp_npes, get_mpp_comm, mpp_max
  use mpp_mod,  only : mpp_error, FATAL, WARNING, lowercase
  use pio,      only : PIO_init, pio_createfile, PIO_initdecomp
  use pio,      only : pio_createfile, pio_openfile, pio_file_is_open
  use pio,      only : File_desc_t, var_desc_t, iosystem_desc_t, IO_desc_t
  use pio,      only : pio_write, pio_clobber, pio_nowrite
  use pio,      only : PIO_put_att, pio_inquire, PIO_get_local_array_size
  use pio,      only : pio_iotype_netcdf, pio_iotype_pnetcdf
  use pio,      only : PIO_set_log_level
  use pio,      only : PIO_DOUBLE, PIO_REAL, PIO_INT
  use pio,      only : PIO_64BIT_OFFSET, PIO_64BIT_DATA
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

  ! below variables (maxl3, maxl4, lens3, lens4, ndlo) are used to match iodesc pointers
  integer, parameter :: maxl3 = 5 ! max number of lengths for dim=3 for iodesc bookkeeping
  integer, parameter :: maxl4 = 5 ! max number of lengths for dim=4 for iodesc bookkeeping
  integer, dimension(maxl3) :: lens3 = 0 ! array of dimension lengths for 3rd dimension
  integer, dimension(maxl4) :: lens4 = 0 ! array of dimension lengths for 4th dimension
  integer, parameter :: ndlo = 5 ! number of potential dimension layouts. see get_dlo function.

  ! IO descriptions
  type ioDesc_t
     type (IO_desc_t), pointer :: ptr => NULL()
  end type ioDesc_t
  type (ioDesc_t), dimension(npos,ndlo,maxl3,maxl4) :: ioDesc_i
  type (ioDesc_t), dimension(npos,ndlo,maxl3,maxl4) :: ioDesc_r
  type (ioDesc_t), dimension(npos,ndlo,maxl3,maxl4) :: ioDesc_d

  contains

  subroutine mpp_pio_init()
    ! local
    integer :: unit_begin, unit_end, unit_nml, io_status, unit
    logical :: opened
    integer :: pe, npes, localcomm
    integer :: numAggregator = 0 !TODO
    integer :: ierr

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

    !ierr = PIO_set_log_level(1)

    print *, "initialized PIO: ", localcomm

  end subroutine mpp_pio_init

  subroutine mpp_pio_stage_ioDesc(basetype_nf, domain, ioDesc, pos, ndim, time_axis_index, dlen)
    integer,          intent(in) :: basetype_nf
    type(domain2D),   intent(in) :: domain
    type (IO_desc_t), pointer    :: ioDesc
    integer,          intent(in) :: pos ! EAST=3, NORTH=5, CENTER=7, CORNER=8
    integer,          intent(in) :: ndim
    integer,          intent(in) :: time_axis_index
    integer, dimension(:), pointer :: dlen
    ! local
    type (IO_desc_t), pointer :: ioDesc_wrk
    integer :: ni, nj
    integer :: local_size
    integer :: global_size
    integer :: nig, njg, nk
    integer :: is, ie, js, je
    integer :: isd, ied, jsd, jed
    integer :: ism, iem, jsm, jem
    integer :: ioff, joff
    integer :: i,j,k,l,n
    integer, dimension(:), allocatable :: dof
    logical :: decomp_init_needed
    integer :: dlo ! dimension layout index

    ioDesc_wrk => null()

    if (associated(ioDesc)) return

    ! check if cell position is valid
    select case (pos)
      case(EAST)
      case(NORTH)
      case(CENTER)
      case(CORNER)
      case default
        print *, "field position: ", pos
        call mpp_error(FATAL,'mpp_pio_stage_ioDesc - Unknown field position')
    end select

    ! determine which ioDesc is the corresponding one
    call get_ioDesc_ptr(pos, ndim, dlen, basetype_nf, time_axis_index, ioDesc_wrk)

    decomp_init_needed = .false.
    if (.not. associated(ioDesc_wrk)) then
      allocate(ioDesc_wrk)
      decomp_init_needed = .true.
    endif

    ! If not already done, initialize the corresponding ioDesc
    if (decomp_init_needed) then

      call mpp_get_compute_domain( domain, is, ie, js, je, xsize = ni, ysize = nj, position=pos)
      call mpp_get_global_domain ( domain, xsize = nig, ysize = njg, position=pos)
      call mpp_get_data_domain   ( domain, isd, ied, jsd, jed, position=pos )
      call mpp_get_memory_domain ( domain, ism, iem, jsm, jem, position=pos )

      ! determine the index offsets:
      ioff = 1 - is
      joff = 1 - js
      call mpp_max(ioff)
      call mpp_max(joff)

      if (ndim == 2) then
        global_size = nig*njg
        local_size = ni*nj
        allocate(dof(local_size))

        n = 1
        do j=js+joff,je+joff
          do i=is+ioff,ie+ioff
            dof(n) = i + (j-1)*nig
            n = n+1
          enddo
        enddo

        ! sanity check:
        if (any(dof<1) .or. any(dof>global_size)) then
          print *, "global size: ", global_size, &
                   "minval(dof) ", minval(dof),  &
                   "maxval(dof) ", maxval(dof)
          call mpp_error(FATAL,'error in dof construction')
        endif

        ! initialize decomposition
        call PIO_initdecomp(pio_iosystem, basetype_nf, (/nig, njg/), dof, ioDesc_wrk)

      else if (ndim == 3) then
        ! construct dof
        global_size = nig*njg*dlen(3)
        local_size = ni*nj*dlen(3)
        allocate(dof(local_size))

        n = 1
        do k=1,dlen(3)
          do j=js+joff,je+joff
            do i=is+ioff,ie+ioff
              dof(n) = i + (j-1)*nig + (k-1)*nig*njg
              n = n+1
            enddo
          enddo
        enddo

        ! sanity check:
        if (any(dof<1) .or. any(dof>global_size)) then
          print *, "global size: ", global_size, &
                   "minval(dof) ", minval(dof),  &
                   "maxval(dof) ", maxval(dof)
          call mpp_error(FATAL,'error in dof construction')
        endif

        ! initialize decomposition
        call PIO_initdecomp(pio_iosystem, basetype_nf, (/nig, njg, dlen(3)/), dof, ioDesc_wrk)

      else if (ndim == 4) then
        ! construct dof
        global_size = nig*njg*dlen(3)*dlen(4)
        local_size = ni*nj*dlen(3)*dlen(4)
        allocate(dof(local_size))

        n = 1
        do l=1,dlen(4)
          do k=1,dlen(3)
            do j=js+joff,je+joff
              do i=is+ioff,ie+ioff
                dof(n) = i + (j-1)*nig + (k-1)*nig*njg
                n = n+1
              enddo
            enddo
          enddo
        enddo

        ! sanity check:
        if (any(dof<1) .or. any(dof>global_size)) then
          print *, "global size: ", global_size, &
                   "minval(dof) ", minval(dof),  &
                   "maxval(dof) ", maxval(dof)
          call mpp_error(FATAL,'error in dof construction')
        endif

        ! initialize decomposition
        call PIO_initdecomp(pio_iosystem, basetype_nf, (/nig, njg, dlen(3), dlen(4)/), dof, ioDesc_wrk)

      else
        call mpp_error(FATAL,'mpp_pio_decomp_init - Unsupported number of dimensions')
      endif
    else
      ! Check if the size of the decomposition to be adopted is compatible with the size of the field:
      ! to be adopted is the right size:
      call mpp_get_compute_domain( domain, is, ie, js, je, xsize = ni, ysize = nj, position=pos)
      select case (ndim)
        case(2)
          local_size = ni*nj
        case(3)
          local_size = ni*nj*dlen(3)
        case(4)
          local_size = ni*nj*dlen(3)*dlen(4)
        case default
          call mpp_error(FATAL,'Incompatible ndim in mpp_pio_stage_ioDesc.')
      end select

      if (local_size /= PIO_get_local_array_size(ioDesc_wrk))then
          print *, mpp_pe(), "local_size:", local_size, "decomp.size:", PIO_get_local_array_size(ioDesc_wrk)
          print *, mpp_pe(), "ni:",ni, "nj:",nj, "pos:",pos
          call mpp_error(FATAL,'incompatible local_size and ioDesc.')
      endif
    endif

    ioDesc => ioDesc_wrk

  end subroutine mpp_pio_stage_ioDesc

  function mpp_pio_openfile(file_desc, file_name, action_flag, ncid)
    type(File_desc_t),  intent(inout)         :: file_desc
    character(len=*),   intent(in)            :: file_name
    integer,            intent(in)            :: action_flag
    integer,            intent(out)           :: ncid
    ! local
    integer :: mpp_pio_openfile, stat
    integer :: nmode
    logical :: file_exists
    integer :: nDimensions

    print *, "OPENING >>>>>>>>>>>>>>>> ", trim(file_name)

    inquire(file=trim(file_name), exist=file_exists)

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

    if (trim(pio_netcdf_format) == "64bit_offset" .or. trim(pio_netcdf_format) == "64BIT_OFFSET") then
      nmode = ior(nmode, PIO_64BIT_OFFSET)
    else if (trim(pio_netcdf_format) == "64bit_data" .or. trim(pio_netcdf_format) == "64BIT_DATA") then
      nmode = ior(nmode, PIO_64BIT_DATA)
    else
      call mpp_error(FATAL,'Cannot determine pio_netcdf_format')
    endif

    if (pio_file_is_open(file_desc)) then
      mpp_pio_openfile = 0
      return
    endif

    if (action_flag == MPP_WRONLY .or. action_flag == MPP_OVERWR ) then
      stat = pio_createfile(pio_iosystem, file_desc, pio_iotype, trim(file_name), nmode)
    else if (action_flag == MPP_RDONLY) then
      stat = pio_openfile(pio_iosystem, file_desc, pio_iotype, trim(file_name), nmode)
    else
      print *, "NOT_IMPLEMENTED ", __FILE__, __LINE__
      call mpp_error(FATAL,'TODO - NOT_IMPLEMENTED')
    endif

    ncid = file_desc%fh

    mpp_pio_openfile = stat
  end function mpp_pio_openfile

  ! determine dimension layout index:
  ! ( we need to determine the dimension layout, because we may need to treat, for instance,  two 3D fields
  ! differently, e.g., when the third dimension of one field is z-coordinate  while the third dimension of
  ! another field may be time.
  integer function get_dlo(ndim, time_axis_index)
    integer, intent(in) :: ndim
    integer, intent(in) :: time_axis_index

    if (ndim==2) then
      get_dlo = 1
    elseif (ndim==3 .and. time_axis_index==-1) then
      get_dlo = 2
    elseif (ndim==3 .and. time_axis_index/=-1) then
      get_dlo = 3
    elseif (ndim==4 .and. time_axis_index==-1) then
      get_dlo = 4
    elseif (ndim==4 .and. time_axis_index/=-1) then
      get_dlo = 5
    else
      call mpp_error(FATAL,'get_dlo - Cannot determine dimension configuration')
    endif

    return
  end function get_dlo

  subroutine get_ioDesc_ptr(pos, ndim, dlen, basetype_nf, time_axis_index, ioDesc_ptr)
    integer, intent(in) :: pos
    integer, intent(in) :: ndim
    integer, dimension(:), pointer :: dlen
    integer, intent(in) :: basetype_nf
    integer, intent(in) :: time_axis_index
    type (IO_desc_t), pointer  :: ioDesc_ptr
    !local
    integer :: len3
    integer :: len4
    integer :: i, ilen3, ilen4
    integer :: dlo ! dimension layout index

    if (associated(ioDesc_ptr)) then
      call mpp_error(FATAL,'ioDesc to be set already associated')
    endif

    dlo = get_dlo(ndim, time_axis_index)

    len3=1; len4=1;
    if (ndim>=3) len3 = dlen(3)
    if (ndim>=4) len4 = dlen(4)

    ! Determine the dimension length index for dimension3
    ilen3 = 0
    if (len3==1) then
      ilen3 = 1
    else if (len3>1) then
      ! check if this size was entered before
      do i=2, maxl3
        if (lens3(i) == len3) ilen3 = i
      enddo

      ! ilen3 was not entered before. Enter it:
      if (ilen3==0) then
        if (.not. any(lens3(2:) == 0) ) then
          call mpp_error(FATAL,'No available entry left in lens3. May need to increase maxl3')
        endif
        do i=2,maxl3
          if (lens3(i) == 0 ) then
            lens3(i) = len3
            ilen3 = i
            exit
          endif
        enddo
      endif
    else
      call mpp_error(FATAL,'Invalid dimension length')
    endif

    ! Determine the dimension length index for dimension4
    ilen4 = 0
    if (len4==1) then
      ilen4 = 1
    else if (len4>1) then
      ! check if this size was entered before
      do i=2, maxl4
        if (lens4(i) == len4) ilen4 = i
      enddo

      ! ilen4 was not entered before. Enter it:
      if (ilen4==0) then
        if (.not. any(lens4(2:) == 0) ) then
          call mpp_error(FATAL,'No available entry left in lens4. May need to increase maxl4')
        endif
        do i=2,maxl4
          if (lens4(i) == 0 ) then
            lens4(i) = len4
            ilen4 = i
            exit
          endif
        enddo
      endif
    else
      call mpp_error(FATAL,'Invalid dimension length')
    endif

    select case (basetype_nf)
      case(NF_INT)
        ioDesc_ptr => ioDesc_i(pos,dlo,ilen3,ilen4)%ptr
      case(NF_REAL)
        ioDesc_ptr => ioDesc_r(pos,dlo,ilen3,ilen4)%ptr
      case(NF_DOUBLE)
        ioDesc_ptr => ioDesc_d(pos,dlo,ilen3,ilen4)%ptr
      case default
        call mpp_error(FATAL,'get_ioDesc - Unknown kind')
    end select

  end subroutine get_ioDesc_ptr

#endif
end module mpp_pio_mod
