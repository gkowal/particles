!!******************************************************************************
!!
!!  This file is part of the PARTICLES source code, a program to integrate
!!  test particles in magnetohydrodynamical simulations.
!!
!!  Copyright (C) 2018 Grzegorz Kowal <grzegorz@amuncode.org>
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!******************************************************************************
!!
!! module: FITSIO
!!
!!  This module handles FITS files.
!!
!!******************************************************************************
!
module fitsio

  implicit none

  interface write_data
    module procedure write_data_integer_1d
    module procedure write_data_integer_2d
    module procedure write_data_real_1d
    module procedure write_data_real_2d
    module procedure write_data_real_3d
  end interface

  character(len=128), save :: datadir  = './'
  character(len= 32), save :: uxname   = 'velx'
  character(len= 32), save :: uyname   = 'vely'
  character(len= 32), save :: uzname   = 'velz'
  character(len= 32), save :: bxname   = 'magx'
  character(len= 32), save :: byname   = 'magy'
  character(len= 32), save :: bzname   = 'magz'
  character(len=128), save :: uxfname  = './velx.fits'
  character(len=128), save :: uyfname  = './vely.fits'
  character(len=128), save :: uzfname  = './velz.fits'
  character(len=128), save :: bxfname  = './magx.fits'
  character(len=128), save :: byfname  = './magy.fits'
  character(len=128), save :: bzfname  = './magz.fits'

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_FITSIO:
! ----------------------------
!
!   Subroutine initialized module.
!
!
!===============================================================================
!
  subroutine initialize_fitsio()

    use parameters, only : get_parameter_string

    implicit none

    logical           :: info
!
!-------------------------------------------------------------------------------
!
    call get_parameter_string("datadir"        , datadir)
    call get_parameter_string("x-velocity-name", uxname)
    call get_parameter_string("y-velocity-name", uyname)
    call get_parameter_string("z-velocity-name", uzname)
    call get_parameter_string("x-magnetic-name", bxname)
    call get_parameter_string("y-magnetic-name", byname)
    call get_parameter_string("z-magnetic-name", bzname)

! add bar to the data path, if necessary
!
    if (index(datadir, '/', .true.) /= len_trim(datadir)) then
      write(datadir,"(a,'/')") trim(datadir)
    end if

! check if files exist
!
    call check_if_exists(datadir, uxname, uxfname, info)
    call check_if_exists(datadir, uyname, uyfname, info)
    call check_if_exists(datadir, uzname, uzfname, info)
    call check_if_exists(datadir, bxname, bxfname, info)
    call check_if_exists(datadir, byname, byfname, info)
    call check_if_exists(datadir, bzname, bzfname, info)

!-------------------------------------------------------------------------------
!
  end subroutine initialize_fitsio
!
!===============================================================================
!
! subroutine FINALIZE_FITSIO:
! --------------------------
!
!   Subroutine finalizes module.
!
!
!===============================================================================
!
  subroutine finalize_fitsio()

    implicit none
!
!-------------------------------------------------------------------------------
!

!-------------------------------------------------------------------------------
!
  end subroutine finalize_fitsio
!
!===============================================================================
!
! subroutine CHECK_IF_EXISTS:
! --------------------------
!
!   Subroutine checks if the file exists and returns its full path.
!
!
!===============================================================================
!
  subroutine check_if_exists(ddir, cname, fname, info)

    implicit none

    character(len=*)  , intent(in)  :: ddir, cname
    character(len=128), intent(out) :: fname
    logical           , intent(out) :: info
!
!-------------------------------------------------------------------------------
!
    write(fname,"(a,a,'.fits')") trim(ddir), trim(cname)
    inquire(file = fname, exist = info)
    if (.not. info) then
      write(fname,"(a,a,'.fits.gz')") trim(ddir), trim(cname)
      inquire(file = fname, exist = info)
      if (.not. info) then
        write(*,"('Data file for ',a,' does not exists under ',a,'!')")        &
                trim(cname), trim(ddir)
        stop
      end if
    end if

!-------------------------------------------------------------------------------
!
  end subroutine check_if_exists
!
!===============================================================================
!
! subroutine GET_DIMENSIONS:
! -------------------------
!
!   Subroutine reads the dimensions of data files.
!
!
!===============================================================================
!
  subroutine get_dimensions(dm)

    implicit none

! subroutine arguments
!
    integer, dimension(3), intent(out) :: dm

! FITS variables
!
    integer :: status, iunit, bsize, bpix, naxes
!
!-------------------------------------------------------------------------------
!
    status = 0
    call ftgiou(iunit, status)
    call ftopen(iunit, uxfname, 0, bsize, status)
    call ftgipr(iunit, 3, bpix, naxes, dm, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)

!-------------------------------------------------------------------------------
!
  end subroutine get_dimensions
!
!===============================================================================
!
! subroutine READ_FIELD:
! ---------------------
!
!   Subroutine reads the data of selected field.
!
!
!===============================================================================
!
  subroutine read_field(var, qty)

    implicit none

! subroutine arguments
!
    character(len=*)                  , intent(in)    :: var
    real(kind=PREC) , dimension(:,:,:), intent(inout) :: qty

! local variables
!
    logical               :: info
    character(len=128)    :: fname
    integer               :: status, iunit, bsize, bpix, naxes, nn
    integer, dimension(3) :: dm

! allocatable arrays
!
    real(kind=4), dimension(:,:,:), allocatable :: arr4
    real(kind=8), dimension(:,:,:), allocatable :: arr8
!
!-------------------------------------------------------------------------------
!
    select case(trim(var))
    case("ux", "velx", "x-velocity")
      fname = uxfname
    case("uy", "vely", "y-velocity")
      fname = uyfname
    case("uz", "velz", "z-velocity")
      fname = uzfname
    case("bx", "magx", "x-magnetic")
      fname = bxfname
    case("by", "magy", "y-magnetic")
      fname = byfname
    case("bz", "magz", "z-magnetic")
      fname = bzfname
    end select

    dm = 1
    status = 0
    call ftgiou(iunit, status)
    call ftopen(iunit, fname, 0, bsize, status)
    call ftgipr(iunit, 3, bpix, naxes, dm, status)
    nn = product(dm)
    select case(bpix)
    case(-32)
      allocate(arr4(dm(1),dm(2),dm(3)))
      call ftgpve(iunit, 0, 1, nn, 0.0, arr4, info, status)
      qty = real(arr4, kind=PREC)
      deallocate(arr4)
    case(-64)
      allocate(arr8(dm(1),dm(2),dm(3)))
      call ftgpve(iunit, 0, 1, nn, 0.0, arr8, info, status)
      qty = real(arr8, kind=PREC)
      deallocate(arr8)
    end select
    call ftclos(iunit, status)
    call ftfiou(iunit, status)

!-------------------------------------------------------------------------------
!
  end subroutine read_field
!
!===============================================================================
!
! subroutine WRITE_DATA_INTEGER_1D:
! --------------------------------
!
!   Subroutine write the input data of integer type to file.
!
!
!===============================================================================
!
  subroutine write_data_integer_1d(fname, qty, tunit, temperature)

    implicit none

! subroutine arguments
!
    character(len=*)             , intent(in) :: fname
    integer(kind=8), dimension(:), intent(in) :: qty
    real(kind=PREC), optional    , intent(in) :: tunit, temperature

! local variables
!
    logical               :: info
    integer               :: status, iunit, bsize, bpix, naxes, nn
    integer, dimension(1) :: dm
!
!-------------------------------------------------------------------------------
!
    status = 0
    call ftgiou(iunit, status)
    call ftinit(iunit, fname, 1, status)
    dm(1) = size(qty, 1)
    nn    = size(qty)
    bpix  = 64
    call ftphpr(iunit, .true., bpix, 1, dm, 0, 1, .true., status)
    call ftpprk(iunit, 1, 1, nn, qty, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)

!-------------------------------------------------------------------------------
!
  end subroutine write_data_integer_1d
!
!===============================================================================
!
! subroutine WRITE_DATA_INTEGER_2D:
! --------------------------------
!
!   Subroutine write the input data of integer type to file.
!
!
!===============================================================================
!
  subroutine write_data_integer_2d(fname, qty, tunit, temperature)

    implicit none

! subroutine arguments
!
    character(len=*)               , intent(in) :: fname
    integer(kind=8), dimension(:,:), intent(in) :: qty
    real(kind=PREC), optional      , intent(in) :: tunit, temperature

! local variables
!
    logical               :: info
    integer               :: status, iunit, bsize, bpix, naxes, nn
    integer, dimension(2) :: dm
!
!-------------------------------------------------------------------------------
!
    status = 0
    call ftgiou(iunit, status)
    call ftinit(iunit, fname, 1, status)
    dm(1) = size(qty, 1)
    dm(2) = size(qty, 2)
    nn    = size(qty)
    bpix  = 64
    call ftphpr(iunit, .true., bpix, 2, dm, 0, 1, .true., status)
    call ftpprk(iunit, 1, 1, nn, qty, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)

!-------------------------------------------------------------------------------
!
  end subroutine write_data_integer_2d
!
!===============================================================================
!
! subroutine WRITE_DATA_REAL_1D:
! -----------------------------
!
!   Subroutine write the input data to file.
!
!
!===============================================================================
!
  subroutine write_data_real_1d(fname, qty, tunit, temperature)

    implicit none

! subroutine arguments
!
    character(len=*)              , intent(in) :: fname
    real(kind=PREC) , dimension(:), intent(in) :: qty
    real(kind=PREC), optional     , intent(in) :: tunit, temperature

! local variables
!
    logical               :: info
    integer               :: status, iunit, bsize, bpix, naxes, nn
    integer, dimension(1) :: dm
!
!-------------------------------------------------------------------------------
!
    status = 0
    call ftgiou(iunit, status)
    call ftinit(iunit, fname, 1, status)
    dm(1) = size(qty, 1)
    nn    = size(qty)
    bpix = - 8 * PREC
    call ftphpr(iunit, .true., bpix, 1, dm, 0, 1, .true., status)
#if PREC == 4
    if (present(tunit)) then
      call ftpkye(iunit, 'TUNIT', tunit, -8, 'Time unit in seconds', status)
    end if
    if (present(temperature)) then
      call ftpkye(iunit, 'TEMP', temperature, -8, 'Initial thermal temperature', status)
    end if
    call ftppre(iunit, 1, 1, nn, qty, status)
#endif
#if PREC == 8
    call ftpprd(iunit, 1, 1, nn, qty, status)
    if (present(tunit)) then
      call ftpkyd(iunit, 'TUNIT', tunit, -8, 'Time unit in seconds', status)
    end if
    if (present(temperature)) then
      call ftpkyd(iunit, 'TEMP', temperature, -8, 'Initial thermal temperature', status)
    end if
    call ftpprd(iunit, 1, 1, nn, qty, status)
#endif
    call ftclos(iunit, status)
    call ftfiou(iunit, status)

!-------------------------------------------------------------------------------
!
  end subroutine write_data_real_1d
!
!===============================================================================
!
! subroutine WRITE_DATA_REAL_2D:
! -----------------------------
!
!   Subroutine write the input data to file.
!
!
!===============================================================================
!
  subroutine write_data_real_2d(fname, qty, tunit, temperature)

    implicit none

! subroutine arguments
!
    character(len=*)                , intent(in) :: fname
    real(kind=PREC) , dimension(:,:), intent(in) :: qty
    real(kind=PREC), optional       , intent(in) :: tunit, temperature

! local variables
!
    logical               :: info
    integer               :: status, iunit, bsize, bpix, naxes, nn
    integer, dimension(2) :: dm
!
!-------------------------------------------------------------------------------
!
    status = 0
    call ftgiou(iunit, status)
    call ftinit(iunit, fname, 1, status)
    dm(1) = size(qty, 1)
    dm(2) = size(qty, 2)
    nn    = size(qty)
    bpix = - 8 * PREC
    call ftphpr(iunit, .true., bpix, 2, dm, 0, 1, .true., status)
#if PREC == 4
    if (present(tunit)) then
      call ftpkye(iunit, 'TUNIT', tunit, -8, 'Time unit in seconds', status)
    end if
    if (present(temperature)) then
      call ftpkye(iunit, 'TEMP', temperature, -8, 'Initial thermal temperature', status)
    end if
    call ftppre(iunit, 1, 1, nn, qty, status)
#endif
#if PREC == 8
    if (present(tunit)) then
      call ftpkyd(iunit, 'TUNIT', tunit, -8, 'Time unit in seconds', status)
    end if
    if (present(temperature)) then
      call ftpkyd(iunit, 'TEMP', temperature, -8, 'Initial thermal temperature', status)
    end if
    call ftpprd(iunit, 1, 1, nn, qty, status)
#endif
    call ftclos(iunit, status)
    call ftfiou(iunit, status)

!-------------------------------------------------------------------------------
!
  end subroutine write_data_real_2d
!
!===============================================================================
!
! subroutine WRITE_DATA_REAL_3D:
! -----------------------------
!
!   Subroutine write the input data to file.
!
!
!===============================================================================
!
  subroutine write_data_real_3d(fname, qty, tunit, temperature)

    implicit none

! subroutine arguments
!
    character(len=*)                  , intent(in) :: fname
    real(kind=PREC) , dimension(:,:,:), intent(in) :: qty
    real(kind=PREC), optional         , intent(in) :: tunit, temperature

! local variables
!
    logical               :: info
    integer               :: status, iunit, bsize, bpix, naxes, nn
    integer, dimension(3) :: dm
!
!-------------------------------------------------------------------------------
!
    status = 0
    call ftgiou(iunit, status)
    call ftinit(iunit, fname, 1, status)
    dm(1) = size(qty, 1)
    dm(2) = size(qty, 2)
    dm(3) = size(qty, 3)
    nn    = size(qty)
    bpix = - 8 * PREC
    call ftphpr(iunit, .true., bpix, 3, dm, 0, 1, .true., status)
#if PREC == 4
    if (present(tunit)) then
      call ftpkye(iunit, 'TUNIT', tunit, -8, 'Time unit in seconds', status)
    end if
    if (present(temperature)) then
      call ftpkye(iunit, 'TEMP', temperature, -8, 'Initial thermal temperature', status)
    end if
    call ftppre(iunit, 1, 1, nn, qty, status)
#endif
#if PREC == 8
    if (present(tunit)) then
      call ftpkyd(iunit, 'TUNIT', tunit, -8, 'Time unit in seconds', status)
    end if
    if (present(temperature)) then
      call ftpkyd(iunit, 'TEMP', temperature, -8, 'Initial thermal temperature', status)
    end if
    call ftpprd(iunit, 1, 1, nn, qty, status)
#endif
    call ftclos(iunit, status)
    call ftfiou(iunit, status)

!-------------------------------------------------------------------------------
!
  end subroutine write_data_real_3d

!===============================================================================
!
end module
