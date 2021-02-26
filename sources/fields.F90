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
!! module: FIELDS
!!
!!  This module handles magnetohydrodynamic field components.
!!
!!******************************************************************************
!
module fields

  implicit none

  integer, dimension(3) :: dm

  real(kind=PREC), dimension(:,:,:,:), allocatable :: uu, bb

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INITIALIZE_FIELDS:
! ----------------------------
!
!   Subroutine allocates space for fields.
!
!
!===============================================================================
!
  subroutine initialize_fields(uunit, bunit)

    use fitsio, only : get_dimensions, read_field

    implicit none

    real(kind=PREC), intent(in) :: uunit, bunit
!
!-------------------------------------------------------------------------------
!
    call get_dimensions(dm)

    allocate(uu(dm(1), dm(2), dm(3), 3))
    allocate(bb(dm(1), dm(2), dm(3), 3))

    call read_field('ux', uu(:,:,:,1))
    call read_field('uy', uu(:,:,:,2))
    call read_field('uz', uu(:,:,:,3))
    call read_field('bx', bb(:,:,:,1))
    call read_field('by', bb(:,:,:,2))
    call read_field('bz', bb(:,:,:,3))

    uu = uunit * uu
    bb = bunit * bb

!-------------------------------------------------------------------------------
!
  end subroutine initialize_fields
!
!===============================================================================
!
! subroutine FINALIZE_FIELDS:
! --------------------------
!
!   Subroutine deallocates space for fields.
!
!
!===============================================================================
!
  subroutine finalize_fields()

    implicit none
!
!-------------------------------------------------------------------------------
!
    deallocate(uu)
    deallocate(bb)

!-------------------------------------------------------------------------------
!
  end subroutine finalize_fields

!===============================================================================
!
end module
