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
!! module: ACCELERATIONS
!!
!!  This module handles acceleration terms.
!!
!!******************************************************************************
!
module accelerations

  implicit none

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine ACCELERATION:
! -----------------------
!
!   Subroutine calculates the Lorentz force from a uniform magnetic field line.
!
!   Arguments:
!
!     vun    - the unit of velocity in terms of the speed of light, which
!              defines how to rescale the field velocity and how to convert
!              particle speed back to plasma velocity units (input);
!     qom    - the charge over mass coefficient (input);
!     dm     - dimensions of the field components uu and bb (input);
!     uu, bb - the arrays of velocity and magnetic fields (input);
!     vec    - particle state vector, ie. position and velocity (input);
!     acc    - the velocity and acceleration vector (output);
!
!===============================================================================
!
  subroutine acceleration(qom, vun, dm, uu, bb, vec, acc)

    use interpolations, only : interpolate_near, interpolate_lin

    implicit none
    !$acc routine (acceleration) seq

! subroutine arguments
!
    real(kind=PREC)                                , intent(in)  :: qom, vun
    integer        , dimension(3)                  , intent(in)  :: dm
    real(kind=PREC), dimension(dm(1),dm(2),dm(3),3), intent(in)  :: uu, bb
    real(kind=PREC), dimension(6)                  , intent(in)  :: vec
    real(kind=PREC), dimension(6)                  , intent(out) :: acc

! local variables
!
    real(kind=PREC), dimension(3)  :: u, b, v, w
    real(kind=PREC)                :: g
!
!-------------------------------------------------------------------------------
!
    call interpolate_lin(dm(:), uu(:,:,:,:), bb(:,:,:,:), vec(1:3), u(:), b(:))

    g        = sqrt(dot_product(vec(4:6), vec(4:6)) + 1.0d+00)
    v(1:3)   = vec(4:6) / g

    w(1:3)   = v(1:3) - u(1:3) * vun

    acc(1:3) = v(1:3) / vun
    acc(4)   = qom * (w(2) * b(3) - w(3) * b(2))
    acc(5)   = qom * (w(3) * b(1) - w(1) * b(3))
    acc(6)   = qom * (w(1) * b(2) - w(2) * b(1))

!-------------------------------------------------------------------------------
!
  end subroutine acceleration

!===============================================================================
!
end module
