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
!! module: INTERPOLATIONS
!!
!!  This module handles interpolation methods.
!!
!!******************************************************************************
!
module interpolations

  implicit none

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INTERPOLATE_NEAR:
! ---------------------------
!
!   Subroutine interpolates fields at a given position using near-neighbor.
!
!   Arguments:
!
!     dm     - dimensions of the interpolated fields uu and bb (input);
!     uu, bb - the arrays of velocity and magnetic fields (input);
!     pos    - the interpolation position step (input);
!     u , b  - the interpolated fields (output);
!
!===============================================================================
!
  subroutine interpolate_near(dm, uu, bb, pos, u, b)

    implicit none
    !$acc routine (interpolate_near) seq

! subroutine arguments
!
    integer        , dimension(3)                  , intent(in)  :: dm
    real(kind=PREC), dimension(dm(1),dm(2),dm(3),3), intent(in)  :: uu, bb
    real(kind=PREC), dimension(3)                  , intent(in)  :: pos
    real(kind=PREC), dimension(3)                  , intent(out) :: u, b

! local variables
!
    integer, dimension(3) :: p
!
!-------------------------------------------------------------------------------
!
    p      = floor(dm(1:3) * (pos(1:3) - floor(pos(1:3)))) + 1

    u(1:3) = uu(p(1),p(2),p(3),1:3)
    b(1:3) = bb(p(1),p(2),p(3),1:3)

!-------------------------------------------------------------------------------
!
  end subroutine
!
!===============================================================================
!
! subroutine INTERPOLATE_LIN:
! --------------------------
!
!   Subroutine interpolates fields at a given position using trilinear
!   interpolation.
!
!   Arguments:
!
!     dm     - dimensions of the interpolated fields uu and bb (input);
!     uu, bb - the arrays of velocity and magnetic fields (input);
!     pos    - the interpolation position step (input);
!     u , b  - the interpolated fields (output);
!
!===============================================================================
!
  subroutine interpolate_lin(dm, uu, bb, pos, u, b)

    implicit none
    !$acc routine (interpolate_lin) seq

! subroutine arguments
!
    integer        , dimension(3)                  , intent(in)  :: dm
    real(kind=PREC), dimension(dm(1),dm(2),dm(3),3), intent(in)  :: uu, bb
    real(kind=PREC), dimension(3)                  , intent(in)  :: pos
    real(kind=PREC), dimension(3)                  , intent(out) :: u, b

! local variables
!
!     integer                           :: m
    integer        , dimension(3)     :: pl, pr, l, r
    real(kind=PREC), dimension(3)     :: xp, xd, xu
!     real(kind=PREC), dimension(2,2,6) :: t2
!     real(kind=PREC), dimension(2,6)   :: t1

    real(kind=PREC), dimension(8)     :: c
!     real(kind=PREC), dimension(8)     :: q1, c1, p
!     real(kind=PREC), dimension(8,8), parameter   :: b1 = reshape(                        &
!        (/  1.0d+00, -1.0d+00, -1.0d+00, -1.0d+00,  1.0d+00,  1.0d+00,  1.0d+00, -1.0d+00 &
!         ,  0.0d+00,  0.0d+00,  0.0d+00,  1.0d+00,  0.0d+00, -1.0d+00, -1.0d+00,  1.0d+00 &
!         ,  0.0d+00,  0.0d+00,  1.0d+00,  0.0d+00, -1.0d+00, -1.0d+00,  0.0d+00,  1.0d+00 &
!         ,  0.0d+00,  0.0d+00,  0.0d+00,  0.0d+00,  0.0d+00,  1.0d+00,  0.0d+00, -1.0d+00 &
!         ,  0.0d+00,  1.0d+00,  0.0d+00,  0.0d+00, -1.0d+00,  0.0d+00, -1.0d+00,  1.0d+00 &
!         ,  0.0d+00,  0.0d+00,  0.0d+00,  0.0d+00,  0.0d+00,  0.0d+00,  1.0d+00, -1.0d+00 &
!         ,  0.0d+00,  0.0d+00,  0.0d+00,  0.0d+00,  1.0d+00,  0.0d+00,  0.0d+00, -1.0d+00 &
!         ,  0.0d+00,  0.0d+00,  0.0d+00,  0.0d+00,  0.0d+00,  0.0d+00,  0.0d+00,  1.0d+00 /), (/ 8,8 /))
!
!-------------------------------------------------------------------------------
!
    xd = pos(1:3) + 1.0d+00 / dm(1:3)
    xp = dm(1:3) * ( xd(1:3) - floor( xd(1:3)))
    r  = floor(xp) + 1

    xp = dm(1:3) * (pos(1:3) - floor(pos(1:3)))
    pl = floor(xp)
    xd = xp - pl
    xu = 1.0d+00 - xd
    l  = pl + 1

!     write (*,"(a,3es15.5)") '0: ', xd
!     write (*,"(a,3es20.8)") '0: ', uu(l(1),l(2),l(3),1:3)

!     t2(1,1,1:3) = xu(3) * uu(l(1),l(2),l(3),1:3) + xd(3) * uu(l(1),l(2),r(3),1:3)
!     t2(1,2,1:3) = xu(3) * uu(l(1),r(2),l(3),1:3) + xd(3) * uu(l(1),r(2),r(3),1:3)
!     t2(2,1,1:3) = xu(3) * uu(r(1),l(2),l(3),1:3) + xd(3) * uu(r(1),l(2),r(3),1:3)
!     t2(2,2,1:3) = xu(3) * uu(r(1),r(2),l(3),1:3) + xd(3) * uu(r(1),r(2),r(3),1:3)
!
!     t2(1,1,4:6) = xu(3) * bb(l(1),l(2),l(3),1:3) + xd(3) * bb(l(1),l(2),r(3),1:3)
!     t2(1,2,4:6) = xu(3) * bb(l(1),r(2),l(3),1:3) + xd(3) * bb(l(1),r(2),r(3),1:3)
!     t2(2,1,4:6) = xu(3) * bb(r(1),l(2),l(3),1:3) + xd(3) * bb(r(1),l(2),r(3),1:3)
!     t2(2,2,4:6) = xu(3) * bb(r(1),r(2),l(3),1:3) + xd(3) * bb(r(1),r(2),r(3),1:3)
!
!     t1(1,1:6)   = xu(2) * t2(1,1,1:6)            + xd(2) * t2(1,2,1:6)
!     t1(2,1:6)   = xu(2) * t2(2,1,1:6)            + xd(2) * t2(2,2,1:6)
!
!     u(1:3)      = xu(1) * t1(1,1:3)              + xd(1) * t1(2,1:3)
!     b(1:3)      = xu(1) * t1(1,4:6)              + xd(1) * t1(2,4:6)

!     write (*,"(a,3es20.8)") '1: ', u(1:3)
!
!     q1 = (/ 1.0d+00, xd(1), xd(2), xd(3), xd(1) * xd(2), xd(2) * xd(3)         &
!                    , xd(3) * xd(1), xd(1) * xd(2) * xd(3) /)
!     do m = 1, 3
!       p(1:8) = (/ uu(l(1),l(2),l(3),m), uu(l(1),l(2),r(3),m)                   &
!                 , uu(l(1),r(2),l(3),m), uu(l(1),r(2),r(3),m)                   &
!                 , uu(r(1),l(2),l(3),m), uu(r(1),l(2),r(3),m)                   &
!                 , uu(r(1),r(2),l(3),m), uu(r(1),r(2),r(3),m) /)
!       u(m)   = dot_product(q1(:), matmul(b1(:,:), p(:)))
!       p(1:8) = (/ bb(l(1),l(2),l(3),m), bb(l(1),l(2),r(3),m)                   &
!                 , bb(l(1),r(2),l(3),m), bb(l(1),r(2),r(3),m)                   &
!                 , bb(r(1),l(2),l(3),m), bb(r(1),l(2),r(3),m)                   &
!                 , bb(r(1),r(2),l(3),m), bb(r(1),r(2),r(3),m) /)
!       b(m)   = dot_product(q1(:), matmul(b1(:,:), p(:)))
!     end do
!
!     write (*,"(a,3es20.8)") '2: ', u(1:3)
!
!     q1 = (/ 1.0d+00, xd(1), xd(2), xd(3), xd(1) * xd(2), xd(2) * xd(3)         &
!                    , xd(3) * xd(1), xd(1) * xd(2) * xd(3) /)
!     c1 = matmul(q1, b1)
!     do m = 1, 3
!       p(1:8) = (/ uu(l(1),l(2),l(3),m), uu(l(1),l(2),r(3),m)                   &
!                 , uu(l(1),r(2),l(3),m), uu(l(1),r(2),r(3),m)                   &
!                 , uu(r(1),l(2),l(3),m), uu(r(1),l(2),r(3),m)                   &
!                 , uu(r(1),r(2),l(3),m), uu(r(1),r(2),r(3),m) /)
!       u(m)   = dot_product(c1(:), p(:))
!       p(1:8) = (/ bb(l(1),l(2),l(3),m), bb(l(1),l(2),r(3),m)                   &
!                 , bb(l(1),r(2),l(3),m), bb(l(1),r(2),r(3),m)                   &
!                 , bb(r(1),l(2),l(3),m), bb(r(1),l(2),r(3),m)                   &
!                 , bb(r(1),r(2),l(3),m), bb(r(1),r(2),r(3),m) /)
!       b(m)   = dot_product(c1(:), p(:))
!     end do
!
!     write (*,"(a,3es20.8)") '3: ', u(1:3)
!
    c(1:8) = (/ xu(1) * xu(2) * xu(3), xu(1) * xu(2) * xd(3)                   &
              , xu(1) * xd(2) * xu(3), xu(1) * xd(2) * xd(3)                   &
              , xd(1) * xu(2) * xu(3), xd(1) * xu(2) * xd(3)                   &
              , xd(1) * xd(2) * xu(3), xd(1) * xd(2) * xd(3) /)

    u(1:3) = c(1) * uu(l(1),l(2),l(3),1:3) + c(2) * uu(l(1),l(2),r(3),1:3)     &
           + c(3) * uu(l(1),r(2),l(3),1:3) + c(4) * uu(l(1),r(2),r(3),1:3)     &
           + c(5) * uu(r(1),l(2),l(3),1:3) + c(6) * uu(r(1),l(2),r(3),1:3)     &
           + c(7) * uu(r(1),r(2),l(3),1:3) + c(8) * uu(r(1),r(2),r(3),1:3)

    b(1:3) = c(1) * bb(l(1),l(2),l(3),1:3) + c(2) * bb(l(1),l(2),r(3),1:3)     &
           + c(3) * bb(l(1),r(2),l(3),1:3) + c(4) * bb(l(1),r(2),r(3),1:3)     &
           + c(5) * bb(r(1),l(2),l(3),1:3) + c(6) * bb(r(1),l(2),r(3),1:3)     &
           + c(7) * bb(r(1),r(2),l(3),1:3) + c(8) * bb(r(1),r(2),r(3),1:3)
!
!     write (*,"(a,3es20.8)") '4: ', u(1:3)

!-------------------------------------------------------------------------------
!
  end subroutine

!===============================================================================
!
end module
