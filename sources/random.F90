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
!! module: RANDOM
!!
!!  This module handles random number generators.
!!
!!******************************************************************************
!
module random

! declare all module variables as implicit
!
  implicit none

! by default everything is private
!
  private

! declare public subroutines
!
  public :: initialize_random, finalize_random
  public :: uniform, normal, moments

  contains
!
!===============================================================================
!
! subroutine INITIALIZE_RANDOM:
! ----------------------------
!
!   subroutine initializes random number generator;
!
!===============================================================================
!
  subroutine initialize_random()

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
!
  end subroutine initialize_random
!
!===============================================================================
!
! subroutine FINALIZE_RANDOM:
! --------------------------
!
!   subroutine releases memory allocated by random number generator variables;
!
!===============================================================================
!
  subroutine finalize_random()

! declare all variables as implicit
!
    implicit none
!
!-------------------------------------------------------------------------------
!

!-------------------------------------------------------------------------------
!
  end subroutine finalize_random
!
!===============================================================================
!
! subroutine UNIFORM:
! ------------------
!
!   Subroutine fills an input 2D array with uniform random numbers between
!   0.0 and 1.0.
!
!   Arguments:
!
!     arr - the array of numbers (input/output);
!
!===============================================================================
!
  subroutine uniform(amp, arr)

    implicit none

! subroutine arguments
!
    real(kind=PREC)                , intent(in)    :: amp
    real(kind=PREC), dimension(:,:), intent(inout) :: arr
!
!-------------------------------------------------------------------------------
!
    call random_number(arr)
    arr = amp * arr

!-------------------------------------------------------------------------------
!
  end subroutine
!
!===============================================================================
!
! subroutine NORMAL:
! -----------------
!
!   Subroutine fills an input 2D array with normal distribution random numbers
!   with a given mean and sigma.
!
!   Arguments:
!
!     mu    - distribution mean (input);
!     sigma - distribution standard deviation (input);
!     arr   - the array of numbers (input/output);
!
!===============================================================================
!
  subroutine normal(mu, sigma, arr)

    implicit none

! subroutine arguments
!
    real(kind=PREC)                , intent(in)    :: mu, sigma
    real(kind=PREC), dimension(:,:), intent(inout) :: arr

! local variables
!
    integer :: i, j, m, n
    real(kind=PREC), dimension(:,:), allocatable :: u, v, s
!
!-------------------------------------------------------------------------------
!
    m = size(arr,1)
    n = size(arr,2) / 2
    allocate(u(m,n), v(m,n), s(m,n))
    call random_number(u)
    u(:,:) = 2.0d+00 * u(:,:) - 1.0d+00
    call random_number(v)
    v(:,:) = 2.0d+00 * v(:,:) - 1.0d+00
    s(:,:) = u(:,:)**2 + v(:,:)**2
    do j = 1, n
      do i = 1, m
      do while (s(i,j) > 1.0d+00 .or. s(i,j) == 0.0d+00)
          call random_number(u(i,j))
          call random_number(v(i,j))
          u(i,j) = 2.0d+00 * u(i,j) - 1.0d+00
          v(i,j) = 2.0d+00 * v(i,j) - 1.0d+00
          s(i,j) = u(i,j)**2 + v(i,j)**2
        end do
      end do
    end do
    arr(:,   :n) = u(:,:) * sqrt(-2.0d+00 * log(s(:,:)) / s(:,:)) * sigma + mu
    arr(:,n+1: ) = v(:,:) * sqrt(-2.0d+00 * log(s(:,:)) / s(:,:)) * sigma + mu
    deallocate(u, v, s)

!-------------------------------------------------------------------------------
!
  end subroutine
!
!===============================================================================
!
! subroutine MOMENTS:
! ------------------
!
!   Subroutine calculates statistical moments: mean, standard deviation,
!   skewness and kurtosis of the input array.
!
!   Arguments:
!
!     arr - the array of numbers (input);
!     mom - the vector of statistical moments (output);
!
!===============================================================================
!
  subroutine moments(arr, mom)

    implicit none

! subroutine arguments
!
    real(kind=PREC), dimension(:), intent(in)  :: arr
    real(kind=PREC), dimension(4), intent(out) :: mom

! allocatable arrays
!
    real(kind=PREC), dimension(:), allocatable :: tmp
!
!-------------------------------------------------------------------------------
!
    mom(1) = mean(arr)
    allocate(tmp(size(arr)))
    tmp(:) = arr(:) - mom(1)
    mom(2) = sqrt(mean(tmp(:)**2))
    if (mom(2) /= 0.0d+00) tmp(:) = tmp(:) / mom(2)
    mom(3) = mean(tmp(:)**3)
    mom(4) = mean(tmp(:)**4) - 3.0d+00
    deallocate(tmp)

!-------------------------------------------------------------------------------
!
  end subroutine
!
!===============================================================================
!
! function MEAN:
! -------------
!
!   Function calculates the mean of elements of the input array.
!
!   Arguments:
!
!     arr - the array of numbers (input);
!
!===============================================================================
!
  function mean(arr) result(res)

    implicit none

! function arguments
!
    real(kind=PREC), dimension(:), intent(in)  :: arr

! local variables
!
    real(kind=PREC) :: res
!
!-------------------------------------------------------------------------------
!
    res = neumaiersum(arr(:)) / size(arr)

    return

!-------------------------------------------------------------------------------
!
  end function
!
!===============================================================================
!
! function KAHANSUM:
! -----------------
!
!   Function calculates sum of the elements of input array using Kahan summation
!   algorithm in order to reduce the numerical error.
!
!   Arguments:
!
!     arr - the array of numbers (input);
!
!===============================================================================
!
  function kahansum(arr) result(s)

    implicit none

! function arguments
!
    real(kind=PREC), dimension(:), intent(in)  :: arr

! local variables
!
    integer         :: m, i
    real(kind=PREC) :: s, y, t, c
!
!-------------------------------------------------------------------------------
!
    s = 0.0d+00
    c = 0.0d+00
    m = size(arr)
    do i = 1, m
      y = arr(i) - c
      t = s + y
      c = (t - s) - y
      s = t
    end do

    return

!-------------------------------------------------------------------------------
!
  end function
!
!===============================================================================
!
! function NEUMAIERSUM:
! --------------------
!
!   Function calculates sum of the elements of input array using Neumaier
!   summation algorithm in order to reduce the numerical error.
!
!   Arguments:
!
!     arr - the array of numbers (input);
!
!===============================================================================
!
  function neumaiersum(arr) result(s)

    implicit none

! function arguments
!
    real(kind=PREC), dimension(:), intent(in)  :: arr

! local variables
!
    integer         :: m, i
    real(kind=PREC) :: s, t, c
!
!-------------------------------------------------------------------------------
!
    m = size(arr)

    s = arr(1)
    c = 0.0d+00
    do i = 2, m
      t = s + arr(i)
      if (abs(s) >= abs(arr(i))) then
        c = c + (s - t) + arr(i)
      else
        c = c + (arr(i) - t) + s
      end if
      s = t
    end do

    s = s + c

    return

!-------------------------------------------------------------------------------
!
  end function

!===============================================================================
!
end module random
