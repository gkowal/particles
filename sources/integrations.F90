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
!! module: INTEGRATIONS
!!
!!  This module handles integration methods.
!!
!!******************************************************************************
!
module integrations

  implicit none

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine INTEGRATE_RK4:
! ------------------------
!
!   Subroutine integrate particle trajectory using RK4 method.
!
!   Arguments:
!
!     m      - the number of steps to store (input);
!     dm     - domensions of the field components (input);
!     params - the vector of parameters (input);
!     time   - the vector of snapshot times (input);
!     uu, bb - the components of velocity and magnetic fields (input);
!     state  - the particle state vector, position, velocity, and others
!              (input/output);
!
!===============================================================================
!
  subroutine integrate_rk4(m, dm, params, time, uu, bb, state, cnt)

    use accelerations, only : acceleration

    implicit none
    !$acc routine (integrate_rk4) seq

! subroutine arguments
!
    integer                                        , intent(in)    :: m
    integer, dimension(3)                          , intent(in)    :: dm
    real(kind=PREC), dimension(8)                  , intent(in)    :: params
    real(kind=PREC), dimension(m)                  , intent(in)    :: time
    real(kind=PREC), dimension(dm(1),dm(2),dm(3),3), intent(in)    :: uu, bb
    real(kind=PREC), dimension(8,m)                , intent(inout) :: state
    integer(kind=8)                                , intent(out)   :: cnt

! local variables
!
    integer                       :: i
    real(kind=PREC)               :: qom, vun, tmx, dti, dtm, tol
    real(kind=PREC)               :: t, dt, del, htol, fl, fr
    real(kind=PREC), dimension(6) :: s1, s2, s3, s4, s5, ds
    real(kind=PREC), dimension(6) :: k1, k2, k3, k4, k5
    real(kind=PREC), dimension(6) :: ac
!
!-------------------------------------------------------------------------------
!
    qom = params(1)
    vun = params(2)
    tmx = params(3)
    dti = params(4)
    dtm = params(5)
    tol = params(6)

    cnt        = 0
    htol       = 0.5d+00 * tol
    t          = 0.0d+00
    dt         = min(dtm, dti)
    i          = 1
    s1(1:6)    = state(1:6,i)
    state(7,i) = sqrt(dot_product(s1(4:6), s1(4:6)) + 1.0d+00)
    state(8,i) = dt
    i          = i + 1
    call acceleration(qom, vun, dm, uu, bb, s1, ac)

    do while(t < tmx)

      k1(:) = dt * ac(:)

      s2(:) = s1(:) + 0.5d+00 * k1(:)

      call acceleration(qom, vun, dm, uu, bb, s2, ac)

      k2(:) = dt * ac(:)

      s3(:) = s1(:) + 0.5d+00 * k2(:)

      call acceleration(qom, vun, dm, uu, bb, s3, ac)

      k3(:) = dt * ac(:)

      s4(:) = s1(:) + k3(:)

      call acceleration(qom, vun, dm, uu, bb, s4, ac)

      k4(:) = dt * ac(:)

      s5(:) = s1(:) + (2.0d+00 * (k2(:) + k3(:)) + k1(:) + k4(:)) / 6.0d+00

      call acceleration(qom, vun, dm, uu, bb, s5, ac)

      k4(:) = k4(:) - dt * ac(:)

      del   = sqrt(dot_product(k4, k4)) / 6.0d+00

      if (del <= tol) then
        cnt   = cnt + 1
        t     = t + dt
        do while (t >= time(i) .and. i <= m)
          fr = (t - time(i)) / dt
          fl = 1.0d+00 - fr
          ds = s5(1:6) - s1(1:6)
          state(1:6,i) = fr * s1(1:6) + fl * s5(1:6)                           &
                       + fl * fr * ((k1(1:6) - ds(1:6)) * fr                   &
                                  + (ds(1:6) - k4(1:6)) * fl)
          state(  7,i) = sqrt(dot_product(state(4:6,i), state(4:6,i)) + 1.0d+00)
          state(  8,i) = dt
          i  = i + 1
        end do
        s1(:) = s5(:)
      else
        ac(:) = k1(:) / dt
      end if

      dt = min(2.0d+00 * dt, dtm, dt * sqrt(sqrt(htol / max(1.0d-16, del))))

    end do

!-------------------------------------------------------------------------------
!
  end subroutine
!
!===============================================================================
!
! subroutine INTEGRATE_RK5:
! ------------------------
!
!   Subroutine integrate particle trajectory using 5th order RK-Cash-Karp
!   method with adaptive time step.
!
!   Arguments:
!
!     m      - the number of steps to store (input);
!     dm     - domensions of the field components (input);
!     params - the vector of parameters (input);
!     time   - the vector of snapshot times (input);
!     uu, bb - the components of velocity and magnetic fields (input);
!     state  - the particle state vector, position, velocity, and others
!              (input/output);
!
!===============================================================================
!
  subroutine integrate_rk5(m, dm, params, time, uu, bb, state, cnt)

    use accelerations, only : acceleration

    implicit none
    !$acc routine (integrate_rk5) seq

! subroutine arguments
!
    integer                                        , intent(in)    :: m
    integer, dimension(3)                          , intent(in)    :: dm
    real(kind=PREC), dimension(8)                  , intent(in)    :: params
    real(kind=PREC), dimension(m)                  , intent(in)    :: time
    real(kind=PREC), dimension(dm(1),dm(2),dm(3),3), intent(in)    :: uu, bb
    real(kind=PREC), dimension(8,m)                , intent(inout) :: state
    integer(kind=8)                                , intent(out)   :: cnt

! local variables
!
    integer                         :: i, j, n
    real(kind=PREC)                 :: qom, vun, tmx, dti, dtm, tol
    real(kind=PREC)                 :: t, dt, dtn, del, fl, fr
    real(kind=PREC), dimension(6)   :: si, sf, er, sr, ds
    real(kind=PREC), dimension(6)   :: ac
    real(kind=PREC), dimension(6,6) :: s, k

! parameters
!
    real(kind=PREC), dimension(5,5), parameter :: b = reshape(                 &
            (/ 1.10592d+05, 4.14720d+04, 1.65888d+05,-1.12640d+05, 1.63100d+04 &
             , 0.0d+00    , 1.24416d+05,-4.97664d+05, 1.38240d+06, 1.89000d+05 &
             , 0.0d+00    , 0.0d+00    , 6.63552d+05,-1.43360d+06, 2.30000d+04 &
             , 0.0d+00    , 0.0d+00    , 0.0d+00    , 7.16800d+05, 2.21375d+05 &
             , 0.0d+00    , 0.0d+00    , 0.0d+00    , 0.0d+00    , 3.41550d+04 &
             /) / 5.5296d+05, (/ 5, 5 /))
    real(kind=PREC), dimension(6)  , parameter :: c = (/  9.361d+03  , 0.0d+00 &
             , 3.85d+04   , 2.0125d+04 , 0.0d+00    , 2.7648d+04 /)            &
             / 9.5634d+04
    real(kind=PREC), dimension(6)  , parameter :: d = (/ -1.40162d+05, 0.0d+00 &
             , 6.094d+05  ,-1.114925d+06,-6.30729d+05, 1.276416d+06 /)         &
             / 3.2643072d+07
    real(kind=PREC), parameter :: eps = tiny(qom)
!
!-------------------------------------------------------------------------------
!
    qom = params(1)
    vun = params(2)
    tmx = params(3)
    dti = params(4)
    dtm = params(5)
    tol = params(6)

    cnt        = 0
    t          = 0.0d+00
    dt         = min(dtm, dti)
    n          = 1
    si(1:6)    = state(1:6,n)
    state(7,n) = sqrt(dot_product(si(4:6), si(4:6)) + 1.0d+00)
    state(8,n) = dt
    n          = n + 1

    do while(t < tmx)

      do j = 1, 6
        s(j,1:6) = si(1:6)
        do i = 1, j - 1
          s(j,1:6) = s(j,1:6) + b(j-1,i) * k(i,1:6)
        end do
        call acceleration(qom, vun, dm, uu, bb, s(j,1:6), ac)
        k(j,1:6) = dt * ac(1:6)
      end do

      do i = 1, 6
        sf(i) = si(i) + sum(c(1:6) * k(1:6,i))
        er(i) =         sum(d(1:6) * k(1:6,i))
      end do

      sr(1:6) = abs(si(1:6)) + abs(k(1,1:6)) + eps

      del = maxval(abs(er / sr)) / tol

      if (del > 1.0d+00) then
        dtn = 9.0d-01 * dt * del**(-0.25d+00)
        dt  = min(dtm, max(abs(dtn), 1.0d-01 * dt))
      else
        cnt   = cnt + 1
        t     = t + dt
        call acceleration(qom, vun, dm, uu, bb, s(6,1:6), ac)
        k(6,1:6) = dt * ac(1:6)
        do while (t >= time(n) .and. n <= m)
          fr = (t - time(n)) / dt
          fl = 1.0d+00 - fr
          ds = sf(1:6) - si(1:6)
          state(1:6,n) = fr * si(1:6) + fl * sf(1:6)                           &
                       + fl * fr * ((k(1,1:6) - ds(1:6)) * fr                  &
                                  + (ds(1:6) - k(6,1:6)) * fl)
          state(  7,n) = sqrt(dot_product(state(4:6,n), state(4:6,n)) + 1.0d+00)
          state(  8,n) = dt
          n  = n + 1
        end do
        si(:) = sf(:)
        dtn = 9.0d-01 * dt * del**(-0.20d+00)
        dt  = min(dtm, dtn, 5.0d+00 * dt)
      end if
    end do

!-------------------------------------------------------------------------------
!
  end subroutine
!
!===============================================================================
!
! subroutine INTEGRATE_SI4:
! ------------------------
!
!   Subroutine integrate particle trajectory using implicit 4th order
!   Gauss-Legendre-RK method with adaptive time step.
!
!   Arguments:
!
!     m      - the number of steps to store (input);
!     dm     - domensions of the field components (input);
!     params - the vector of parameters (input);
!     time   - the vector of snapshot times (input);
!     uu, bb - the components of velocity and magnetic fields (input);
!     state  - the particle state vector, position, velocity, and others
!              (input/output);
!
!===============================================================================
!
  subroutine integrate_si4(m, dm, params, time, uu, bb, state, cnt)

    use accelerations, only : acceleration

    implicit none
    !$acc routine (integrate_si4) seq

! subroutine arguments
!
    integer                                        , intent(in)    :: m
    integer, dimension(3)                          , intent(in)    :: dm
    real(kind=PREC), dimension(8)                  , intent(in)    :: params
    real(kind=PREC), dimension(m)                  , intent(in)    :: time
    real(kind=PREC), dimension(dm(1),dm(2),dm(3),3), intent(in)    :: uu, bb
    real(kind=PREC), dimension(8,m)                , intent(inout) :: state
    integer(kind=8)                                , intent(out)   :: cnt

! local variables
!
    integer                         :: i, j, n
    real(kind=PREC)                 :: qom, vun, tmx, dti, dtm, tol
    real(kind=PREC)                 :: t, dt, dtn, del, fl, fr
    real(kind=PREC), dimension(6)   :: si, sf, er, sr, ds
    real(kind=PREC), dimension(6)   :: ac, ki, kf
    real(kind=PREC), dimension(6,6) :: s, k

    real(kind=PREC), dimension(2,6) :: z, zn

! parameters
!
    real(kind=PREC), parameter :: b1  = 5.0d-01, bh = 2.5d-01
    real(kind=PREC), parameter :: cs  = sqrt(3.0d+00)
    real(kind=PREC), parameter :: ch  = cs / 6.0d+00
    real(kind=PREC), dimension(2,2), parameter :: a = reshape(                 &
                                (/ bh, bh + ch, bh - ch, bh /), (/ 2, 2 /))
    real(kind=PREC), dimension(2)  , parameter :: c = (/  - cs, cs /)
    real(kind=PREC), dimension(2)  , parameter :: d =                          &
                        (/ 3.0d+00 + sqrt(2.7d+01), 3.0d+00 - sqrt(2.7d+01) /)
    real(kind=PREC), parameter :: eps = tiny(qom)

! ! local parameters
! !
!     real(kind=PREC), parameter :: c1  = b1 - ch, c2  = b1 + ch
!     real(kind=PREC), parameter :: a11 = bh, a12 = bh - ch, a21 = bh + ch, a22 = bh
!     real(kind=PREC), parameter :: d1   = - sqrt(3.0d+00), d2   = - d1
!
!-------------------------------------------------------------------------------
!
    qom = params(1)
    vun = params(2)
    tmx = params(3)
    dti = params(4)
    dtm = params(5)
    tol = params(6)

    cnt        = 0
    t          = 0.0d+00
    dt         = min(dtm, dti)
    n          = 1
    si(1:6)    = state(1:6,n)
    state(7,n) = sqrt(dot_product(si(4:6), si(4:6)) + 1.0d+00)
    state(8,n) = dt
    n          = n + 1

    call acceleration(qom, vun, dm, uu, bb, si(1:6), ac)
    z(1,1:6) = 4.0d-01 * dt * ac(1:6)
    z(2,1:6) = 6.0d-01 * dt * ac(1:6)

    do while(t < tmx)

! estimate intermediate states
!
      del = 1.0d+00
      do while(del > 1.0d-12)
        do j = 1, 2
          s(j,1:6) = si(1:6) + z(j,1:6)
          call acceleration(qom, vun, dm, uu, bb, s(j,1:6), ac)
          k(j,1:6) = dt * ac(1:6)
        end do
        do j = 1, 2
          zn(j,1:6) = a(j,1) * k(1,1:6) + a(j,2) * k(2,1:6)
        end do
        del = maxval(abs(zn(1:2,1:6) - z(1:2,1:6)))
        z(1:2,1:6) = zn(1:2,1:6)
      end do

      do i = 1, 6
        sf(i) = si(i) + sum(c(1:2) * z(1:2,i))
        er(i) =         sum(d(1:2) * z(1:2,i))
      end do

      sr(1:6) = abs(si(1:6)) + eps

      del = maxval(abs(er / sr)) / tol
!       print *, dt, del, er, sr

      if (del > 1.0d+00) then
        dtn = 9.0d-01 * dt * del**(-0.33d+00)
        dt  = min(dtm, max(abs(dtn), 1.0d-01 * dt))
      else
        cnt   = cnt + 1
        t     = t + dt
        call acceleration(qom, vun, dm, uu, bb, si(1:6), ac)
        ki(1:6) = dt * ac(1:6)
        call acceleration(qom, vun, dm, uu, bb, sf(1:6), ac)
        kf(1:6) = dt * ac(1:6)
        do while (t >= time(n) .and. n <= m)
          fr = (t - time(n)) / dt
          fl = 1.0d+00 - fr
          ds = sf(1:6) - si(1:6)
          state(1:6,n) = fr * si(1:6) + fl * sf(1:6)                           &
                       + fl * fr * ((ki(1:6) - ds(1:6)) * fr                   &
                                  + (ds(1:6) - kf(1:6)) * fl)
          state(  7,n) = sqrt(dot_product(state(4:6,n), state(4:6,n)) + 1.0d+00)
          state(  8,n) = dt
          n  = n + 1
        end do
        si(:) = sf(:)
        dtn = 9.0d-01 * dt * del**(-0.25d+00)
        dt  = min(dtm, dtn, 5.0d+00 * dt)
      end if
    end do

!-------------------------------------------------------------------------------
!
  end subroutine

!===============================================================================
!
end module
