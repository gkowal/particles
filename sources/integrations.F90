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
    real(kind=PREC), dimension(:)                  , intent(in)    :: params
    real(kind=PREC), dimension(m)                  , intent(in)    :: time
    real(kind=PREC), dimension(dm(1),dm(2),dm(3),3), intent(in)    :: uu, bb
    real(kind=PREC), dimension(8,m)                , intent(inout) :: state
    integer(kind=8), dimension(3)                  , intent(out)   :: cnt

! local variables
!
    integer                         :: i, j, n
    real(kind=PREC)                 :: qom, vun, tmx, dti, dtm, atol, rtol, fcmn, fcmx, safe
    real(kind=PREC)                 :: t, dt, dtn, del, fl, fr
    real(kind=PREC), dimension(6)   :: si, sf, er, sr, ds
    real(kind=PREC), dimension(6)   :: ac
    real(kind=PREC), dimension(6,6) :: s, k
    real(kind=PREC), dimension(3)   :: xl, xs

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
    qom   = params(1)
    vun   = params(2)
    tmx   = params(3)
    dti   = params(4)
    dtm   = params(5)
    atol  = params(6)
    rtol  = params(7)
    fcmn  = params(8)
    fcmx  = params(9)
    safe  = params(10)
    xl(1) = params(12)
    xl(2) = params(13)
    xl(3) = params(14)
    xs(1) = params(15)
    xs(2) = params(16)
    xs(3) = params(17)

    cnt        = 0
    t          = 0.0d+00
    dt         = dti_guess(5, dm, params, uu, bb, state(1:6,1))
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
        call acceleration(qom, vun, dm, xl, xs, uu, bb, s(j,1:6), ac)
        k(j,1:6) = dt * ac(1:6)
      end do

      do i = 1, 6
        sf(i) = si(i) + sum(c(1:6) * k(1:6,i))
        er(i) =         sum(d(1:6) * k(1:6,i))
      end do

! error estimator
!
      sr(1:6) = abs(si(1:6)) + abs(k(1,1:6)) + eps
      del = maxval(abs(er / sr)) / atol

      if (del <= 1.0d+00) then

        t = t + dt

! output
!
        call acceleration(qom, vun, dm, xl, xs, uu, bb, s(6,1:6), ac)
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

! increase counter of accepted steps
!
        cnt(2) = cnt(2) + 1

! new time step
!
        dtn = dt * min(fcmx, safe * del**(-0.20d+00))
      else

! increase counter of rejected steps
!
        cnt(3) = cnt(3) + 1

! new time step
!
        dtn = dt * max(fcmn, safe * del**(-0.25d+00))
      end if

! increase counter of all steps
!
      cnt(1) = cnt(1) + 1

! substitute time step
!
      dt = min(dtn, dtm)
    end do

!-------------------------------------------------------------------------------
!
  end subroutine
!
!===============================================================================
!
! subroutine INTEGRATE_DP853:
! --------------------------
!
!   Subroutine integrate particle trajectory using an explicit Runge-Kutta
!   method of order 8(5,3) due to Dorman & Prince with step control and dense
!   output.
!
!   Arguments:
!
!     m      - the number of steps to store (input);
!     dm     - dimensions of the field components (input);
!     params - the vector of parameters (input);
!     time   - the vector of snapshot times (input);
!     uu, bb - the components of velocity and magnetic fields (input);
!     state  - the particle state vector, position, velocity, and others
!              (input/output);
!     cnt    - the number of computed steps (output);
!
!===============================================================================
!
  subroutine integrate_dp853(m, dm, params, time, uu, bb, state, cnt)

    use accelerations, only : acceleration

    implicit none
    !$acc routine (integrate_dp853) seq

! subroutine arguments
!
    integer                                        , intent(in)    :: m
    integer, dimension(3)                          , intent(in)    :: dm
    real(kind=PREC), dimension(:)                  , intent(in)    :: params
    real(kind=PREC), dimension(m)                  , intent(in)    :: time
    real(kind=PREC), dimension(dm(1),dm(2),dm(3),3), intent(in)    :: uu, bb
    real(kind=PREC), dimension(8,m)                , intent(inout) :: state
    integer(kind=8), dimension(3)                  , intent(out)   :: cnt

! local variables
!
    logical                          :: rejected
    integer                          :: n
    real(kind=PREC)                  :: qom, vun, tmx, dti, dtm, atol, rtol, beta
    real(kind=PREC)                  :: fcmn, fcmx, safe, expo
    real(kind=PREC)                  :: t, dt, dtn, err3, err5, deno, err, errold
    real(kind=PREC)                  :: fl, fr
    real(kind=PREC), dimension(6)    :: si, ss, sf, er, sr, ds
    real(kind=PREC), dimension(10,6) :: k
    real(kind=PREC), dimension(3)    :: xl, xs

! parameters
!
    real(kind=PREC), parameter :: a0201 =  5.26001519587677318785587544488d-02 &
                                , a0301 =  1.97250569845378994544595329183d-02 &
                                , a0302 =  5.91751709536136983633785987549d-02 &
                                , a0401 =  2.95875854768068491816892993775d-02 &
                                , a0403 =  8.87627564304205475450678981324d-02 &
                                , a0501 =  2.41365134159266685502369798665d-01 &
                                , a0503 = -8.84549479328286085344864962717d-01 &
                                , a0504 =  9.24834003261792003115737966543d-01 &
                                , a0601 =  3.70370370370370370370370370370d-02 &
                                , a0604 =  1.70828608729473871279604482173d-01 &
                                , a0605 =  1.25467687566822425016691814123d-01 &
                                , a0701 =  3.71093750000000000000000000000d-02 &
                                , a0704 =  1.70252211019544039314978060272d-01 &
                                , a0705 =  6.02165389804559606850219397283d-02 &
                                , a0706 = -1.75781250000000000000000000000d-02 &
                                , a0801 =  3.70920001185047927108779319836d-02 &
                                , a0804 =  1.70383925712239993810214054705d-01 &
                                , a0805 =  1.07262030446373284651809199168d-01 &
                                , a0806 = -1.53194377486244017527936158236d-02 &
                                , a0807 =  8.27378916381402288758473766002d-03 &
                                , a0901 =  6.24110958716075717114429577812d-01 &
                                , a0904 = -3.36089262944694129406857109825d+00 &
                                , a0905 = -8.68219346841726006818189891453d-01 &
                                , a0906 =  2.75920996994467083049415600797d+01 &
                                , a0907 =  2.01540675504778934086186788979d+01 &
                                , a0908 = -4.34898841810699588477366255144d+01 &
                                , a1001 =  4.77662536438264365890433908527d-01 &
                                , a1004 = -2.48811461997166764192642586468d+00 &
                                , a1005 = -5.90290826836842996371446475743d-01 &
                                , a1006 =  2.12300514481811942347288949897d+01 &
                                , a1007 =  1.52792336328824235832596922938d+01 &
                                , a1008 = -3.32882109689848629194453265587d+01 &
                                , a1009 = -2.03312017085086261358222928593d-02 &
                                , a1101 = -9.37142430085987325717040216580d-01 &
                                , a1104 =  5.18637242884406370830023853209d+00 &
                                , a1105 =  1.09143734899672957818500254654d+00 &
                                , a1106 = -8.14978701074692612513997267357d+00 &
                                , a1107 = -1.85200656599969598641566180701d+01 &
                                , a1108 =  2.27394870993505042818970056734d+01 &
                                , a1109 =  2.49360555267965238987089396762d+00 &
                                , a1110 = -3.04676447189821950038236690220d+00 &
                                , a1201 =  2.27331014751653820792359768449d+00 &
                                , a1204 = -1.05344954667372501984066689879d+01 &
                                , a1205 = -2.00087205822486249909675718444d+00 &
                                , a1206 = -1.79589318631187989172765950534d+01 &
                                , a1207 =  2.79488845294199600508499808837d+01 &
                                , a1208 = -2.85899827713502369474065508674d+00 &
                                , a1209 = -8.87285693353062954433549289258d+00 &
                                , a1210 =  1.23605671757943030647266201528d+01 &
                                , a1211 =  6.43392746015763530355970484046d-01

    real(kind=PREC), parameter :: b01 =  5.42937341165687622380535766363d-02   &
                                , b06 =  4.45031289275240888144113950566d+00   &
                                , b07 =  1.89151789931450038304281599044d+00   &
                                , b08 = -5.80120396001058478146721142270d+00   &
                                , b09 =  3.11164366957819894408916062370d-01   &
                                , b10 = -1.52160949662516078556178806805d-01   &
                                , b11 =  2.01365400804030348374776537501d-01   &
                                , b12 =  4.47106157277725905176885569043d-02

    real(kind=PREC), parameter :: bh1 = 0.244094488188976377952755905512d+00   &
                                , bh2 = 0.733846688281611857341361741547d+00   &
                                , bh3 = 0.220588235294117647058823529412d-01

    real(kind=PREC), parameter :: er01 =  0.1312004499419488073250102996d-01   &
                                , er06 = -0.1225156446376204440720569753d+01   &
                                , er07 = -0.4957589496572501915214079952d+00   &
                                , er08 =  0.1664377182454986536961530415d+01   &
                                , er09 = -0.3503288487499736816886487290d+00   &
                                , er10 =  0.3341791187130174790297318841d+00   &
                                , er11 =  0.8192320648511571246570742613d-01   &
                                , er12 = -0.2235530786388629525884427845d-01
!
!-------------------------------------------------------------------------------
!
    qom   = params(1)
    vun   = params(2)
    tmx   = params(3)
    dti   = params(4)
    dtm   = params(5)
    atol  = params(6)
    rtol  = params(7)
    fcmn  = params(8)
    fcmx  = params(9)
    safe  = params(10)
    beta  = params(11)
    xl(1) = params(12)
    xl(2) = params(13)
    xl(3) = params(14)
    xs(1) = params(15)
    xs(2) = params(16)
    xs(3) = params(17)

    expo       = 2.0d-01 * beta - 1.25d-01
    errold     = 1.0d-04
    rejected   = .false.

    cnt        = 0
    t          = 0.0d+00
    dt         = dti_guess(8, dm, params, uu, bb, state(1:6,1))
    n          = 1
    si(1:6)    = state(1:6,n)
    state(7,n) = sqrt(dot_product(si(4:6), si(4:6)) + 1.0d+00)
    state(8,n) = dt
    n          = n + 1

    do while(t < tmx)

! the 12 steps
!
      ss(1:6)  = si(1:6)
      call acceleration(qom, vun, dm, xl, xs, uu, bb, ss, k(1,1:6))

      ss(1:6)  = si(1:6) + dt * a0201 * k(1,1:6)
      call acceleration(qom, vun, dm, xl, xs, uu, bb, ss, k(2,1:6))

      ss(1:6)  = si(1:6) + dt * (a0301 * k(1,1:6) + a0302 * k(2,1:6))
      call acceleration(qom, vun, dm, xl, xs, uu, bb, ss, k(3,1:6))

      ss(1:6)  = si(1:6) + dt * (a0401 * k(1,1:6) + a0403 * k(3,1:6))
      call acceleration(qom, vun, dm, xl, xs, uu, bb, ss, k(4,1:6))

      ss(1:6)  = si(1:6) + dt * (a0501 * k(1,1:6) + a0503 * k(3,1:6)           &
                               + a0504 * k(4,1:6))
      call acceleration(qom, vun, dm, xl, xs, uu, bb, ss, k(5,1:6))

      ss(1:6)  = si(1:6) + dt * (a0601 * k(1,1:6) + a0604 * k(4,1:6)           &
                               + a0605 * k(5,1:6))
      call acceleration(qom, vun, dm, xl, xs, uu, bb, ss, k(6,1:6))

      ss(1:6)  = si(1:6) + dt * (a0701 * k(1,1:6) + a0704 * k(4,1:6)           &
                               + a0705 * k(5,1:6) + a0706 * k(6,1:6))
      call acceleration(qom, vun, dm, xl, xs, uu, bb, ss, k(7,1:6))

      ss(1:6)  = si(1:6) + dt * (a0801 * k(1,1:6) + a0804 * k(4,1:6)           &
                               + a0805 * k(5,1:6) + a0806 * k(6,1:6)           &
                               + a0807 * k(7,1:6))
      call acceleration(qom, vun, dm, xl, xs, uu, bb, ss, k(8,1:6))

      ss(1:6)  = si(1:6) + dt * (a0901 * k(1,1:6) + a0904 * k(4,1:6)           &
                               + a0905 * k(5,1:6) + a0906 * k(6,1:6)           &
                               + a0907 * k(7,1:6) + a0908 * k(8,1:6))
      call acceleration(qom, vun, dm, xl, xs, uu, bb, ss, k(9,1:6))

      ss(1:6)  = si(1:6) + dt * (a1001 * k(1,1:6) + a1004 * k(4,1:6)           &
                               + a1005 * k(5,1:6) + a1006 * k(6,1:6)           &
                               + a1007 * k(7,1:6) + a1008 * k(8,1:6)           &
                               + a1009 * k(9,1:6))
      call acceleration(qom, vun, dm, xl, xs, uu, bb, ss, k(10,1:6))

      ss(1:6)  = si(1:6) + dt * (a1101 * k(1,1:6) + a1104 * k(4,1:6)           &
                               + a1105 * k(5,1:6) + a1106 * k(6,1:6)           &
                               + a1107 * k(7,1:6) + a1108 * k(8,1:6)           &
                               + a1109 * k(9,1:6) + a1110 * k(10,1:6))
      call acceleration(qom, vun, dm, xl, xs, uu, bb, ss, k(2,1:6))

      ss(1:6)  = si(1:6) + dt * (a1201 * k(1,1:6) + a1204 * k(4,1:6)           &
                               + a1205 * k(5,1:6) + a1206 * k(6,1:6)           &
                               + a1207 * k(7,1:6) + a1208 * k(8,1:6)           &
                               + a1209 * k(9,1:6) + a1210 * k(10,1:6)          &
                               + a1211 * k(2,1:6))
      call acceleration(qom, vun, dm, xl, xs, uu, bb, ss, k(3,1:6))

      k(4,1:6) = b01 * k( 1,1:6) + b06 * k( 6,1:6) + b07 * k( 7,1:6)           &
               + b08 * k( 8,1:6) + b09 * k( 9,1:6) + b10 * k(10,1:6)           &
               + b11 * k( 2,1:6) + b12 * k( 3,1:6)
      sf(1:6)  = si(1:6) + dt * k(4,1:6)

! error estimation
!
      sr(1:6) = atol + rtol * max(abs(si(1:6)), abs(sf(1:6)))
      er(1:6) = k(4,1:6) - bh1 * k(1,1:6) - bh2 * k(9,1:6) - bh3 * k(3,1:6)
      err3 = sum((er(:) / sr(:))**2) ! 3rd order error
      er(1:6) = er01 * k(1,1:6) + er06 * k(6,1:6) + er07 * k( 7,1:6)           &
              + er08 * k(8,1:6) + er09 * k(9,1:6) + er10 * k(10,1:6)           &
              + er11 * k(2,1:6) + er12 * k(3,1:6)
      err5 = sum((er(:) / sr(:))**2) ! 5th order error
      deno = err5 + 1.0d-02 * err3
      if (deno <= 0.0d+00) then
        err = abs(dt) * err5 / sqrt(6.0d+00)
      else
        err = abs(dt) * err5 / sqrt(6.0d+00 * deno)
      end if

      if (err <= 1.0d+00) then

        t = t + dt

! output
!
        call acceleration(qom, vun, dm, xl, xs, uu, bb, si(1:6), k(1,1:6))
        call acceleration(qom, vun, dm, xl, xs, uu, bb, sf(1:6), k(2,1:6))
        do while (t >= time(n) .and. n <= m)
          fr = (t - time(n)) / dt
          fl = 1.0d+00 - fr
          ds = sf(1:6) - si(1:6)
          state(1:6,n) = fr * si(1:6) + fl * sf(1:6)                           &
                       + fl * fr * ((dt*k(1,1:6) - ds(1:6)) * fr               &
                                  + (ds(1:6) - dt*k(2,1:6)) * fl)
          state(  7,n) = sqrt(dot_product(state(4:6,n), state(4:6,n)) + 1.0d+00)
          state(  8,n) = dt
          n  = n + 1
        end do
        si(:) = sf(:)

! increase counter of accepted steps
!
        cnt(2) = cnt(2) + 1

! new time step
!
        dtn   = max(fcmn, safe * err**expo * errold**beta)
        if (rejected) then
          dtn = dt * min(dtn, 1.0d+00)
        else
          dtn = dt * min(dtn, fcmx)
        end if
        errold = max(err, 1.0d-04)

! set rejected flag
!
        rejected = .false.

      else

! increase counter of rejected steps
!
        cnt(3) = cnt(3) + 1

! new time step
!
        dtn = dt * max(fcmn, safe * err**expo)

! set rejected flag
!
        rejected = .true.

      end if

! increase counter of all steps
!
      cnt(1) = cnt(1) + 1

! substitute time step
!
      dt = min(dtn, dtm)
    end do

!-------------------------------------------------------------------------------
!
  end subroutine
!
!===============================================================================
!
! function DTI_GUESS:
! ------------------
!
!   Function estimates the initial time step for a given integration order.
!
!   Arguments:
!
!     iord   - the order of integration (input);
!     dm     - dimensions of the field components (input);
!     params - the vector of parameters (input);
!     uu, bb - the components of velocity and magnetic fields (input);
!     si     - the initial particle state vector: position and velocity (input);
!
!===============================================================================
!
  function dti_guess(iord, dm, params, uu, bb, si) result(dti)

    use accelerations, only : acceleration

    implicit none
    !$acc routine (dti_guess) seq

! subroutine arguments
!
    integer                                        , intent(in) :: iord
    integer, dimension(3)                          , intent(in) :: dm
    real(kind=PREC), dimension(16)                 , intent(in) :: params
    real(kind=PREC), dimension(dm(1),dm(2),dm(3),3), intent(in) :: uu, bb
    real(kind=PREC), dimension(6)                  , intent(in) :: si

! local variables
!
    real(kind=PREC)                  :: qom, vun, dti, dtf, dtm, atol, rtol
    real(kind=PREC)                  :: d0, d1, d2
    real(kind=PREC), dimension(6)    :: sr, sf, ki, kf
    real(kind=PREC), dimension(3)    :: xl, xs

! parameters
!
!
!-------------------------------------------------------------------------------
!
    qom   = params(1)
    vun   = params(2)
    dtm   = params(5)
    atol  = params(6)
    rtol  = params(7)
    xl(1) = params(12)
    xl(2) = params(13)
    xl(3) = params(14)
    xs(1) = params(15)
    xs(2) = params(16)
    xs(3) = params(17)

    call acceleration(qom, vun, dm, xl, xs, uu, bb, si(1:6), ki(1:6))

    sr(1:6) = atol + rtol * abs(si(1:6))

    d0  = sqrt(sum((si(1:6) / sr(1:6))**2) / 6.0d+00)
    d1  = sqrt(sum((ki(1:6) / sr(1:6))**2) / 6.0d+00)
    if (min(d0, d1) < 1.0d-10) then
      dti = 1.0d-06
    else
      dti = 1.0d-01 * (d0 / d1)
    end if

    sf(1:6) = si(1:6) + dti * ki(1:6)
    call acceleration(qom, vun, dm, xl, xs, uu, bb, sf(1:6), kf(1:6))
    d2  = max(d1, sqrt(sum(((sf(1:6) - si(1:6)) / sr(1:6))**2) / 6.0d+00) / dti)
    if (d2 <= 1.0d-15) then
      dtf = max(1.0d-06, 1.0d-03 * dti)
    else
      dtf = (1.0d-02 / d2)**(1.0d+00 / iord)
    end if
    dti = min(1.0d+02 * dti, dtf, dtm)

    return

!-------------------------------------------------------------------------------
!
  end function

!===============================================================================
!
end module
