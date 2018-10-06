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
!! program: PARTICLES
!!
!!  This is the main subroutine of the PARTICLES code.
!!
!!******************************************************************************
!
program particles

  use fields      , only : initialize_fields, finalize_fields
  use fields      , only : dm, uu, bb
  use fitsio      , only : initialize_fitsio, finalize_fitsio, write_data
  use integrations, only : integrate_rk5
  use parameters  , only : read_parameters
  use parameters  , only : get_parameter_integer, get_parameter_real, get_parameter_string
  use random      , only : initialize_random, finalize_random                  &
                         , uniform, normal, moments

  implicit none

  real(kind=PREC), parameter :: pi = 3.141592653589793d+00
  real(kind=PREC), parameter :: c  = 2.99792458d+08
  real(kind=PREC), parameter :: kb = 1.38064852d-23
  real(kind=PREC), parameter :: mp = 1.672621898d-27
  real(kind=PREC), parameter :: me = 9.10938356d-31
  real(kind=PREC), parameter :: e  = 1.6021766208d-19
  real(kind=PREC), parameter :: fc = (1.0d+00 * PREC) / 1024**3

  character(len=32) :: simulation = 'newtonian'
  character(len=32) :: ptype      = 'proton'
  character(len=32) :: tscale     = 'linear'
  character(len=32) :: method     = 'rk5'
  integer           :: nparticles = 1000
  integer           :: nsteps     = 100
  real(kind=PREC)   :: dtini      = 1.0d-04
  real(kind=PREC)   :: dtmax      = 1.0d-03
  real(kind=PREC)   :: tmin       = 1.0d-03
  real(kind=PREC)   :: tmax       = 1.0d+00
  real(kind=PREC)   :: tolerance  = 1.0d-06
  real(kind=PREC)   :: temp       = 1.0d+03
  real(kind=PREC)   :: vth        = 1.0d-03
  real(kind=PREC)   :: valf       = 1.0d-03
  real(kind=PREC)   :: bfac       = 1.1209982432795858d+01
  real(kind=PREC)   :: dens       = 1.0d+00
  real(kind=PREC)   :: bunit      = 1.0d+00
  real(kind=PREC)   :: tunit      = 1.0d+00
  real(kind=PREC)   :: lunit      = 1.0d+00
  real(kind=PREC)   :: lscale     = 1.0d+00
  real(kind=PREC)   :: gamma      = 5.0d+00 / 3.0d+00
  real(kind=PREC)   :: qom        = 9.5788332241480675d+03

  character(len=4), dimension(6) :: label = (/ 'posx', 'posy', 'posz',         &
                                               'momx', 'momy', 'momz' /)

  character(len=32) :: fname, stmp
  integer           :: i, iret
  integer(kind=8)   :: n
  integer(kind=8)   :: t1, t2, dt, count_rate, count_max
  real              :: secs
  real(kind=PREC)   :: ltmn, ltmx
  real(kind=PREC), dimension(8)                  :: params
  real(kind=PREC), dimension(4)                  :: moms

! OpenMP variables
!
!$  integer :: np, npar, omp_get_num_threads, omp_get_thread_num

! allocatable variables
!
  integer(kind=8), dimension(:)    , allocatable :: counter
  real(kind=PREC), dimension(:)    , allocatable :: time
  real(kind=PREC), dimension(:,:,:), allocatable :: state
!
!-------------------------------------------------------------------------------
!
!$omp parallel
!$  npar = omp_get_num_threads()
!$omp end parallel
    write (*,"(80('-'))")
    write (*,"(19('='),7x,'Particles algorithm started',8x,19('='))")
    write (*,"(19('='),5x,'Copyright (C) 2018 Grzegorz Kowal',4x,19('='))")
#ifdef _OPENACC
    write (*,"(19('='),9x,'OpenACC support enabled',10x,19('='))")
#endif /* _OPENACC */
!$  write (*,"(19('='),10x,'OpenMP:',i5,1x,'threads  ',10x,19('='))") npar

  iret = 0
  call read_parameters(iret)

  call get_parameter_string ("simulation"   , simulation)
  call get_parameter_string ("method"       , method    )
  call get_parameter_string ("particle_type", ptype     )
  call get_parameter_string ("time_scale"   , tscale    )
  call get_parameter_integer("nparticles"   , nparticles)
  call get_parameter_integer("nsteps"       , nsteps    )
  call get_parameter_real   ("bunit"        , bunit     )
  call get_parameter_real   ("tunit"        , tunit     )
  call get_parameter_real   ("temperature"  , temp      )
  call get_parameter_real   ('density'      , dens      )
  call get_parameter_real   ("tmin"         , tmin      )
  call get_parameter_real   ("tmax"         , tmax      )
  call get_parameter_real   ("dtini"        , dtini     )
  call get_parameter_real   ("dtmax"        , dtmax     )
  call get_parameter_real   ("tolerance"    , tolerance )

  valf = 7.2757116026881272d+00 * bunit / sqrt(dens)
  vth  = sqrt(kb * temp / mp) / c

  write(*,*)
  write(*,"(' Fields:')")
  select case(trim(simulation))
  case('relativistic', 'rel', 'r', 'srmhd', 'srhd')
    write(*,"('   simulation    =  relativistic')")
    lscale = 1.0d+00
    lunit  = tunit * c
  case default
    write(*,"('   simulation    =  newtonian')")
    lscale = 1.0d+00 / valf
    lunit  = tunit * valf * c
  end select
  write(*,"('   B unit        = ', 1es22.15, ' [Gs]')"    ) bunit
  write(*,"('   length unit   = ', 1es22.15, ' [m]')"     ) lunit
  write(*,"('   length scale  = ', 1es22.15, ' [1/L]')"   ) lscale
  write(*,"('   density       = ', 1es22.15, ' [mₚ/cm³]')") dens
  write(*,"('   Alfvén speed  = ', 1es22.15, ' [c]')"     ) valf

  write(*,*)
  write(*,"(' Particles:')")
  write(stmp,"(i12)") nparticles
  write(*,"('   nparticles    = ', 1x, a)") trim(adjustl(stmp))
  select case(trim(ptype))
  case('proton', 'p')
    write(*,"('   particle type =  proton')")
    qom =   1.0d-04 * e * bunit * tunit / mp
  case('electron', 'e')
    write(*,"('   particle type =  electron')")
    qom = - 1.0d-04 * e * bunit * tunit / me
  case default
    write(*,"('   particle type =  unknown')")
  end select
  write(*,"('   length unit   = ', 1es22.15, ' [m]')"     ) c * tunit
  write(*,"('   time unit     = ', 1es22.15, ' [s]')"     ) tunit
  write(*,"('   q/m           = ', 1es22.15, ' [1/s]')"   ) qom
  write(*,"('   temperature   = ', 1es22.15, ' [K]')"     ) temp
  write(*,"('   thermal speed = ', 1es22.15, ' [c]')"     ) vth
  qom = qom * bfac

  write(*,*)
  write(*,"(' Methods:')")
  select case(trim(method))
  case('rk5')
    write(*,"('   integration   =  5th order Runge-Kutta-Cash-Karp')")
  end select
  write(*,"('   dtini         = ', 1es12.5)") dtini
  write(*,"('   dtmax         = ', 1es12.5)") dtmax
  write(*,"('   tolerance     = ', 1es12.5)") tolerance

  write(*,*)
  write(*,"(' Output:')")
  select case(trim(tscale))
  case('lin', 'linear')
    write(*,"('   time scale    =  linear')")
    nsteps = nsteps + 1
  case('log', 'logarithmic', 'log10')
    write(*,"('   time scale    =  logarithmic')")
    nsteps = int(log10(tmax) - log10(tmin)) * nsteps + 2
  case default
    write(*,"('   time scale    =  unknown')")
  end select
  write(stmp,"(i12)") nsteps
  write(*,"('   nsteps        = ', 1x, a)"  ) trim(adjustl(stmp))
  write(*,"('   tmin          = ', 1es12.5)") tmin
  write(*,"('   tmax          = ', 1es12.5)") tmax

  call system_clock(count_max=count_max, count_rate=count_rate)

  call initialize_random()
  call initialize_fitsio()
  call initialize_fields()

  allocate(counter(nparticles))
  allocate(time(nsteps))
  allocate(state(8,nsteps,nparticles))

  counter(:)   = 0
  time(:)      = 0.0d+00
  state(:,:,:) = 0.0d+00

  select case(trim(tscale))
  case('lin', 'linear')
    do n = 2, nsteps
      time(n) = (tmax * (n - 1)) / (nsteps - 1)
    end do
  case('log', 'logarithmic', 'log10')
    ltmn = dlog10(tmin)
    ltmx = dlog10(tmax)
    do n = 2, nsteps
      time(n) = 1.0d+01**(((nsteps - n) * ltmn + (n - 2) * ltmx) / (nsteps - 2))
    end do
  end select

  write(*,*)
  write(*,"(' Memory requirements:')")
  write(stmp,"(1f12.3)") fc * (size(uu) + size(bb))
  write(*,"('   Plasma data   =  ', a, ' GB')") trim(adjustl(stmp))
  write(stmp,"(1f12.3)") fc * size(state)
  write(*,"('   Particle data =  ', a, ' GB')") trim(adjustl(stmp))

! generate initial positions and velocity components
!
  call system_clock(t1)
  call uniform(1.0d+00, state(1:3,1,1:nparticles))
  call normal(0.0d+00, vth, state(4:6,1,1:nparticles))
  call system_clock(t2)
  dt = t2 - t1
  secs = real(dt) / real(count_rate)
  write(*,*)
  write(stmp,"(1f16.4)") secs
  write(*,"(' Initial state generated in ',a,' seconds.')") trim(adjustl(stmp))

! check statistical moments of velocity components
!
  write(*,*)
  write(*,"(' Statistics of the initial state:')")
  write(*,"(6x,4a24)") 'mean', 'standard deviation', 'skewness', 'kurtosis'
  do n = 1, 6
    call moments(state(n,1,:), moms(1:4))
    write(*,"(3x,a,': ',4es24.15)") label(n), moms(1:4)
  end do

  write(*,*)
  write(*,"(' Integrating trajectories...')")

  params(1) = qom
  params(2) = valf
  params(3) = tmax
  params(4) = dtini
  params(5) = dtmax
  params(6) = tolerance

  call system_clock(t1)

  select case(trim(method))
  case('rk5')
!$omp parallel do
!$acc data copy(state) copyin(dm,uu,bb,time)
!$acc kernels
    do n = 1, nparticles
      call integrate_rk5(nsteps, dm, params, time, uu, bb, state(:,:,n), counter(n))
    end do
!$acc end kernels
!$acc end data
  end select

  call system_clock(t2)

  dt = t2 - t1
  secs = real(dt) / real(count_rate)
  write(stmp,"(1f16.4)") secs
  write(*,"(' Integration completed in ',a,' seconds.')") trim(adjustl(stmp))

  write(*,*)
  write(*,"(' Storing data...')")

  call system_clock(t1)
  call write_data('!counter.fits', counter)
  call write_data('!time.fits'   , time , tunit = tunit)
  call write_data('!state.fits'  , state, temperature = temp)
  call system_clock(t2)

  dt = t2 - t1
  secs = real(dt) / real(count_rate)
  write(stmp,"(1f16.4)") secs
  write(*,"(' Data stored in ',a,' seconds.')") trim(adjustl(stmp))

  deallocate(counter)
  deallocate(time)
  deallocate(state)

  call finalize_fields()
  call finalize_fitsio()
  call finalize_random()

end program
