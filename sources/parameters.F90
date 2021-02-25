!!******************************************************************************
!!
!!  This file is part of the PARTICLES source code, a program to integrate
!!  test particles in magnetohydrodynamical simulations.
!!
!!  Copyright (C) 2008-2021 Grzegorz Kowal <grzegorz@amuncode.org>
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
!! module: PARAMETERS
!!
!!  This module handles runtime parameters by reading them from a parameter
!!  file and distributing among all processes.
!!
!!******************************************************************************
!
module parameters

! module variables are not implicit by default
!
  implicit none

! MODULE INTERFACES:
! =================
!
  interface get_parameter
    module procedure get_parameter_integer
    module procedure get_parameter_real
    module procedure get_parameter_string
  end interface

! MODULE PARAMETERS:
! =================
!
! module parameters determining the name and value field lengths, and the
! maximum string length
!
  integer, parameter        :: nlen = 64, vlen = 128, mlen = 256

! the name of the parameter file
!
  character(len=mlen), save :: fname   = './params.in'

! a file handler to the parameter file
!
  integer(kind=4)    , save :: punit   = 10

! the number of parameters stored in the parameter file
!
  integer            , save :: nparams = 0

! allocatable arrays to store parameter names and values
!
  character(len=nlen), dimension(:), allocatable, save :: pnames
  character(len=vlen), dimension(:), allocatable, save :: pvalues

! by default everything is private
!
  private

! declare public subroutines
!
  public :: read_parameters, finalize_parameters
  public :: get_parameter_file, get_parameter

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
  contains
!
!===============================================================================
!
! subroutine READ_PARAMETERS:
! --------------------------
!
!   Subroutine checks if the parameter file exists, checks its length,
!   allocates structures to store all parameters provided in the parameters
!   file and read parameters into these structures.
!
!   Arguments:
!
!     verbose - the flag determining if the subroutine should be verbose;
!     status  - the return flag of the procedure execution status;
!
!   Note:
!
!     There is a possibility to specify customized input file by adding a
!     command line option -i or --input followed by the name of the input file.
!
!===============================================================================
!
  subroutine read_parameters(verbose, status)

! include external procedures and variables
!
    use iso_fortran_env, only : error_unit

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    logical, intent(in)  :: verbose
    integer, intent(out) :: status

! local variables
!
    character(len=mlen) :: opt, arg
    integer             :: l
    logical             :: info
!
!-------------------------------------------------------------------------------
!
    status = 0

! parse the command line to check if a different parameter file has been
! provided
!
    do l = 1, command_argument_count()
      call get_command_argument(l, arg)
      if (trim(arg) == '-i' .or. trim(arg) == '--input') then
        opt = trim(arg)
        call get_command_argument(l + 1, arg)
        if (trim(arg) /= '') then
          fname = trim(arg)
        else
          if (verbose) then
            write(error_unit,*) "The option '" // trim(opt) //                 &
                              "' requires an argument! Exiting..."
          end if
          status = 112
          return
        end if
      end if
    end do

! check if the file exists
!
    inquire(file = fname, exist = info)

    if (info) then

! obtain the number of parameters stored in the file
!
      call get_parameters_number(status)

      if (status > 0) return

! if the parameter file is empty, print an error and quit the subroutine
!
      if (nparams <= 0) then
        write(error_unit,*) "The parameter file '" // trim(fname)              &
                                                   // "' is empty! Exiting..."
        status = 110

        return
      end if

! allocate arrays to store the parameter names and values as string variables
!
      allocate(pnames (nparams))
      allocate(pvalues(nparams))

! get the parameter names and values and copy them to the corresponding arrays
!
      call get_parameters(status)

    else

      write(error_unit,*) "The parameter file '" // trim(fname)                &
                                             // "' does not exist! Exiting..."
      status = 111

    end if

!-------------------------------------------------------------------------------
!
  end subroutine read_parameters
!
!===============================================================================
!
! subroutine FINALIZE_PARAMETERS:
! ------------------------------
!
!   Subroutine releases memory used by arrays in this module.
!
!
!===============================================================================
!
  subroutine finalize_parameters()

! local variables are not implicit by default
!
    implicit none
!
!-------------------------------------------------------------------------------
!
    if (allocated(pnames) ) deallocate(pnames)
    if (allocated(pvalues)) deallocate(pvalues)

!-------------------------------------------------------------------------------
!
  end subroutine finalize_parameters
!
!===============================================================================
!
! subroutine GET_PARAMETER_FILE:
! -----------------------------
!
!   Subroutine returns the full path to the parameter file.
!
!   Arguments:
!
!     pfile  - the parameter full file path;
!     status - the status value, 0 for success;
!
!===============================================================================
!
  subroutine get_parameter_file(pfile, status)

! import external procedures
!
    use iso_fortran_env, only : error_unit

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(out) :: pfile
    integer         , intent(out) :: status

! local variables
!
    character(len=mlen) :: tfile

! local parameters
!
    character(len=*), parameter :: loc = 'PARAMETERS::get_parameter_file()'
!
!-------------------------------------------------------------------------------
!
    status = 0
    if (len(pfile) <= mlen) then
      write(pfile,"(a)") trim(adjustl(fname))
    else
      write(error_unit,"('[',a,']: ',a)") trim(loc)                            &
                      , "Parameter file path too long for subroutine argument!"
      write(tfile,"(a)") trim(adjustl(fname))
      write(pfile,"(a)") tfile(1:len(pfile))
      status = 1
    end if

!-------------------------------------------------------------------------------
!
  end subroutine get_parameter_file
!
!===============================================================================
!
! subroutine GET_PARAMETERS_NUMBER:
! --------------------------------
!
!   Subroutine scans the input file and accounts the number of parameters
!   stored in it.
!
!   Arguments:
!
!     status - the status value, 0 for success;
!
!===============================================================================
!
  subroutine get_parameters_number(status)

! import external procedures
!
    use iso_fortran_env, only : error_unit

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(out) :: status

! local variable to store the line content
!
    character(len=mlen) :: line

! local parameters
!
    character(len=*), parameter :: loc = 'PARAMETERS::get_parameters_number()'
!
!-------------------------------------------------------------------------------
!
! reset the number of parameters
!
    nparams = 0

! open the parameter file
!
    open(newunit = punit, file = fname, err = 30)

! read the line
!
10  read(punit, fmt = "(a)", end = 20) line

! if the line is empty or it's a comment, skip the counting
!
    if ((len_trim(line) == 0)                                                  &
                           .or. index(trim(adjustl(line)), '#') == 1) go to 10

! increase the number of parameters
!
    nparams = nparams + 1

! go to the next line
!
    go to 10

! close the file
!
20  close(punit)

! quit the subroutine without printing any errors since everything went fine
!
    return

! print a massage if an error occurred
!
30  write(error_unit,"('[',a,']: ',a)") trim(loc)                              &
                    , "Cannot open the parameter file '" // trim(fname) // "'!"

! set the status flag
!
    status = 112

!-------------------------------------------------------------------------------
!
  end subroutine get_parameters_number
!
!===============================================================================
!
! subroutine GET_PARAMETERS:
! -------------------------
!
!   Subroutine scans the input file, reads parameter names and values, and
!   stores them in module arrays.
!
!   Arguments:
!
!     status - the status value, 0 for success;
!
!===============================================================================
!
  subroutine get_parameters(status)

! import external procedures
!
    use iso_fortran_env, only : error_unit

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    integer, intent(out) :: status

! the parameter counter
!
    integer :: np, nl

! local variables to store the line content, the parameter name and value
!
    character(len=256) :: line, name, value

! local parameters
!
    character(len=*), parameter :: loc = 'PARAMETERS::get_parameters_number()'
!
!-------------------------------------------------------------------------------
!
! initialize the parameter counter
!
    np = 1
    nl = 0

! open the parameter file
!
    open(newunit = punit, file = fname, err = 30)

! read the line
!
10  read(punit, fmt = "(a)", end = 20) line

! increase the line number
!
    nl = nl + 1

! if the line is empty or it's a comment, skip the counting
!
    if ((len_trim(line) == 0)                                                  &
                           .or. index(trim(adjustl(line)), '#') == 1) go to 10

! parse the line to get parameter name and value
!
    call parse_line(line, name, value, status)

! check if the line was parsed successfuly
!
    if (status > 0) then
      write(error_unit,"('[',a,']: ',a)") trim(loc)                            &
                      , "Wrong parameter format in '"                          &
                                       // trim(adjustl(fname)) // "'."
      write (error_unit,"('[',a,']: Line',i4,' : ',a)")                        &
                                          trim(loc), nl, trim(line)
      go to 30
    end if

! fill the arrays of parameter names and values
!
    pnames (np) = name (1:nlen)
    pvalues(np) = value(1:vlen)

! increase the parameter counter
!
    np = np + 1

! go to the next line
!
    go to 10

! close the file
!
20  close(punit)

! quit the subroutine without printing any errors since everything went fine
!
    return

! print a massage if an error occurred
!
30  write(error_unit,"('[',a,']: ',a)") trim(loc)                              &
                    , "Cannot open the parameter file '" // trim(fname) // "'!"

! set the status flag
!
    status = 140

!-------------------------------------------------------------------------------
!
  end subroutine get_parameters
!
!===============================================================================
!
! subroutine PARSE_LINE:
! ---------------------
!
!   Subroutine extracts the parameter name and value from the input line.
!
!   Arguments:
!
!     line   - the input line containing the parameter information;
!     name   - the extracted name of the parameter;
!     value  - the extracted value of the parameter;
!     status - the status value, 0 for success;
!
!===============================================================================
!
  subroutine parse_line(line, name, value, status)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(in)    :: line
    character(len=*), intent(inout) :: name, value
    integer         , intent(out)   :: status

! local indices to store positions in the input line
!
    integer :: l, p, c, i, j = 1
!
!-------------------------------------------------------------------------------
!
! reset the status flag
!
    status = 0

! get the length of line
!
    l = len_trim(line)

! find the indices of '=' and '#' in the line
!
    p = index(line, '=')
    c = index(line, '#')
    i = index(line, '"')
    if (i > 0) then
      j = index(line, '"', back = .true.)
    else
      i = index(line, "'")
      if (i > 0) then
        j = index(line, "'", back = .true.)
      end if
    end if

! remove the length of the in-line comment from the length of line
!
    if (c > 0) l = c - 1
    if (i > 0 .and. j > 0) then
     i = i + 1
     j = j - 1
    else
     i = p + 1
     j = l
    end if

! limit the indices, so we don't overrun the variable memory
!
    p = min(p, nlen)
    j = min(j, vlen + i)

! extract the parameter name
!
    name  = trim(adjustl(line(1:p-1)))

! extract the parameter value
!
    value = trim(adjustl(line(i:j)))

! check possible errors in formatting
!
    if (p <= 2 .or. len_trim(name) == 0 .or. len_trim(value) == 0) status = 1

!-------------------------------------------------------------------------------
!
  end subroutine parse_line
!
!===============================================================================
!
! subroutine GET_PARAMETER_INTEGER:
! --------------------------------
!
!   Subroutine reads a given parameter name and returns its integer value.
!
!   Arguments:
!
!     name  - the input parameter name;
!     value - the output integer value of parameter;
!
!===============================================================================
!
  subroutine get_parameter_integer(name, value)

! import external procedures
!
    use iso_fortran_env, only : error_unit

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(in)    :: name
    integer         , intent(inout) :: value

! local parameter counter
!
    integer :: np

! local parameters
!
    character(len=*), parameter :: loc = 'PARAMETERS::get_parameter_integer()'
!
!-------------------------------------------------------------------------------
!
! find the selected parameter
!
    np = 1
    do while (np <= nparams)
      if (name == pnames(np)) then
        read(pvalues(np), err = 100, fmt = *) value
      end if
      np = np + 1
    end do

    return

100 write(error_unit,"('[',a,']: ',a)") trim(loc)                              &
                    , "Wrong format of the parameter '" // trim(name) //       &
                      "' or the value is too small or too large!"

!-------------------------------------------------------------------------------
!
  end subroutine get_parameter_integer
!
!===============================================================================
!
! subroutine GET_PARAMETER_REAL:
! -----------------------------
!
!   Subroutine reads a given parameter name and returns its real value.
!
!   Arguments:
!
!     name  - the input parameter name;
!     value - the output real value of parameter;
!
!===============================================================================
!
  subroutine get_parameter_real(name, value)

! import external procedures
!
    use iso_fortran_env, only : error_unit

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(in)    :: name
    real(kind=8)    , intent(inout) :: value

! local parameter counter
!
    integer :: np

! local parameters
!
    character(len=*), parameter :: loc = 'PARAMETERS::get_parameter_real()'
!
!-------------------------------------------------------------------------------
!
! find the selected parameter
!
    np = 1
    do while (np <= nparams)
      if (name == pnames(np)) then
        read(pvalues(np), err = 100, fmt = *) value
      end if
      np = np + 1
    end do

    return

100 write(error_unit,"('[',a,']: ',a)") trim(loc)                              &
                    , "Wrong format of the parameter '" // trim(name) //       &
                      "' or the value is too small or too large!"

!-------------------------------------------------------------------------------
!
  end subroutine get_parameter_real
!
!===============================================================================
!
! subroutine GET_PARAMETER_STRING:
! -------------------------------
!
!   Subroutine reads a given parameter name and returns its string value.
!
!   Arguments:
!
!     name  - the input parameter name;
!     value - the output string value of parameter;
!
!===============================================================================
!
  subroutine get_parameter_string(name, value)

! local variables are not implicit by default
!
    implicit none

! subroutine arguments
!
    character(len=*), intent(in)    :: name
    character(len=*), intent(inout) :: value

! local parameter counters
!
    integer :: np, nl
!
!-------------------------------------------------------------------------------
!
! get the length of the output string
!
    nl = min(vlen, len(value))

! find the selected parameter
!
    np = 1
    do while (np <= nparams)
      if (name == pnames(np)) then
        value = pvalues(np)(1:nl)
      end if
      np = np + 1
    end do

!-------------------------------------------------------------------------------
!
  end subroutine get_parameter_string

!===============================================================================
!
end module parameters
