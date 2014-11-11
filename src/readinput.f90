!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all right reserved
!
! Version: :VERSION:
! Date:    :DATE:
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build:   :BUILD:
! Author:  :AUTHOR:
! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This file is part of felixsim.
!
!  felixsim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  felixsim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with felixsim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: ReadInputFile.f90,v 1.23 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Calls the varoius functions which read in all the required data/parameters
!(doxygen)
SUBROUTINE ReadInput(IErr)

!Module call
  USE MyNumbers

  USE MPI
  USE MyMPI

IMPLICIT NONE

INTEGER(IKIND):: IErr

!Calling all functions which read felix.inp, felix.sca and felix.cif
!(felix.hkl too depending on user preference)
!ensure all input files are in working directory

  !felix.inp
  CALL ReadInpFile(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadInput(", my_rank, ") error",IErr, &
          " in ReadInpFile()"
     RETURN
  ENDIF

  !felix.hkl
  CALL ReadHklFile(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadInput(", my_rank, ") error",IErr, &
          "in ReadHklFile()"
     RETURN
  ENDIF

  !felix.sca
  CALL ReadScaFile(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadInput(", my_rank, ") error",IErr, &
          " in ReadScaFile()"
     RETURN
  ENDIF

  !felix.cif
  CALL ReadCifFile(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadInput(", my_rank, ") error",IErr, &
          "in ReadCifFile()"
     RETURN
  ENDIF

END SUBROUTINE ReadInput
