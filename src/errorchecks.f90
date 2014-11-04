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
! $Id: FelixSim.f90,v 1.89 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE ErrorChecks(ProgramInCurrently, ProgramWithIssue, IErr)

    USE MyNumbers

    USE MPI
    USE MyMPI

    IMPLICIT NONE

    CHARACTER*40 ProgramInCurrently, ProgramWithIssue

    INTEGER(IKIND) IErr

    !Only checks if IErr is anything other than zero,
    !Will be expanded to include individual error codes
    !Eventually case structure built up so for example
    !will say - "in ALLOCATION etc ..." for a certain IErr number

    IF (IErr.NE.ZERO) THEN
       PRINT*,ProgramInCurrently,"(", my_rank, ") error", IErr,&
            "in ", ProgramWithIssue

       !Closes MPI Interface 
       CALL MPI_Finalize(IErr)

       IF( IErr.NE.ZERO ) THEN
          PRINT*,"ErrorChecks(", my_rank, ") error ", IErr, " in MPI_Finalize()"
          STOP
       ENDIF
       
       !Clean Shutdown
       STOP
    ENDIF

END SUBROUTINE ErrorChecks

