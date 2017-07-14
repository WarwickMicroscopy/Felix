!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-17, all rights reserved
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
!  Felix is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  Felix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with Felix.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: FelixSim.f90,v 1.89 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ErrorChecks(ProgramInCurrently, ProgramWithIssue,IDamage,IErr)
  
  USE MyNumbers
  USE ErrorCodes
  
  USE MPI
  USE MyMPI
  USE IPara
  
  USE IConst
  
  IMPLICIT NONE
  
  CHARACTER(*) ProgramInCurrently, ProgramWithIssue
  CHARACTER*100 :: myrankstring, IErrString
  
  INTEGER(IKIND) IErr, IDamage

!!$  write my rank and Ierr into strings
  WRITE(myrankstring,*) my_rank
  WRITE(IErrstring,*) IErr
  
!!$  Only checks if IErr is anything other than zero,
!!$  Will be expanded to include individual error codes
!!$  Eventually case structure built up so for example
!!$  will say - "in ALLOCATION etc ..." for a certain IErr number
  IF (IErr.NE.ZERO) THEN
     
     SELECT CASE (IDamage)
     CASE(IWarning)

        SELECT CASE (IErr)

        CASE DEFAULT

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
        END SELECT
        
     CASE(IPotError)
        
        SELECT CASE (IErr)
        CASE(IReflectionMismatch)
           IF (my_rank.EQ.0) THEN
              PRINT*,TRIM(ADJUSTL(ProgramInCurrently)) //"(" // TRIM(ADJUSTL(myrankString)) // &
                   ") Potential Error " // TRIM(ADJUSTL(IErrString)) // " in "// TRIM(ADJUSTL(ProgramWithIssue))
              PRINT*, "Number of Reflections quantised in each Laue Zone does not match the total Reflections in the system"
              PRINT*, "Please Contact: Keith.Evans@Warwick.ac.uk or a.j.m.hubert@warwick.ac.uk for help,"
              PRINT*, "Problem is very likely to be a source code bug, we want to know about those!"
              PRINT*, "To continue wih the simulation, please switch the RAcceptanceAngle in the input file to zero,"
              PRINT*, "and run the simulation again"
           
              IErr=IReflectionMismatch
           END IF
           
        CASE DEFAULT
           
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
        END SELECT
        
     CASE(ICritError)
        
        SELECT CASE (IErr)             
        CASE DEFAULT

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
        END SELECT
     END SELECT
  END IF
  
  
END SUBROUTINE ErrorChecks
