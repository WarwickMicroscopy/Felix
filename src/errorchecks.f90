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

SUBROUTINE ErrorChecks(ProgramInCurrently, ProgramWithIssue,IDamage,IErr)
  
  USE MyNumbers
  
  USE MPI
  USE MyMPI
  USE IPara
  
  USE IConst
  
  IMPLICIT NONE
  
  CHARACTER(*) ProgramInCurrently, ProgramWithIssue
  
  INTEGER(IKIND) IErr, IDamage
  
  !Only checks if IErr is anything other than zero,
  !Will be expanded to include individual error codes
  !Eventually case structure built up so for example
  !will say - "in ALLOCATION etc ..." for a certain IErr number
  IF (IErr.NE.ZERO) THEN
     
     SELECT CASE (IDamage)
     CASE(IWarning)

        SELECT CASE (IErr)
!!$           CASE(654)

!!$              IF (my_rank.EQ.0) THEN
!!$                 PRINT*,"DEBUG MESSAGE: ",ProgramInCurrently,"(", my_rank, ") Warning:", IErr,&
!!$                      "in ", ProgramWithIssue, " error in optional variable passing"
!!$                 PRINT*, "Check type of variable name matches the variable type in Call statement "
!!$                 PRINT*, "for example: IMinWeakBeams has to have Integer type variable - the compiler should pick this up "
!!$                 PRINT*, "if you are here it is likely the I is missing off of the IMinWeakBeams MessageVariable "
!!$              END IF
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
