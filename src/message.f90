!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all rights reserved
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


!Subroutine reads in a message the user wants to display to the screen, The IWriteFLAG
!is then compared to the PriorityFLAG which determines whether it will get printed out or not.
MODULE WriteToScreen
CONTAINS

SUBROUTINE Message(ProgramName,IPriorityFLAG,IErr,MessageVariable,RVariable,IVariable,CVariable,MessageString)

  USE MyNumbers
  USE IPara

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
   
  CHARACTER(*), INTENT (IN), OPTIONAL ::  MessageVariable, MessageString

  REAL(RKIND),INTENT(IN), OPTIONAL :: RVariable 
  INTEGER(IKIND),INTENT (IN), OPTIONAL ::IVariable
  COMPLEX(CKIND),INTENT (IN), OPTIONAL :: CVariable

  CHARACTER(*),INTENT (IN) :: ProgramName
  INTEGER(IKIND) :: IErr,IPriorityFLAG
  
  CHARACTER*30 VariableString, my_rank_string

  IF (my_rank == 0) THEN
     WRITE(my_rank_string,'(I1)') my_rank 
     IF (PRESENT(RVariable)) THEN
        WRITE(VariableString,'(F30.16)') RVariable
     ELSE IF (PRESENT(IVariable)) THEN
        WRITE(VariableString,'(I10.1)') IVariable
     ELSE IF (PRESENT(CVariable)) THEN
        WRITE(VariableString,'(F30.16)') CVariable  
     ELSE
        VariableString = ""
     END IF
  END IF

  IF (IWriteFLAG.GE.100) THEN
     IDebugFLAG = IWriteFLAG
     IWriteFLAG = IDebugFLAG - 100
  END IF
!!$  If IPriorityFLAG is over 100 (Debug messaging) below won't execute
 
  IF (IPriorityFLAG .LT. 100) THEN 
 
     ! Checks if MessageVariable has been read into the function
     IF (PRESENT(MessageVariable).AND.PRESENT(MessageString)) THEN
        
        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
             .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
           PRINT*,ProgramName,"( ",TRIM(my_rank_string)," ) ", &
                TRIM(MessageVariable)," = ",TRIM(VariableString),"  ",TRIM(MessageString)
        END IF
        
     ELSE IF (PRESENT(MessageVariable)) THEN
       
        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
             .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
           PRINT*,ProgramName,"( ",TRIM(my_rank_string)," ) ",TRIM(MessageVariable)," = ",TRIM(VariableString)
        END IF

     ELSE IF (PRESENT(MessageString)) THEN
        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
             .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
           PRINT*,ProgramName,"( ",TRIM(my_rank_string)," ) ", TRIM(MessageString)
        END IF
        
        
     ELSE
        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
             .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
           PRINT*,ProgramName,"( ",TRIM(my_rank_string)," ) "
        END IF
     END IF

!END IF
        
!!$        !If it has, we then check to see which type of variable needs to be printed
!!$        !out to the screen
!!$        IF (PRESENT(RVariable).AND.PRESENT(MessageString)) THEN
!!$           
!!$           !Check for mismatch
!!$           IF (SCAN(MessageVariable,'R').EQ.0.AND.IWriteFLAG.GE.100) THEN
!!$              !654 is error code associated with the misinterpreted variable
!!$              CALL ErrorChecks("Message","Message",654)
!!$           ENDIF
!!$
!!$           IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
!!$                .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
!!$              PRINT*,ProgramName,"(",my_rank,") ",MessageVariable," =",RVariable
!!$           END IF
!!$           
!!$        ELSE IF (PRESENT(IVariable).AND.PRESENT(MessageString)) THEN
!!$           
!!$           !Check for mismatch
!!$           IF (SCAN(MessageVariable,'I').EQ.0.AND.IWriteFLAG.GE.100) THEN
!!$              !654 is error code associated with the misinterpreted variable
!!$              CALL ErrorChecks("Message","Message",654)
!!$           ENDIF
!!$           
!!$           IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
!!$                .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
!!$              PRINT*,ProgramName,"(",my_rank,") ",MessageVariable," =",IVariable
!!$           END IF
!!$           
!!$        ELSE IF (PRESENT(CVariable)) THEN
!!$           IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
!!$                .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
!!$              PRINT*,ProgramName,"(",my_rank,") ",MessageVariable," =",CVariable
!!$           END IF
!!$           
!!$           !For character variable
!!$        ELSE IF (PRESENT(SVariable).AND.PRESENT(MessageVariable)) THEN
!!$
!!$           IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
!!$                .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
!!$              PRINT*,ProgramName,"(",my_rank,") ",MessageVariable, SVariable,MessageString
!!$           END IF
!!$           
!!$        ELSE 
!!$           IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
!!$                .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
!!$              PRINT*,ProgramName,"(",my_rank,") ",MessageVariable,SVariable,RVariable,IVariable
!!$           END IF
!!$        
!!$        END IF
           
!!$           !If there is only a message to print out, below is executed
!!$     ELSE IF (PRESENT(MessageString)) THEN
!!$        
!!$        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
!!$             .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
!!$           PRINT*,ProgramName,"(",my_rank,") ",MessageString
!!$        END IF
!!$     
!!$        !if there is no message or variable - the function is just printed out
!!$        !to identify it is entering that program
!!$     ELSE
!!$        
!!$        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
!!$             .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
!!$           PRINT*,ProgramName,"(",my_rank,") "
!!$        END IF
!!$     END IF

!!$     !below only executes if message is debug message and set as debug mode
  
  ELSE IF(IDebugFLAG .GE. 100) THEN

 ! Checks if MessageVariable has been read into the function
     IF (PRESENT(MessageVariable).AND.PRESENT(MessageString)) THEN
        
        IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
             .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. 2) THEN
           PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(my_rank_string)," ) ", &
                TRIM(MessageVariable)," = ",TRIM(VariableString),"  ",TRIM(MessageString)
        END IF
        
     ELSE IF (PRESENT(MessageVariable)) THEN
       
        IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
             .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. 2) THEN
           PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(my_rank_string)," ) ",TRIM(MessageVariable)," = ",TRIM(VariableString)
        END IF

     ELSE IF (PRESENT(MessageString)) THEN
        IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
             .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. 2) THEN
           PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(my_rank_string)," ) ", TRIM(MessageString)
        END IF
        
        
     ELSE
        IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
             .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. 2) THEN
           PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(my_rank_string)," ) "
        END IF
     END IF
  END IF





     
!!$     ! Checks if MessageVariable has been read into the function
!!$     IF (PRESENT(MessageVariable)) THEN
!!$        
!!$        !If it has, we then check to see which type of variable needs to be printed
!!$        !out to the screen
!!$        IF (PRESENT(RVariable)) THEN
!!$           
!!$           !Check for mismatch
!!$           IF (SCAN(MessageVariable,'R').EQ.0) THEN
!!$              !654 is error code associated with the misinterpreted variable
!!$              CALL ErrorChecks("Message","Message",654)
!!$           ENDIF
!!$
!!$           IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) & 
!!$                .OR.IWriteFLAG.GE.110.AND.ISoftwareMode .LT. 2) THEN
!!$              PRINT*,ProgramName,"(",my_rank,") ",MessageVariable," =",RVariable
!!$           END IF
!!$           
!!$        ELSE IF (PRESENT(IVariable)) THEN
!!$
!!$           !Check for mismatch
!!$           IF (SCAN(MessageVariable,'I').EQ.0) THEN
!!$              !654 is error code associated with the misinterpreted variable
!!$              CALL ErrorChecks("Message","Message",654)
!!$           ENDIF
!!$
!!$           IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) & 
!!$                .OR.IWriteFLAG.GE.110.AND.ISoftwareMode .LT. 2) THEN
!!$              PRINT*,ProgramName,"(",my_rank,") ",MessageVariable," =",IVariable
!!$           END IF
!!$           
!!$        ELSE IF (PRESENT(CVariable)) THEN
!!$           IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
!!$                .OR.IWriteFLAG.GE.110.AND.ISoftwareMode .LT. 2) THEN
!!$              PRINT*,ProgramName,"(",my_rank,") ",MessageVariable," =",CVariable
!!$           END IF
!!$        END IF
!!$     
!!$     
!!$        !If There is a message with the FunctionName then the Message is Printed with its message
!!$     ELSE IF (PRESENT(MessageString)) THEN
!!$        
!!$        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
!!$             .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
!!$           PRINT*,ProgramName,"(",my_rank,")",MessageString
!!$        END IF
!!$    
!!$        !if there is no message or variable - the function is just printed out
!!$        !to identify it is entering that program
!!$     ELSE
!!$        
!!$        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) & 
!!$             .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
!!$           PRINT*,ProgramName,"(",my_rank,")"
!!$        END IF
!!$     END IF
!!$  END IF

END SUBROUTINE Message
   
END MODULE WriteToScreen
  
