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
!is then compared to the priorityFLAG which determines whether it will get printed out or not.
MODULE ScreenWrite
CONTAINS
SUBROUTINE Message(ProgramName,IPriorityFLAG,IErr,MessageVariableIn,RVariableIn,IVariableIn,CVariableIn,MessageStringIn)

  USE MyNumbers
  USE IPara

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
   
  CHARACTER(*), INTENT (IN), OPTIONAL ::  MessageVariableIn, MessageStringIn 

  REAL(RKIND),INTENT(IN), OPTIONAL :: RVariableIn 
  INTEGER(IKIND),INTENT (IN), OPTIONAL ::IVariableIn
  COMPLEX(CKIND),INTENT (IN), OPTIONAL :: CVariableIn

  
  REAL(RKIND):: RVariable 
  INTEGER(IKIND) :: IVariable
  COMPLEX(CKIND) :: CVariable

  CHARACTER(*),INTENT (IN) :: ProgramName
  CHARACTER(LEN=LEN(MessageVariableIn)) :: MessageVariable
  CHARACTER(LEN=LEN(MessageStringIn)) :: MessageString
  INTEGER(IKIND) :: IErr,IPriorityFLAG


  LOGICAL :: R,U

  PRINT*," IPriorityFLAG   ",IPriorityFLAG
  PRINT*,"IErr &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&7", IErr

   PRINT*,"ProgramName = ", ProgramName
  ! PRINT*,MessageVariable
  
   ! R = PRESENT(MessageVariableIn)
  ! PRINT*,"R present?", R

   U = PRESENT(RVariableIn)
   PRINT*,"U present?", U

 ! PRINT*,"MessageVariable = ", MessageVariable

  ! Checks if MessageVariable has been read into the function
  IF (PRESENT(MessageVariableIn)) THEN
     MessageVariable = MessageVariableIn
    
     !If it has, we then check to see which type of variable needs to be printed
     !out to the screen
     IF (PRESENT(RVariableIn)) THEN
        RVariable = RVariableIn
        !Check for mismatch
       ! IF (SCAN(MessageVariable,'R').EQ.0.AND.IDebugMODE.EQ.1) THEN
           !654 is error code associated with the misinterpreted variable
        !   CALL ErrorChecks("Message","Message",654)
        
        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
           PRINT*,ProgramName,"(",my_rank,") ",MessageVariable," =",RVariable
        END IF
           
     ELSE IF (PRESENT(IVariableIn)) THEN

        IVariable = IVariableIn

        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
           PRINT*,ProgramName,"(",my_rank,") ",MessageVariable," =",IVariable
        END IF
           PRINT*," here"
     ELSE IF (PRESENT(CVariableIn)) THEN
        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
           PRINT*,ProgramName,"(",my_rank,") ",MessageVariable," =",CVariable
        END IF
     END IF
     
     
     !If There is a message with the FunctionName then the Message is Printed with its message
  ELSE IF (PRESENT(MessageStringIn)) THEN
     MessageString = MessageStringIn
     IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
        PRINT*,ProgramName,"(",my_rank,")",MessageString
     END IF
    
     !if there is no message or variable - the function is just printed out
     !to identify it is entering that program
  ELSE
     PRINT*," here.........................................................."
     IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
        PRINT*,ProgramName,"(",my_rank,")"
     END IF
  END IF

END SUBROUTINE Message
   
END MODULE
     


  


!old way
 ! INTEGER(IKIND) MessageLen
  !PRINT*,MessageString,LEN(MessageString)

 !find position of left Parenthesis of MessageString - then Split strings to input my_rank
 ! ILeftParenthesisPosition = SCAN(MessageString,"(") 
 ! LeftMessage = MessageString(1:ILeftParenthesisPosition)
 ! RightMessage = MessageString(ILeftParenthesisPosition + 1 : ILeftParenthesisPosition + 1)
  
  !PRINT*,MessageVariable

 ! IF (IMessageFLAG.GE.1) THEN
   !  !find '=' sign of Message Variable - cut at this point to remove random memory rubbish
   !  ILeftParenthesisPosition = SCAN(MessageVariable,"=")
   !  MessageVariableShort = MessageVariable(1:ILeftParenthesisPosition)
 ! END IF

!   ! SELECT CASE (IMessageFlag)
     
!   !Prints Just Variable name 
!   CASE(0)

!      IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
!         PRINT*,MessageString,"(",my_rank,")"
!      END IF
     
!   !Prints Message with variable
!   CASE(1)
    
!      IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
!         PRINT*,MessageString,"(",my_rank,") ",MessageVariable," =",RVariable
!      END IF
        
!   CASE DEFAULT
!      !Error message call here from subroutine errorchecks.f90,
!      !with associated error flag
!   END SELECT
  
! END SUBROUTINE Message
  
