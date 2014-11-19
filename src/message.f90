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
SUBROUTINE Message(MessageString,MessageVariable,RVariable,IMessageFLAG,IPriorityFLAG,IErr)

  USE MyNumbers
  USE IPara

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  
 
  CHARACTER(*), INTENT (IN) :: MessageString, MessageVariable

  !not sure how long the longest variable name is
  !could well need to be changed
  !CHARACTER*32 :: LeftMessage,RightMessage,MessageVariableShort

  INTEGER(IKIND) :: IMessageFLAG,IErr,IPriorityFLAG,ILeftParenthesisPosition
  REAL(RKIND) :: RVariable
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

 ! PRINT*,"HERE?",MessageVariableShort

  SELECT CASE (IMessageFlag)
     
  !Prints Just Variable name 
  CASE(0)

     IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
        PRINT*,MessageString,"(",my_rank,")"
     END IF
     
  !Prints Message with variable
  CASE(1)
    
     IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
        PRINT*,MessageString,"(",my_rank,") ",MessageVariable," =",RVariable
     END IF
        
  CASE DEFAULT
     !Error message call here from subroutine errorchecks.f90,
     !with associated error flag
  END SELECT
  
END SUBROUTINE Message
  
  
