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
SUBROUTINE Message(MessageString,MessageVariable,IMessageFLAG,IPriorityFLAG,IErr)

  USE MyNumbers
  USE IPara

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  
 
  CHARACTER*52  MessageString

  !not sure how long the longest variable name is
  !could well need to be changed
  CHARACTER*32 :: MessageVariable,LeftMessage,RightMessage

  INTEGER(IKIND) :: IMessageFLAG,IErr,ILeftParenthesisPosition,IPriorityFLAG

!  INTEGER(IKIND) MessageLen
  !PRINT*,MessageString,LEN(MessageString)

  !find position of left Parenthesis - then Split strings
  ILeftParenthesisPosition = SCAN(MessageString,"(") 
  LeftMessage = MessageString(1:ILeftParenthesisPosition)
  RightMessage = MessageString(ILeftParenthesisPosition + 1 : ILeftParenthesisPosition + 1)


  SELECT CASE (IMessageFlag)
     
  !Prints Variable name as well as Message
  CASE(0)
     
     IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
        PRINT*,TRIM(LeftMessage),my_rank,TRIM(RightMessage),MessageVariable
     END IF

  !Prints Message without variable
  CASE(1)
     IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN
        PRINT*,TRIM(LeftMessage),my_rank,TRIM(RightMessage)
     END IF
     
        
  CASE DEFAULT
     !Error message call here from subroutine errorchecks.f90,
     !with associated error flag
  END SELECT
  
END SUBROUTINE Message
  
  
