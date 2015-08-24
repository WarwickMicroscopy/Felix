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

SUBROUTINE Message(ProgramName,IPriorityFLAG,IErr,MessageVariable,RVariable,IVariable,RVector,CVariable,MessageString)

  USE MyNumbers
  USE IPara

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
   
  CHARACTER(*), INTENT (IN), OPTIONAL ::  &
       MessageVariable, MessageString

  REAL(RKIND),INTENT(IN), OPTIONAL :: &
       RVariable, RVector(:) 
  INTEGER(IKIND),INTENT (IN), OPTIONAL :: &
       IVariable
  COMPLEX(CKIND),INTENT (IN), OPTIONAL :: &
       CVariable

  INTEGER(IKIND) :: &
       ISizeVector

  CHARACTER(*),INTENT (IN) :: &
       ProgramName
  INTEGER(IKIND) :: &
       IErr,IPriorityFLAG,ind
  
  CHARACTER*100 SVariable, SVariableTemp, my_rank_string,DebugString
  CHARACTER*30 SVariableString3DVector(THREEDIM)
  CHARACTER*30,DIMENSION(:), ALLOCATABLE :: &
       SVariableVector

!!$  Converts my_rank to string and then either the RVariable, IVariable, CVariable to a string
     WRITE(my_rank_string,'(I6.1)') my_rank 
     IF (PRESENT(RVariable).AND.IPriorityFLAG.LT.100) THEN
           WRITE(SVariable,'(F15.3)') RVariable
     ELSE IF (PRESENT(IVariable)) THEN
        WRITE(SVariable,'(I10.1)') IVariable
     ELSE IF (PRESENT(CVariable)) THEN
        WRITE(SVariable,'(F30.16)') CVariable

!!$  Converts A vector input into a string - finds the size of the input vector,
!!$  and allocates appropriately
     ELSE IF (PRESENT(RVector)) THEN
        ILenVector=SIZE(RVector)
        ALLOCATE(&
             SVariableVector(ILenVector),&
             STAT=IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"felixrefine (", my_rank, ") error in Allocation() of IIterativeVariableUniqueIDs"
           RETURN
        ENDIF
        
        DO ind=1,ILenVector
           WRITE(SVariableStringVector(ind),'F15.3)') RVector(ind)
        END DO

        SVariable=""
        SVariableTemp=""
        DO ind=1,ILenVector
           SVariableTemp=//"  "//TRIM(ADJUSTL(SVariableStringVector(ind)))//"  "
           SVariable=
        !DO ind=1,3
        !   WRITE(SVariableStringVector(ind),'(F15.3)') RVector(ind)
        !   !SVariable=SVariable //"  "//SVariableStringVector(ind)
        !END DO
        !SVariable="  "//TRIM(ADJUSTL(SVariableStringVector(1)))//"  "// &
        !     TRIM(ADJUSTL(SVariableStringVector(2)))//"  "//TRIM(ADJUSTL(SVariableStringVector(3)))
     ELSE
        SVariable = ""
     END IF

!!$  If IWriteFLAG is set to over 100 - IDebugFLAG is activated, IWriteFLAG set back to normal setting

  IF (IWriteFLAG.GE.100) THEN
     IDebugFLAG = IWriteFLAG
     IWriteFLAG = IDebugFLAG - 100
  END IF

  
!!$  If IPriorityFLAG is over 100 (Debug messaging) below won't execute
!!$  Prints out specified variation of message (dependent on presence of variables), to the screen
  IF (IPriorityFLAG .LT. 100) THEN 
 
     ! Checks if MessageVariable & MessageString has been read into the function
     IF (PRESENT(MessageVariable).AND.PRESENT(MessageString)) THEN
        !Prints out message
        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
             .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. IRefineSwitch) THEN
           PRINT*,ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ", &
                TRIM(MessageVariable)," = ",TRIM(ADJUSTL(SVariable)),"  ",TRIM(ADJUSTL(MessageString))
        END IF
        
     ELSE IF (PRESENT(MessageVariable)) THEN
       
        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
             .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. IRefineSwitch) THEN
           PRINT*,ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ",&
                TRIM(ADJUSTL(MessageVariable))," = ",TRIM(ADJUSTL(SVariable))
        END IF
        
     ELSE IF (PRESENT(MessageString)) THEN

        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
             .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. IRefineSwitch) THEN
           PRINT*,ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ", TRIM(ADJUSTL(MessageString))
        END IF
        
     ELSE
        IF((IPriorityFLAG.LE.IWriteFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
             .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. IRefineSwitch) THEN
           PRINT*,ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) "
        END IF
     END IF

!!$-----------------------------------------------------------------------------
!!$  below only executes if message is a debug message and set in debug mode
!!$  Debug messages are printed out here
  ELSE IF(IPriorityFLAG .GE. 100) THEN

!!$     Prints out reals with greater precision
     IF(PRESENT(RVariable)) THEN
        WRITE(SVariable,'(F30.16)') RVariable
     END IF

     ! Checks if MessageVariable & MessageString has been read into the function
     IF (PRESENT(MessageVariable).AND.PRESENT(MessageString)) THEN
        !Prints out message
        IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
             .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. IRefineSwitch) THEN
           PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ", &
                TRIM(ADJUSTL(MessageVariable))," = ",TRIM(ADJUSTL(SVariable)),"  ",TRIM(ADJUSTL(MessageString))
        END IF
        
     ELSE IF (PRESENT(MessageVariable)) THEN
       
        IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
             .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. IRefineSwitch) THEN
           PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ",TRIM(ADJUSTL(MessageVariable)),&
                " = ",TRIM(ADJUSTL(SVariable))
        END IF

     ELSE IF (PRESENT(MessageString)) THEN

        IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
             .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. IRefineSwitch) THEN
           PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ", TRIM(ADJUSTL(MessageString))
        END IF
        
        
     ELSE

        IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
             .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. IRefineSwitch) THEN
           PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) "
        END IF
     END IF
  END IF


END SUBROUTINE Message
   
END MODULE WriteToScreen
  
