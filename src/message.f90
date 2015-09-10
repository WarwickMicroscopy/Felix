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

  SUBROUTINE Message(ProgramName,IPriorityFLAG,IErr,MessageVariable,RVariable,IVariable, &
       RVector,RMatrix,CVariable,MessageString)

    USE MyNumbers
    USE IPara

    USE MPI
    USE MyMPI

    IMPLICIT NONE

    CHARACTER(*), INTENT (IN), OPTIONAL ::  &
         MessageVariable, MessageString
    REAL(RKIND),INTENT(IN), OPTIONAL :: &
         RVariable, RVector(:), RMatrix(:,:) 
    INTEGER(IKIND),INTENT (IN), OPTIONAL :: &
         IVariable
    COMPLEX(CKIND),INTENT (IN), OPTIONAL :: &
         CVariable

    INTEGER(IKIND) :: &
         ISizeVector,ISizeMatrixX,ISizeMatrixY

    CHARACTER(*),INTENT (IN) :: &
         ProgramName
    INTEGER(IKIND) :: &
         IErr,IPriorityFLAG,ind,jnd

    CHARACTER*100 &
         SVariable, SVariableTemp,SVariableOld, my_rank_string,DebugString,Sind
    CHARACTER*30,DIMENSION(:), ALLOCATABLE :: &
         SVariableVector
    CHARACTER*30,DIMENSION(:,:), ALLOCATABLE :: &
         SVariableMatrix

    INTEGER(IKIND) :: &
         IMatrixPresentSwitch,ILengthofLine,IMaxLengthIndicator,ILineBreaks

    !Variable that switches message subroutine to matrix printing mode 
    IMatrixPresentSwitch=0

!!$  Converts my_rank to string and then either the RVariable, IVariable, CVariable to a string
    WRITE(my_rank_string,'(I6.1)') my_rank 
    IF (PRESENT(RVariable).AND.IPriorityFLAG.LT.100) THEN
       WRITE(SVariable,'(F15.3)') RVariable
    ELSE IF (PRESENT(IVariable)) THEN
       WRITE(SVariable,'(I10.1)') IVariable
    ELSE IF (PRESENT(CVariable)) THEN
       WRITE(SVariable,'(F30.16)') CVariable

!!$  Converts a vector input into a string - finds the size of the input vector,
!!$  and allocates appropriately. Stores in a long string SVariable
    ELSE IF (PRESENT(RVector)) THEN
       ISizeVector=SIZE(RVector)
       ALLOCATE(&
            SVariableVector(ISizeVector),&
            STAT=IErr)
       IF( IErr.NE.0 ) THEN
          PRINT*,"felixrefine (", my_rank, ") error in Allocation() of IIterativeVariableUniqueIDs"
          RETURN
       ENDIF

       DO ind=1,ISizeVector
          WRITE(SVariableVector(ind),'(F15.3)') RVector(ind)
       END DO

       SVariableOld=""
       SVariableTemp=""
       DO ind=1,ISizeVector
          SVariableTemp="  "//TRIM(ADJUSTL(SVariableVector(ind)))//"  "
          SVariable=TRIM(ADJUSTL(SVariableOld))//" "// TRIM(ADJUSTL(SVariableTemp))//"  "
          SVariableOld=SVariable
          SVariableTemp=""
       END DO

    ELSE IF (PRESENT(RMatrix)) THEN
       IMatrixPresentSwitch=1
       ISizeMatrixX=SIZE(RMatrix,1,IKIND)
       ISizeMatrixY=SIZE(RMatrix,2,IKIND)
       ALLOCATE(&
            SVariableMatrix(ISizeMatrixX,ISizematrixY), &
            STAT=IErr)
       IF( IErr.NE.0 ) THEN
          PRINT*,"felixrefine (", my_rank, ") error in Allocation() of IIterativeVariableUniqueIDs"
          RETURN
       ENDIF

       DO ind=1,ISizeMatrixX
          DO jnd=1,ISizeMatrixY
             WRITE(SVariableMatrix(ind,jnd),'(F15.3)') RMatrix(ind,jnd)
          END DO
       END DO

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
    IF (IPriorityFLAG .LT. 100.AND.IMatrixPresentSwitch.EQ.0) THEN 

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
    ELSE IF(IPriorityFLAG .GE. 100.AND.IMatrixPresentSwitch.EQ.0) THEN

!!$     Prints out Reals with greater precision
       IF(PRESENT(RVariable)) THEN
          WRITE(SVariable,'(F30.16)') RVariable
       END IF

       ! Checks if MessageVariable & MessageString has been read into the function
       IF (PRESENT(MessageVariable).AND.PRESENT(MessageString)) THEN
          !Prints out message
          IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
               .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. IRefineSwitch) THEN
             PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ", &
                  TRIM(ADJUSTL(MessageVariable))," = ",TRIM(ADJUSTL(SVariable)),"  ", &
                  TRIM(ADJUSTL(MessageString))
          END IF

       ELSE IF (PRESENT(MessageVariable)) THEN

          IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
               .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. IRefineSwitch) THEN
             PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ", &
                  TRIM(ADJUSTL(MessageVariable))," = ",TRIM(ADJUSTL(SVariable))
          END IF

       ELSE IF (PRESENT(MessageString)) THEN

          IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
               .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. IRefineSwitch) THEN
             PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ", & 
                  TRIM(ADJUSTL(MessageString))
          END IF

       ELSE

          IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
               .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. IRefineSwitch) THEN
             PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) "
          END IF
       END IF

!!$-----------------------------------------------------------------------------
!!$  Below executes in the special case of printing out a matrix - *only for DebugMode*
!!$  We need to loop through each row of the matrix and Print to the screen in felix 
!!$  message format.
    ELSE IF(IPriorityFLAG.GE.100.AND.IMatrixPresentSwitch.EQ.1) THEN

       IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
            .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. IRefineSwitch) THEN

!!$       Loop over rows (ISizeMatrixX) and columns (ISizeMatrixY) - concatenate each row into 
!!$       a dummy variable: SVariable, if matrix is too long (row wise) counter counts how many 
!!$       Line Breaks are required
          DO ind=1,ISizeMatrixX
             SVariableOld=""
             SVariableTemp=""
             ILineBreaks=1
             DO jnd =1,ISizeMatrixY
!!$             Could be an input defined variable
!!$             ************************
                ILengthofLine=10
!!$             ************************
                IMaxLengthIndicator=MOD(jnd,ILengthofline)
                IF(IMaxLengthIndicator.EQ.0) THEN
                   ILineBreaks=ILineBreaks+1
                END IF
                SVariableTemp="  "//TRIM(ADJUSTL(SVariableMatrix(ind,jnd)))//"  "
                SVariable=TRIM(ADJUSTL(SVariableOld))//" "// TRIM(ADJUSTL(SVariableTemp))//"  "
                SVariableOld=SVariable
                SVariableTemp=""
             END DO

!!$          We need to decipher the number of Characters in one line, this will determine the
!!$          Upper limit of the substrings
             INumCharactersinFinalLine=MOD(LEN(SVariable),ILineBreaks-1)
             INumFullLineCharacters=LEN(SVariable)-INumCharactersinFinalLine
             INumCharactersinOneLine=INumFullLineCharacters-(ILineBreaks-1)

             WRITE(Sind,'(I8.1)') ind

!!$             Checks if MessageVariable & MessageString has been read into the function
!!$             Three Conditions - row is one line long or less, row is 2 lines long or row is
!!$             greater than three lines long. The Print Statement is added for each case, this
!!$             provides a consistent matrix format no matter the length of the row
             IF (PRESENT(MessageVariable).AND.PRESENT(MessageString)) THEN
                IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
                     .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. IRefineSwitch) THEN
                   IF(ISizeMatrixY.LE.ILengthofLine) THEN
                      PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ", &
                           TRIM(ADJUSTL(MessageVariable)),"(",TRIM(ADJUSTL(Sind)),":) = ", &
                           TRIM(ADJUSTL(SVariable)),"  ",TRIM(ADJUSTL(MessageString))
                   ELSE
                      PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ", &
                           TRIM(ADJUSTL(MessageVariable)),"(",TRIM(ADJUSTL(Sind)),":) = ", &
                           TRIM(ADJUSTL(SVariable(1:INumCharactersinOneLine))),"  ",TRIM(ADJUSTL(MessageString))
                      IF(ILineBreaks.GT.2) THEN
                         DO knd =1,ILineBreaks-1
                            PRINT*,"DBG_MESSAGE:                            ", &
                                 TRIM(ADJUSTL(SVariable((knd*INumCharactersinOneLine)+1:(knd+1)*INumCharactersinOneLine))), &
                                 "  ",TRIM(ADJUSTL(MessageString))
                         END DO
                      END IF
                      PRINT*,"DBG_MESSAGE:                            ", &
                           TRIM(ADJUSTL(SVariable((ILineBreaks*INumCharactersinOneLine)+1:LEN(SVariable)))), &
                           "  ",TRIM(ADJUSTL(MessageString))
                   END IF
                END IF

             ELSE IF (PRESENT(MessageVariable)) THEN

                IF((IPriorityFLAG.LE.IDebugFLAG.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.IRefineSwitch) &
                     .OR.IDebugFLAG.GE.110.AND.ISoftwareMode .LT. IRefineSwitch) THEN
                   IF(ISizeMatrixY.LE.ILengthofLine) THEN
                      PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ", &
                           TRIM(ADJUSTL(MessageVariable)),"(",TRIM(ADJUSTL(Sind)),":) = ", &
                           TRIM(ADJUSTL(SVariable)),"  ",TRIM(ADJUSTL(MessageString))
                   ELSE
                      PRINT*,"DBG_MESSAGE: ",ProgramName,"( ",TRIM(ADJUSTL(my_rank_string))," ) ", &
                           TRIM(ADJUSTL(MessageVariable)),"(",TRIM(ADJUSTL(Sind)),":) = ", &
                           TRIM(ADJUSTL(SVariable(1:INumCharactersinOneLine))),"  ",TRIM(ADJUSTL(MessageString))
                      IF(ILineBreaks.GT.2) THEN
                         DO knd =1,ILineBreaks-1
                            PRINT*,"DBG_MESSAGE:                            ", &
                                 TRIM(ADJUSTL(SVariable((knd*INumCharactersinOneLine)+1:(knd+1)*INumCharactersinOneLine))), &
                                 "  ",TRIM(ADJUSTL(MessageString))
                         END DO
                      END IF
                      PRINT*,"DBG_MESSAGE:                            ", &
                           TRIM(ADJUSTL(SVariable((ILineBreaks*INumCharactersinOneLine)+1:LEN(SVariable)))), &
                           "  ",TRIM(ADJUSTL(MessageString))
                   END IF
                END IF
             END IF
          END DO
       END IF
    END IF

  END SUBROUTINE Message

END MODULE WriteToScreen

