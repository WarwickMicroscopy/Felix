!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!  along with Felix.  IF not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
MODULE message_mod
!>--------------------------------------------------------------------------------------------
!>  major-authors: Jacob Richardson (2017)
!>--------------------------------------------------------------------------------------------
!>  MODULE OVERVIEW:
!>
!>    This MODULE directly CONTAINS everything for message,
!>      which is a major SUBROUTINE used throughout felix. It replaces the FORTRAN print
!>      statement and provides formatted output chosen at run-time.
!>     Although there is nothing wrong with the FORTRAN print statement
!>     and this module just makes everything more complicated and harder to debug.
!>     
!>
!>    message_mod also accesses l_alert_mod:
!>      l_alert is a simple but important function, used extensively for error handling.
!>
!>    message_mod also conatins print_end_time used to neatly print the time elapsed.
!>--------------------------------------------------------------------------------------------
!>  message( ... ) FEATURES:
!>
!>    IWriteFLAG in felix.inp is the run-time input used to select the output mode of Felix.
!>    message uses IWriteFLAG to decide whether a message should be printed.
!>
!>    - message prints various different variable TYPEs with consistent formatting
!>    - It accommodates for MPI parallel programming used in Felix by printing on one 'core'
!>    - It prints a tiered-structure to distinguish key messages from more niche ones
!>    - It speeds up debugging and testing, by being very quick and easy to use                   
!>--------------------------------------------------------------------------------------------
!>
!>  message( MsgPriority, MsgTag, text_to_print, variable_to_print )
!>
!>    MsgPriority         DERIVED TYPE ( msg_priority_TYPE ) = LS, LM, LL or LXL
!>    (optional)
!>                This determines how important it is to print this message. Is it
!>                stating key operations of Felix or more particular, finer details?                       
!>             
!>                LS  : high priority (small level of output) for important statements & values
!>                LM  : medium prioriity (medium level of output)
!>                LL  : low priority (large output) for fines details & thorough variable output
!>                LXL : extremely low priority large output) & niche details
!>
!>                DEFAULT = LS
!>
!>                At run-time this priority will be compared against the selected output mode
!>                given by IWriteFLAG in felix.inp. This will decide whether a message is
!>                is important enough to be printed. Additionally, Low priority  
!>                messages are indented further to the right to produce a tiered-output
!>                for clarity and readability.  
!>
!>    MsgTag              DERIVED TYPE ( msg_tag_TYPE )
!>    (optional)         
!>                This allows for messages to be grouped together with tags. At run-time you 
!>                can select that only key messages and messages with a certain tag are printed.
!>
!>    SMainMsg            TYPE: String (1D CHARACTER array)
!>
!>                This is some text to print to the terminal.
!>                Any variable values will be printed after this text. 
!>                        
!>    variable_to_print   TYPE: COMPLEX, INTEGER, REAL, LOGICAL, String 
!>    (optional)          DIMENSION: scalar, vector or matrix
!>        
!>                The value of this variable will be printed to the terminal with
!>                a format depending upon its TYPE and DIMENSION. E.g. a REAL matrix
!>                will be printed on multiple lines in columns with fixed width decimal notation.                
!>                       
!>--------------------------------------------------------------------------------------------
!>
!>  Using IWriteFLAG to choose the terminal output mode:
!>
!>    In felix.inp, IWriteFLAG can be assigned different INTEGER values to affect which 
!>    messages are printed. 
!>
!>    The digit of the INTEGER affects the minimum priority of the message printed:
!>      
!>        0 = LS,   2 = LM,   4 = LL,   8 = LXL   
!>
!>        sidenote: 8 rarely used, it prints every message and creates a very verbose output.
!>        LXL messages are usually accessed via a MsgTag instead.
!>
!>    The tens component of the INTEGER relates to the MsgTag. 
!>        
!>        0   - no_tag  : No indentation or tiered-structure and no messages have no_tag tag
!>                        so only priority is considered.
!>        10  - mode1   : This turns on tiered-structure and no messages have mode1 tag 
!>                        so only priority is considered.
!>        30  - dbg3    : Corresponds to IWriteFLAG = 3 IN THE OLD SYSTEM, many print
!>                        statements in the code previously only printed when IWriteFLAG = 3.
!>                        This will print all messages with dbg3 tag and priority messages
!>                        depending upon the digit.
!>
!>        For the full list of MsgTag's and their corresponding IWriteFLAG ten's value,
!>        refer to SetMessageMode() below.
!>
!>    IWriteFLAG examples: 
!>        2  will print all LS, LM messages
!>        32 will print all LS, LM, and dbg3 tagged messages
!>        96 will print all LS, LM, LL, and dbg14 tagged messages 
!>
!>--------------------------------------------------------------------------------------------

  USE l_alert_mod             ! grants access to felix's main error handling
  USE MyNumbers               ! necesary for IKIND, RKIND etc.
  USE MyMPI, ONLY : my_rank   ! necesary for my_rank - used to print on core 0 only

  !?? JR this grants MyNumbers, USE MyNumbers throughout code is for clarity not necessity 

  !-------------------------------------------------------------------------------------------
  ! HUGE INTERFACE to handle various variables TYPEs and optional arguments

  INTERFACE message
    MODULE PROCEDURE message1String 
    MODULE PROCEDURE messageRMatrix        
    MODULE PROCEDURE messageCMatrix       
    MODULE PROCEDURE messageIMatrix
    MODULE PROCEDURE message2Strings 
    MODULE PROCEDURE messageLogical      
    MODULE PROCEDURE message2Strings2Integers    ! prints text, INTEGER, more text, INTEGER 
    ! interfaces for scalars and vector variables, which then go to matrix messages
    MODULE PROCEDURE messageRVector        
    MODULE PROCEDURE messageCVector         
    MODULE PROCEDURE messageIVector         
    MODULE PROCEDURE messageRScalar            
    MODULE PROCEDURE messageCScalar
    MODULE PROCEDURE messageIScalar
    ! interfaces for optional arguments
    MODULE PROCEDURE message1String2            ! suffix :
    MODULE PROCEDURE message1String3            ! 2 = without MsgTag argument
    MODULE PROCEDURE messageRMatrix2            ! 3 = without MsgTag & MsgPriority argument
    MODULE PROCEDURE messageRMatrix3
    MODULE PROCEDURE messageCMatrix2
    MODULE PROCEDURE messageCMatrix3
    MODULE PROCEDURE messageIMatrix2
    MODULE PROCEDURE messageIMatrix3
    MODULE PROCEDURE message2Strings2
    MODULE PROCEDURE message2Strings3
    MODULE PROCEDURE messageLogical2      
    MODULE PROCEDURE messageLogical3
    MODULE PROCEDURE message2Strings2Integers2
    MODULE PROCEDURE message2Strings2Integers3     
    MODULE PROCEDURE messageRVector2 
    MODULE PROCEDURE messageRVector3 
    MODULE PROCEDURE messageCVector2
    MODULE PROCEDURE messageCVector3
    MODULE PROCEDURE messageIVector2
    MODULE PROCEDURE messageIVector3       
    MODULE PROCEDURE messageRScalar2          
    MODULE PROCEDURE messageRScalar3          
    MODULE PROCEDURE messageCScalar2
    MODULE PROCEDURE messageCScalar3
    MODULE PROCEDURE messageIScalar2
    MODULE PROCEDURE messageIScalar3
  END INTERFACE message

  !--------------------------------------------------------------------
  TYPE MsgPriorities
    LOGICAL :: LState               ! set to .TRUE. or .FALSE. at run-time by IWriteFLAG
    CHARACTER(50) :: SInitialMsg    ! used for priority indenting - tiered structure
  END TYPE

  TYPE MsgTags
    LOGICAL :: LState               ! set to .TRUE. or .FALSE. at run-time by IWriteFLAG
    INTEGER(IKIND) :: INumberID     ! comapred with IWriteFLAG to set state to .TRUE. or .FALSE
  END TYPE
  !--------------------------------------------------------------------

  TYPE(MsgPriorities) :: LS, LM, LL, LXL
  TYPE(MsgTags) :: no_tag, dbg7, dbg3, dbg6, dbg14, mode1

  CHARACTER*5 :: SIndentChars, SSpaces ! used for tiered-structure 
  LOGICAL,private :: LPrintThisCore ! used to print on all cores, usually = .FALSE.  
  INTEGER(IKIND) :: IClockRate ! used with print_end_timer
  
CONTAINS

  !--------------------------------------------------------------------
  ! setup message (on felix start-up, before IWriteFLAG read-in)
  !--------------------------------------------------------------------
  
  SUBROUTINE InitialiseMessage  
    ! initially set indentation and initial characters to be blank
    SSpaces = "" 
    SIndentChars   = ""
    LS%SInitialMsg  = ""; LM%SInitialMsg  = ""; LL%SInitialMsg  = ""; LXL%SInitialMsg = ""

    ! intially only set LS = .TRUE. before IWriteFLAG read-in
    LS%LState = .TRUE.; 
    LM%LState = .FALSE.;    LL%LState = .FALSE.;    LXL%LState = .FALSE.
    no_tag%LState=.FALSE.;  dbg7%LState = .FALSE.;  dbg3%LState = .FALSE.;
    dbg6%LState = .FALSE.;  dbg14%LState = .FALSE.; mode1%LState = .FALSE.
    LPrintThisCore = .FALSE.   
  END SUBROUTINE

  !--------------------------------------------------------------------
  ! read-in IWriteFLAG and setup message
  !--------------------------------------------------------------------

  SUBROUTINE SetMessageMode ( IWriteFLAG, IErr )
    
    INTEGER(IKIND), INTENT(out):: IErr
    INTEGER(IKIND), INTENT(IN) :: IWriteFLAG
    INTEGER(IKIND) :: IPrio ! short for priority, acts as IWriteFLAG
    IPrio = IWriteFLAG    
  
    !--------------------------------------------------------------------
    ! Use IWriteFLAG tens component to turn-on a MsgTag
    !--------------------------------------------------------------------

    SELECT CASE (IPrio) 
      CASE (0:9)
        no_tag%LState=.FALSE.  
      CASE (10:19)
        mode1%LState = .TRUE.
        IPrio = IPrio - 10
        SSpaces = " "
        LS%SInitialMsg  = "  "; LM%SInitialMsg  = "     "; 
        LL%SInitialMsg  = "        "; LXL%SInitialMsg = "            "
      CASE (70:79)
        dbg7%LState = .TRUE.
        IPrio = IPrio - 70
      CASE (30:39)
        dbg3%LState = .TRUE.
        IPrio = IPrio - 30
      CASE (60:69)
        dbg6%LState = .TRUE.
        IPrio = IPrio - 60
      CASE (90:99)
        dbg14%LState = .TRUE.
        IPrio = IPrio - 90
      CASE DEFAULT  ! error (minor) - tens component not recognised
        IErr=1
        IF(l_alert(IErr,"set_message_mod_mode",&
              "IWriteFLAG tens component not recognised check felix.inp")) RETURN
    END SELECT
    ! After the above, IPrio is now just a single digit INTEGER

    ! switch lower priorities on depending upon IWriteFLAG digit compoent
    IF (IPrio >= 2) LM%LState = .TRUE.;
    IF (IPrio >= 4) LL%LState = .TRUE.;
    IF (IPrio >= 8) LXL%LState = .TRUE.;

  END SUBROUTINE

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE AllowMessageOnThisCore
    LPrintThisCore = .TRUE.
    CALL messageIScalar3("MESSAGE FROM THIS CORE AS WELL, rank =",my_rank)
  END SUBROUTINE  !?? not currently used, but may be useful

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !---------------------------------------------------------------------
  ! PrintEndTime
  !---------------------------------------------------------------------  

  ! compare start time to current and print time-passed
  SUBROUTINE PrintEndTime(MsgPriority, IStartTime, STaskName)
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: STaskName
    INTEGER(IKIND), INTENT(IN) :: IStartTime
    INTEGER(IKIND) :: Ihours,Iminutes,Iseconds,ICurrentTime,INameLength
    REAL(RKIND) :: Rduration
    CHARACTER(100) :: SPrintString,Sfmt

    INameLength=LEN(STaskName)+14
    WRITE(Sfmt,*) '(A',INameLength,',A1,I4,A5,I3,A6,I3,A4)'
    Sfmt=TRIM(ADJUSTL(Sfmt))
    CALL system_clock(ICurrentTime)
    ! converts ticks from system clock into seconds
    Rduration = REAL(ICurrentTime-IStartTime)/REAL(IClockRate)
    Ihours = FLOOR(Rduration/3600.0d0)
    Iminutes = FLOOR(MOD(Rduration,3600.0d0)/60.0d0)
    Iseconds = INT(MOD(Rduration,3600.0d0)-Iminutes*60)
    WRITE(SPrintString,FMT=Sfmt) STaskName," completed in ",Ihours," hrs ",Iminutes," mins ",Iseconds," sec"
    SPrintString=TRIM(ADJUSTL(SPrintString))
    CALL message1String2(MsgPriority,SPrintString)

  END SUBROUTINE PrintEndTime

  !   EXAMPLE USAGE - use intrinisic system_clock, 'IStartTime2' local variable 
  !
  !   CALL SYSTEM_CLOCK( IStartTime2 )
  !   CALL Absorption (IErr)
  !   IF(l_alert(IErr,"felixrefine","Absorption")) CALL abort
  !   CALL print_end_time( LM, IStartTime2, "Absorption" )
  !
  !   EXAMPLE TERMINAL OUTPUT:
  !   ----> Absorption completed in   0 hrs  0 mins  2 sec

  !---------------------------------------------------------------------
  ! Main simple message
  !---------------------------------------------------------------------

  SUBROUTINE message1String ( MsgPriority, MsgTag, SString )

    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SString

    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      WRITE(*,'(a,a)') TRIM(MsgPriority%SInitialMsg)//SSpaces, TRIM(SString)
    END IF
  
  END SUBROUTINE message1String

  !---------------------------------------------------------------------
  ! Main REAL/complex/INTEGER matrix printing - used by vector and scalar printing
  !---------------------------------------------------------------------

  SUBROUTINE messageRMatrix ( MsgPriority, MsgTag, SMainMsg, RMatrix )
    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND),DIMENSION(:,:),INTENT(IN) :: RMatrix
    INTEGER(IKIND) :: i
    CHARACTER(1000) :: SFormatting, SPrintString

    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN

      IF( SIZE(RMatrix,1).EQ.1 .AND. SIZE(RMatrix,2).LE.4 ) THEN
        WRITE(SFormatting,'(a,i3,a)') '(',SIZE(RMatrix,2),'(1x,F6.2))'
        WRITE(SPrintString,SFormatting) RMatrix(:,:)
        CALL message1String ( MsgPriority, MsgTag, SMainMsg//TRIM(SPrintString) )
      ELSE
        WRITE(SFormatting,'(a,i3,a)') '(',SIZE(RMatrix,2),'(1x,F6.2))'
        CALL message1String ( MsgPriority, MsgTag, SMainMsg )
        DO i =1,SIZE(RMatrix,1)
          WRITE(SPrintString,SFormatting) RMatrix(i,:)
          CALL message1String ( MsgPriority, MsgTag, TRIM(SPrintString) )
        END DO
      END IF

    END IF

  END SUBROUTINE messageRMatrix

  ! This currently outputs in a form suitable for mathematica
  SUBROUTINE messageCMatrix ( MsgPriority, MsgTag, SMainMsg, CMatrix )
    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND),DIMENSION(:,:),INTENT(IN) :: CMatrix
    INTEGER(IKIND) :: i
    CHARACTER(10000) :: SFormatting, SPrintString

    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN

      IF(SIZE(CMatrix,2).GE.2) THEN ! multiple columns
        WRITE(SFormatting,'(a,i3,a)') '("{"',SIZE(CMatrix,2)-1,'(F10.4,1x,sp,F10.4,"i,",4x),(F10.4,1x,sp,F10.4,"i",4x)"},")'
      ELSE ! one column matrix / scalar
        WRITE(SFormatting,'(a)') '(F10.4,1x,sp,F10.4,"i",4x)'
      END IF

      IF( SIZE(CMatrix,1).EQ.1 .AND. SIZE(CMatrix,2).LE.6 ) THEN    
        WRITE(SPrintString,SFormatting) CMatrix(:,:)  
        CALL message1String ( MsgPriority, MsgTag, SMainMsg//TRIM(SPrintString) )
      ELSE
        CALL message1String ( MsgPriority, MsgTag, SMainMsg )
        DO i =1,SIZE(CMatrix,1)
          WRITE(SPrintString,SFormatting) CMatrix(i,:)
          CALL message1String ( MsgPriority, MsgTag, TRIM(SPrintString) )
        END DO
      END IF

    END IF

  END SUBROUTINE messageCMatrix

  SUBROUTINE messageIMatrix ( MsgPriority, MsgTag, SMainMsg, IMatrix )
    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND),DIMENSION(:,:),INTENT(IN) :: IMatrix
    INTEGER(IKIND) :: i
    CHARACTER(1000) :: SFormatting, SPrintString

    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN

      IF( SIZE(IMatrix,1).EQ.1 .AND. SIZE(IMatrix,2).LE.15 ) THEN
        WRITE(SFormatting,'(a,i3,a)') '(',SIZE(IMatrix,2),'(1x,i4))'
        WRITE(SPrintString,SFormatting) IMatrix(:,:)
        CALL message1String ( MsgPriority, MsgTag, SMainMsg//TRIM(SPrintString) )
      ELSE
        WRITE(SFormatting,'(a,i3,a)') '(',SIZE(IMatrix,2),'(i4,1x))'
        CALL message1String ( MsgPriority, MsgTag, SMainMsg )
        DO i =1,SIZE(IMatrix,1)
          WRITE(SPrintString,SFormatting) IMatrix(i,:)
          CALL message1String ( MsgPriority, MsgTag, TRIM(SPrintString) )
        END DO
      END IF

    END IF

  END SUBROUTINE messageIMatrix

  !---------------------------------------------------------------------
  ! Other key message types
  !---------------------------------------------------------------------

  SUBROUTINE message2Strings ( MsgPriority, MsgTag, SString1, SString2 )

    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SString1, SString2

    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      WRITE(*,'(a,a,a,a,a)') TRIM(MsgPriority%SInitialMsg)//SSpaces, SString1,"'",SString2,"'"
    END IF
  
  END SUBROUTINE message2Strings

  SUBROUTINE messageLogical ( MsgPriority, MsgTag, SMainMsg, LLogicalVariable )  
  
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    LOGICAL, INTENT(IN) :: LLogicalVariable

    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      WRITE(*,'(a,a,1x,l1)') TRIM(MsgPriority%SInitialMsg)//SSpaces, SMainMsg, LLogicalVariable
    END IF

  END SUBROUTINE messageLogical

  ! text, integer, extra-text, extra-integer
  SUBROUTINE message2Strings2Integers ( MsgPriority, MsgTag, text1, int1, text2, int2)

    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: text1, text2
    INTEGER(IKIND), INTENT(IN) :: int1, int2

    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      WRITE(*,'(a,a,i4.3,a,i4.3)') TRIM(MsgPriority%SInitialMsg)//SSpaces, text1, int1, text2, int2
    END IF
  
  END SUBROUTINE message2Strings2Integers

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  ! MAIN INTERFACES BELOW
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! vector and scalar interfaces
  !---------------------------------------------------------------------

  SUBROUTINE messageRVector ( MsgPriority, MsgTag, SMainMsg, RVector )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND),DIMENSION(:),INTENT(IN) :: RVector
    CALL messageRMatrix( MsgPriority, MsgTag, SMainMsg, RESHAPE(RVector,[1,SIZE(RVector,1)]) )
  END SUBROUTINE messageRVector

  SUBROUTINE messageCVector ( MsgPriority, MsgTag, SMainMsg, CVector )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND),DIMENSION(:),INTENT(IN) :: CVector
    CALL messageCMatrix( MsgPriority, MsgTag, SMainMsg, RESHAPE(CVector,[1,SIZE(CVector,1)]) )
  END SUBROUTINE messageCVector

  SUBROUTINE messageIVector ( MsgPriority, MsgTag, SMainMsg, IVector )  
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND),DIMENSION(:),INTENT(IN) :: IVector
    CALL messageIMatrix( MsgPriority, MsgTag, SMainMsg, RESHAPE(IVector,[1,SIZE(IVector,1)]) )
  END SUBROUTINE messageIVector

  SUBROUTINE messageRScalar ( MsgPriority, MsgTag, SMainMsg, RScalar )
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND), INTENT(IN) :: RScalar
    CALL messageRMatrix ( MsgPriority, MsgTag, SMainMsg, RESHAPE([RScalar],[1,1]) )  
  END SUBROUTINE messageRScalar

  SUBROUTINE messageCScalar ( MsgPriority, MsgTag, SMainMsg, CScalar )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND),INTENT(IN) :: CScalar
    CALL messageCMatrix ( MsgPriority, MsgTag, SMainMsg, RESHAPE([CScalar],[1,1]) )
  END SUBROUTINE messageCScalar

  SUBROUTINE messageIScalar ( MsgPriority, MsgTag, SMainMsg, IScalar )
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND), INTENT(IN) :: IScalar
    CALL messageIMatrix ( MsgPriority, MsgTag, SMainMsg, RESHAPE([IScalar],[1,1]) )
  END SUBROUTINE messageIScalar

  !---------------------------------------------------------------------
  ! optional argument interfaces
  !---------------------------------------------------------------------

  SUBROUTINE message1String2 ( MsgPriority, SMainMsg )
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    CALL message1String ( MsgPriority, no_tag, SMainMsg )  
  END SUBROUTINE message1String2

  SUBROUTINE message1String3 ( SMainMsg )
    CHARACTER(*), INTENT(IN) :: SMainMsg
    CALL message1String ( LS, no_tag, SMainMsg )  
  END SUBROUTINE message1String3

  SUBROUTINE messageRMatrix2 ( MsgPriority, SMainMsg, RMatrix )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND),DIMENSION(:,:),INTENT(IN) :: RMatrix
    CALL messageRMatrix ( MsgPriority, no_tag, SMainMsg, RMatrix )
  END SUBROUTINE messageRMatrix2

  SUBROUTINE messageRMatrix3 ( SMainMsg, RMatrix )    
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND),DIMENSION(:,:),INTENT(IN) :: RMatrix
    CALL messageRMatrix ( LS, no_tag, SMainMsg, RMatrix )
  END SUBROUTINE messageRMatrix3  

  SUBROUTINE messageCMatrix2 ( MsgPriority, SMainMsg, cmatrix )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND), DIMENSION(:,:), INTENT(IN) :: cmatrix
    CALL messageCMatrix ( MsgPriority, no_tag, SMainMsg, cmatrix )
  END SUBROUTINE messageCMatrix2

  SUBROUTINE messageCMatrix3 ( SMainMsg, cmatrix )    
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND), DIMENSION(:,:), INTENT(IN) :: cmatrix
    CALL messageCMatrix ( LS, no_tag, SMainMsg, cmatrix )
  END SUBROUTINE messageCMatrix3

  SUBROUTINE messageIMatrix2 ( MsgPriority, SMainMsg, imatrix )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND),DIMENSION(:,:),INTENT(IN) :: imatrix
    CALL messageIMatrix ( MsgPriority, no_tag, SMainMsg, imatrix )
  END SUBROUTINE messageIMatrix2

  SUBROUTINE messageIMatrix3 ( SMainMsg, imatrix )    
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND),DIMENSION(:,:),INTENT(IN) :: imatrix
    CALL messageIMatrix ( LS, no_tag, SMainMsg, imatrix )
  END SUBROUTINE messageIMatrix3

  SUBROUTINE message2Strings2 ( MsgPriority, SString1, SString2 )
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SString1, SString2
    CALL message2Strings ( MsgPriority, no_tag, SString1, SString2 )
  END SUBROUTINE message2Strings2

  SUBROUTINE message2Strings3 ( SString1, SString2 )
    CHARACTER(*), INTENT(IN) :: SString1, SString2
    CALL message2Strings ( LS, no_tag, SString1, SString2 )
  END SUBROUTINE message2Strings3

  SUBROUTINE messageLogical2 ( MsgPriority, SMainMsg, LLogicalVariable )   
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    LOGICAL, INTENT(IN) :: LLogicalVariable
    CALL messageLogical ( MsgPriority, no_tag, SMainMsg, LLogicalVariable )
  END SUBROUTINE messageLogical2

  SUBROUTINE messageLogical3 ( SMainMsg, LLogicalVariable )   
    CHARACTER(*), INTENT(IN) :: SMainMsg
    LOGICAL, INTENT(IN) :: LLogicalVariable
    CALL messageLogical ( LS, no_tag, SMainMsg, LLogicalVariable )
  END SUBROUTINE messageLogical3

  SUBROUTINE message2Strings2Integers2 ( MsgPriority, text1, int1, text2, int2)
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: text1, text2
    INTEGER(IKIND), INTENT(IN) :: int1, int2
    CALL message2Strings2Integers ( MsgPriority, no_tag, text1, int1, text2, int2)
  END SUBROUTINE message2Strings2Integers2

  SUBROUTINE message2Strings2Integers3 ( text1, int1, text2, int2)
    CHARACTER(*), INTENT(IN) :: text1, text2
    INTEGER(IKIND), INTENT(IN) :: int1, int2
    CALL message2Strings2Integers ( LS, no_tag, text1, int1, text2, int2)
  END SUBROUTINE message2Strings2Integers3

  SUBROUTINE messageRVector2 ( MsgPriority, SMainMsg, RVector )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND),DIMENSION(:),INTENT(IN) :: RVector
    CALL messageRVector ( MsgPriority, no_tag, SMainMsg, RVector )
  END SUBROUTINE messageRVector2

  SUBROUTINE messageRVector3 ( SMainMsg, RVector )    
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND),DIMENSION(:),INTENT(IN) :: RVector
    CALL messageRVector ( LS, no_tag, SMainMsg, RVector )
  END SUBROUTINE messageRVector3

  SUBROUTINE messageCVector2 ( MsgPriority, SMainMsg, CVector )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND), DIMENSION(:), INTENT(IN) :: CVector
    CALL messageCVector ( MsgPriority, no_tag, SMainMsg, CVector )
  END SUBROUTINE messageCVector2

  SUBROUTINE messageCVector3 ( SMainMsg, CVector )    
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND), DIMENSION(:), INTENT(IN) :: CVector
    CALL messageCVector( LS, no_tag, SMainMsg, CVector )
  END SUBROUTINE messageCVector3

  SUBROUTINE messageIVector2 ( MsgPriority, SMainMsg, IVector )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND),DIMENSION(:),INTENT(IN) :: IVector
    CALL messageIVector ( MsgPriority, no_tag, SMainMsg, IVector )
  END SUBROUTINE messageIVector2

  SUBROUTINE messageIVector3 ( SMainMsg, IVector )    
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND),DIMENSION(:),INTENT(IN) :: IVector
    CALL messageIVector ( LS, no_tag, SMainMsg, IVector )
  END SUBROUTINE messageIVector3

  SUBROUTINE messageRScalar2 ( MsgPriority, SMainMsg, RScalar )
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND), INTENT(IN) :: RScalar
    CALL messageRScalar ( MsgPriority, no_tag, SMainMsg, RScalar )  
  END SUBROUTINE messageRScalar2

  SUBROUTINE messageRScalar3 ( SMainMsg, RScalar )
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND), INTENT(IN) :: RScalar
    CALL messageRScalar ( LS, no_tag, SMainMsg, RScalar )  
  END SUBROUTINE messageRScalar3

  SUBROUTINE messageCScalar2 ( MsgPriority, SMainMsg, CScalar )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND),INTENT(IN) :: CScalar
    CALL messageCScalar ( MsgPriority, no_tag, SMainMsg, CScalar )
  END SUBROUTINE messageCScalar2

  SUBROUTINE messageCScalar3 ( SMainMsg, CScalar ) 
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND),INTENT(IN) :: CScalar
    CALL messageCScalar ( LS, no_tag, SMainMsg, CScalar )
  END SUBROUTINE messageCScalar3

  SUBROUTINE messageIScalar2 ( MsgPriority, SMainMsg, IScalar )
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND), INTENT(IN) :: IScalar
    CALL messageIScalar ( MsgPriority, no_tag, SMainMsg, IScalar )  
  END SUBROUTINE messageIScalar2

  SUBROUTINE messageIScalar3 ( SMainMsg, IScalar )
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND), INTENT(IN) :: IScalar
    CALL messageIScalar ( LS, no_tag, SMainMsg, IScalar )  
  END SUBROUTINE messageIScalar3

END MODULE message_mod
