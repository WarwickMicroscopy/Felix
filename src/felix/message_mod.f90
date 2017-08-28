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
!>
!>  MODULE OVERVIEW:
!>
!>    This MODULE directly CONTAINS everything for message,
!>      which is a major SUBROUTINE used throughout felix. It replaces the FORTRAN print
!>      statement and provides formatted output chosen at run-time.
!>
!>    message_mod also accesses l_alert_mod:
!>      l_alert is a simple but important function, used extensively for error handling.
!>
!>    message_mod also conatins print_end_time used to neatly print the time elapsed.
!>      
!>
!>--------------------------------------------------------------------------------------------
!>
!>  message( ... ) FEATURES:
!>
!>    IWriteFLAG in felix.inp is the run-time input used to select the output mode of Felix.
!>    message uses IWriteFLAG to decide whether a message should be printed.
!>
!>    - message prints various different variable TYPEs with consistent formatting
!>    - It accommodates for MPI parallel programming used in Felix by printing on one 'core'
!>    - It prints a tiered-structure to distinguish key messages from more niche ones
!>    - It speeds up debugging and testing, by being very quick and easy to use
!>                    
!>--------------------------------------------------------------------------------------------
!>
!>  message( MsgPriority, MsgTag, text_to_print, variable_to_print )
!>
!>    MsgPriority     
!>    (optional)
!>
!>            DERIVED TYPE ( msg_priority_TYPE ) = LS, LM, LL or LXL
!>
!>            This determines how important it is to print this message. Is it 
!>            stating key operations of Felix or more particular, finer details?              
!>             
!>            LS  : high priority (small level of output) for important statements & values
!>            LM  : medium prioriity (medium level of output)
!>            LL  : low priority (large output) for fines details & thorough variable output
!>            LXL : extremely low priority large output) & niche details
!>
!>            DEFAULT = LS
!>
!>            At run-time this priority will be compared against the selected output mode
!>            given by IWriteFLAG in felix.inp. This will decide whether a message is
!>            is important enough to be printed. Additionally, Low priority  
!>            messages are indented further to the right to produce a tiered-output
!>            for clarity and readability.  
!>
!>    MsgTag       
!>    (optional)
!>              
!>            DERIVED TYPE ( msg_tag_TYPE )
!>              
!>            This allows for messages to be grouped together with tags. At run-time you 
!>            can select that only key messages and messages with a certain tag are printed.
!>
!>    text_to_print
!>     
!>            TYPE: String (1D CHARACTER array)
!>
!>            This is some text to print to the terminal.
!>            Any variable values will be printed after this text. 
!>                        
!>    variable_to_print 
!>    (optional)
!> 
!>            TYPE: complex, INTEGER, REAL, LOGICAL, String 
!>            DIMENSION: scalar, vector or matrix
!>
!>            The value of this variable will be printed to the terminal with
!>            a format depending upon its TYPE and DIMENSION. E.g. a REAL matrix
!>            will be printed on multiple lines in columns with scientific notation.                
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
!>        0   - no_tag  : no specific MsgTag's are switched on, only priority is considered
!>        30  - dbg3    : corresponds to IWriteFLAG = 3 IN THE OLD SYSTEM, many print
!>                        statements in the code previously only printed when IWriteFLAG = 3
!>
!>        For the full list of MsgTag's and their corresponding IWriteFLAG ten's value,
!>        refer to init_message_LOGICALs below.
!>
!>    IWriteFLAG examples: 
!>        2  will print all LS, LM messages
!>        32 will print all LS, LM, and dbg3 tagged messages
!>        96 will print all LS, LM, LL, and dbg14 tagged messages 
!>
!>--------------------------------------------------------------------------------------------
!>
!>  message( ... ) EXAMPLES  
!>
!>--------------------------------------------------------------------------------------------

!   TO-DO LIST (SCRUFFY):

!?? fix message matrices
!?? implicit no_tag & private
!?? consider current procedure state
!?? currently using derived TYPE to force particular priority/dbg use
!?? standardise format top level, SSpaces and INTEGERs...
!?? consider all variables as 2 dimensional matrices?
!?? consider could make dbg initialise select case more concise with array of derived
!?? rearrange to have terminal_error modules, msg_output_group examples, 
!?? rename to MsgPriority, update MsgPriority for clarity high vs. low priority
!?? priority for top level refinement details
!?? adding new MsgTags
!?? generic -> optional arguments not implimented...


!---------------------------------------------------------------------------------------------

  USE l_alert_mod             ! grants access to felix's main error handling
  USE MyNumbers               ! necesary for IKIND, RKIND etc.
  USE MyMPI, ONLY : my_rank   ! necesary for my_rank - used to print on core 0 only

  !?? JR this grants MyNumbers, USE MyNumbers throughout code is for clariry...

!---------------------------------------------------------------------------------------------

! HUGE INTERFACE to handle various variables TYPEs and optional arguments

  INTERFACE message

    MODULE PROCEDURE message1String 

    MODULE PROCEDURE messageRMatrix        
    MODULE PROCEDURE messageCMatrix       
    MODULE PROCEDURE messageIMatrix

    MODULE PROCEDURE message2Strings 
    MODULE PROCEDURE messageLogical  
      
    ! special versions of message
    MODULE PROCEDURE message_alongside_ir   ! prints INT & REAL matrix alongside
    MODULE PROCEDURE message_alongside_ir2  ! prints INT matrix & REAL vector alongside
    MODULE PROCEDURE message_alongside_ic   ! prints INT & complex matrix alongside
    MODULE PROCEDURE messageInteger_isi     ! prints text, INTEGER, more text, INTEGER 

    ! optional argument interfaces for special versions of message
    MODULE PROCEDURE messageInteger_isi2

    ! interfaces used for scalars, vectors and optional arguments
    MODULE PROCEDURE messageRVector        
    MODULE PROCEDURE messageCVector         ! suffix :
    MODULE PROCEDURE messageIVector         ! 2 = without MsgTag argument
    MODULE PROCEDURE messageReal            ! 3 = without MsgTag & MsgPriority argument  
    MODULE PROCEDURE messageReal2          
    MODULE PROCEDURE messageReal3          
    MODULE PROCEDURE messageRVector2 
    MODULE PROCEDURE messageRVector3
    MODULE PROCEDURE messageRMatrix2
    MODULE PROCEDURE messageRMatrix3
    MODULE PROCEDURE messageComplex
    MODULE PROCEDURE messageComplex2
    MODULE PROCEDURE messageComplex3
    MODULE PROCEDURE messageCVector2
    MODULE PROCEDURE messageCVector3
    MODULE PROCEDURE messageCMatrix2
    MODULE PROCEDURE messageCMatrix3
    MODULE PROCEDURE messageInteger
    MODULE PROCEDURE messageInteger2
    MODULE PROCEDURE messageInteger3
    MODULE PROCEDURE messageIVector2
    MODULE PROCEDURE messageIVector3
    MODULE PROCEDURE messageIMatrix2
    MODULE PROCEDURE messageIMatrix3
    MODULE PROCEDURE message2Strings2
    MODULE PROCEDURE message2Strings3
    MODULE PROCEDURE message1String2
    MODULE PROCEDURE message1String3

  END INTERFACE message

  !--------------------------------------------------------------------

  TYPE MsgPriorities
    LOGICAL :: LState               ! set to .true. or .false. at run-time by IWriteFLAG
    CHARACTER(50) :: SInitialMsg    ! used for priority indenting - tiered structure
  END TYPE

  TYPE MsgTags
    LOGICAL :: LState               ! set to .true. or .false. at run-time by IWriteFLAG
    INTEGER(IKIND) :: INumberID     ! comapred with IWriteFLAG to set state to .true. or .false
  END TYPE

  !--------------------------------------------------------------------

  TYPE(MsgPriorities) :: LS, LM, LL, LXL

  TYPE(MsgTags) :: no_tag, dbg7, dbg3, dbg6, dbg14

  LOGICAL,private :: LPrintThisCore ! used to print on all cores, usually = .false.
  
  CHARACTER(:), allocatable :: SIndentSpaces, SSpaces ! used for tiered-structure 

  INTEGER(IKIND) :: IClockRate ! used with print_end_timer
  
CONTAINS

  !--------------------------------------------------------------------
  ! setup message (on felix start-up, before IWriteFLAG read-in)
  !--------------------------------------------------------------------
  
  SUBROUTINE init_message_logicals  
 
    !--------------------------------------------------------------------
    ! MsgTags fixed ID numbers (to compare against IWriteFLAG)
    !--------------------------------------------------------------------

    no_tag%INumberID  = 0
    dbg7%INumberID    = 70
    dbg3%INumberID    = 30
    dbg6%INumberID    = 60
    dbg14%INumberID   = 90 ! NB would use 140, but our integer kind(1) cannot reach 140

    ! NB to add another MsgTag, declare it above, add it here and add it to set_message_mod_mode below

    !--------------------------------------------------------------------
    ! setup tiered-structure
    !--------------------------------------------------------------------

    ! e.g. WRITE(*,formatting) TRIM(MsgPriority%SInitialMsg)//SSpaces, SMainMsg, RVector

    SIndentSpaces   = "" ! leading spaces for all messages

    ! tiered-structure MsgPriority-specific initial message
    LS%SInitialMsg  = SIndentSpaces               
    LM%SInitialMsg  = SIndentSpaces//"   "
    LL%SInitialMsg  = SIndentSpaces//"     "        
    LXL%SInitialMsg = SIndentSpaces//"       "
    ! NB SInitialMsg is TRIM'ed so trailing spaces are pointless 

    SSpaces = " " 
    ! SSpaces are the global SSpaces after SInitialMsg

    !--------------------------------------------------------------------
    ! intially only set LS = .true. before IWriteFLAG read-in
    !--------------------------------------------------------------------

    LS%LState = .true.; 

    LM%LState = .false.;    LL%LState = .false.;   LXL%LState = .false.
    no_tag%LState=.false.;  dbg7%LState = .false.; dbg3%LState = .false.;
    dbg6%LState = .false.;  dbg14%LState = .false.
    LPrintThisCore = .false.   

  END SUBROUTINE

  !--------------------------------------------------------------------
  ! read-in IWriteFLAG and setup message
  !--------------------------------------------------------------------

  SUBROUTINE set_message_mod_mode ( IWriteFLAG, IErr )
    
    INTEGER(IKIND), INTENT(out):: IErr
    INTEGER(IKIND), INTENT(IN) :: IWriteFLAG
    INTEGER(IKIND) :: IPrio ! short for priority, acts as IWriteFLAG
    IPrio = IWriteFLAG    
  
    !--------------------------------------------------------------------
    ! Use IWriteFLAG tens component to turn-on a MsgTag
    !--------------------------------------------------------------------

    IF ( IPrio > 10 ) THEN

      IF ( (dbg7%INumberID <= IPrio) .and. (IPrio  <= dbg7%INumberID+9) ) THEN
        dbg7%LState = .true.
        IPrio = IPrio - dbg7%INumberID
      ELSEIF ( (dbg3%INumberID <= IPrio) .and. (IPrio  <= dbg3%INumberID+9) ) THEN
        dbg3%LState = .true.
        IPrio = IPrio - dbg3%INumberID
      ELSEIF ( (dbg6%INumberID <= IPrio) .and. (IPrio  <= dbg6%INumberID+9) ) THEN
        dbg6%LState = .true.
        IPrio = IPrio - dbg6%INumberID
      ELSEIF ( (dbg14%INumberID <= IPrio) .and. (IPrio  <= dbg14%INumberID+9) ) THEN
        dbg14%LState = .true.
        IPrio = IPrio - dbg14%INumberID
      ELSE  ! error (minor) - tens component not recognised
        IErr=1
        IF(l_alert(IErr,"set_message_mod_mode",&
              "IWriteFLAG tens component not recognised check felix.inp")) RETURN
      END IF

    END IF
    ! After the above, IPrio is now just a single digit INTEGER

    ! switch lower priorities on depending upon IWriteFLAG digit compoent
    IF (IPrio >= 2) LM%LState = .true.;
    IF (IPrio >= 4) LL%LState = .true.;
    IF (IPrio >= 8) LXL%LState = .true.;

  END SUBROUTINE

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE allow_message_on_this_core
    LPrintThisCore = .true.
    CALL messageInteger3("MESSAGE FROM THIS CORE AS WELL, rank =",my_rank)
  END SUBROUTINE  !?? not currently used, but may be useful

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !---------------------------------------------------------------------
  ! print_end_time
  !---------------------------------------------------------------------  

  ! compare start time to current and print time-passed
  SUBROUTINE print_end_time(MsgPriority, Istart_time, STaskName)
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: STaskName
    INTEGER(IKIND), INTENT(IN) :: Istart_time
    INTEGER(IKIND) :: Ihours,Iminutes,Iseconds,ICurrentTime,INameLength
    REAL(RKIND) :: Rduration
    CHARACTER(100) :: SPrintString,Sfmt

    INameLength=LEN(STaskName)+14
    WRITE(Sfmt,*) '(A',INameLength,',A1,I3,A5,I2,A6,I2,A4)'
    Sfmt=TRIM(ADJUSTL(Sfmt))
    CALL system_clock(ICurrentTime)
    ! converts ticks from system clock into seconds
    Rduration = REAL(ICurrentTime-Istart_time)/REAL(IClockRate)
    Ihours = FLOOR(Rduration/3600.0d0)
    Iminutes = FLOOR(MOD(Rduration,3600.0d0)/60.0d0)
    Iseconds = INT(MOD(Rduration,3600.0d0)-Iminutes*60)
    WRITE(SPrintString,FMT=Sfmt) STaskName," completed in ",Ihours," hrs ",Iminutes," mins ",Iseconds," sec"
    SPrintString=TRIM(ADJUSTL(SPrintString))
    CALL message(MsgPriority,SPrintString)

  END SUBROUTINE print_end_time

  !   EXAMPLE USAGE - use intrinisic system_clock, 'IStartTime2' local variable 
  !
  !   CALL SYSTEM_CLOCK( IStartTime2 )
  !   CALL Absorption (IErr)
  !   IF(l_alert(IErr,"felixrefine","Absorption")) CALL abort
  !   CALL print_end_time( LM, IStartTime2, "Absorption" )
  !
  !   EXAMPLE TERMINAL OUTPUT:
  !   @ ---- Absorption completed in   0 hrs  0 mins  2 sec

  !---------------------------------------------------------------------
  ! main REAL/complex/INTEGER matrix printing - used by vector and scalar printing
  !---------------------------------------------------------------------

  SUBROUTINE message1String ( MsgPriority, MsgTag, SString )

    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SString

    ! check priority then print SMainMsg
    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      WRITE(*,'(a,a)') TRIM(MsgPriority%SInitialMsg)//SSpaces, SString
    END IF
  
  END SUBROUTINE message1String

  SUBROUTINE messageRMatrix ( MsgPriority, MsgTag, SMainMsg, RMatrix )
    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND),DIMENSION(:,:),INTENT(IN) :: RMatrix
    INTEGER(IKIND) :: i
    CHARACTER(1000) :: SFormatting, SPrintString

    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      CALL message1String ( MsgPriority, MsgTag, SMainMsg )
      SPrintString = ''
      WRITE(SFormatting,'(a,i3,a)') '(',SIZE(RMatrix,2),'(1x,sp,ES10.3))'
      DO i =1,SIZE(RMatrix,1)
        WRITE(SPrintString,SFormatting) RMatrix(i,:)
        CALL message1String ( MsgPriority, MsgTag, TRIM(SPrintString) )
      END DO
    END IF

  END SUBROUTINE messageRMatrix

  SUBROUTINE messageCMatrix ( MsgPriority, MsgTag, SMainMsg, CMatrix )
    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND),DIMENSION(:,:),INTENT(IN) :: CMatrix ! complex martrix
    INTEGER(IKIND) :: i

    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      CALL message1String ( MsgPriority, MsgTag, SMainMsg )
      SPrintString = ''
      WRITE(SFormatting,'(a,i3,a)') '(',SIZE(CMatrix,2),'(1x,sp,ES8.1,1x,ES8.1,"i"))'
      DO i =1,SIZE(CMatrix,1)
        WRITE(SPrintString,SFormatting) CMatrix(i,:)
        CALL message1String ( MsgPriority, MsgTag, TRIM(SPrintString) )
      END DO
    END IF

  END SUBROUTINE messageCMatrix

  SUBROUTINE messageIMatrix ( MsgPriority, MsgTag, SMainMsg, IMatrix )
    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND),DIMENSION(:,:),INTENT(IN) :: IMatrix ! INTEGER martrix
    INTEGER(IKIND) :: i

    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      CALL message1String ( MsgPriority, MsgTag, SMainMsg )
      SPrintString = ''
      WRITE(SFormatting,'(a,i3,a)') '(',SIZE(IMatrix,2),'(1x,i4.3))'
      DO i =1,SIZE(IMatrix,1)
        WRITE(SPrintString,SFormatting) IMatrix(i,:)
        CALL message1String ( MsgPriority, MsgTag, TRIM(SPrintString) )
      END DO
    END IF

  END SUBROUTINE messageIMatrix

  !---------------------------------------------------------------------
  ! main REAL/complex/INTEGER vector printing - used by matrix and scalar printing
  !---------------------------------------------------------------------





  SUBROUTINE message2Strings ( MsgPriority, MsgTag, SMainMsg, str_variable )

    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg, str_variable

    ! check priority then print SMainMsg & String
    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      WRITE(*,'(a,a,a,a,a)') TRIM(MsgPriority%SInitialMsg)//SSpaces, SMainMsg,"'",str_variable,"'"
    END IF
  
  END SUBROUTINE message2Strings

  ! For printing LOGICALs ! currently no interfaces for optional arguments, matrices, vectors
  SUBROUTINE messageLogical ( MsgPriority, MsgTag, SMainMsg, LOGICAL_var )  
  
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    LOGICAL, INTENT(IN) :: LOGICAL_var ! LOGICAL variable

    ! check priority then print LOGICAL variable
    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      WRITE(*,'(a,a,1x,l1)') TRIM(MsgPriority%SInitialMsg)//SSpaces, SMainMsg, LOGICAL_var
    END IF

  END SUBROUTINE messageLogical



  !---------------------------------------------------------------------
  ! Special message varieties
  !---------------------------------------------------------------------
  
  ! vertically INT matrix beside REAL matrix
  SUBROUTINE message_alongside_ir ( MsgPriority, MsgTag, &
               SMainMsg, imatrix1, rmatrix2 )
    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND),DIMENSION(:,:),INTENT(IN) :: imatrix1 ! 1st matrix - INTEGER
    REAL(RKIND), DIMENSION(:,:), INTENT(IN) :: rmatrix2 ! 2nd matrix - REAL
    CHARACTER(100) :: formatting
    INTEGER(IKIND) :: i

    ! check priority then print two matrices alongside
    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      WRITE(*,'(a,a)') TRIM(MsgPriority%SInitialMsg)//SSpaces, SMainMsg
      WRITE(formatting,'(a,i3.3,a,i3.3,a)') '(a,i3.3,',SIZE(imatrix1,2),'(1x,i4.3),a,',&
            SIZE(rmatrix2,2),'(1x,sp,ES10.3))'
      DO i =1,SIZE(imatrix1,1)
        WRITE(*,formatting) TRIM(MsgPriority%SInitialMsg)//SSpaces,i, imatrix1(i,:)," |",rmatrix2(i,:)
      END DO
    END IF  
  
  END SUBROUTINE message_alongside_ir

  ! vertically INT matrix and REAL vector
  SUBROUTINE message_alongside_ir2 ( MsgPriority, MsgTag, &
               SMainMsg, imatrix, RVector )

    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND),DIMENSION(:,:),INTENT(IN) :: imatrix ! 1st matrix - INTEGER
    REAL(RKIND), DIMENSION(:), INTENT(IN) :: RVector ! 2nd matrix - REAL
    CHARACTER(100) :: formatting
    INTEGER(IKIND) :: i

    CALL message_alongside_ir ( MsgPriority, MsgTag, &
               SMainMsg, imatrix, RESHAPE(RVector,[SIZE(RVector),1]) ) 
  
  END SUBROUTINE message_alongside_ir2

  ! vertically INT matrix beside complex matrix
  SUBROUTINE message_alongside_ic ( MsgPriority, MsgTag, &
               SMainMsg, imatrix1, cmatrix2 )
    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND), DIMENSION(:,:),INTENT(IN) :: imatrix1 ! 1st matrix - INTEGER
    COMPLEX(CKIND), DIMENSION(:,:), INTENT(IN) :: cmatrix2 ! 2nd matrix - REAL
    CHARACTER(100) :: formatting
    INTEGER(IKIND) :: i

    ! check priority then print two matrices alongside
    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      WRITE(*,'(a,a)') TRIM(MsgPriority%SInitialMsg)//SSpaces, SMainMsg
      WRITE(formatting,'(a,i3.3,a,i3.3,a)') '(a,i3.3,',SIZE(imatrix1,2),'(1x,i4.3),a,',&
            SIZE(cmatrix2,2),'(1x,"(",sp,ES8.1," ",ES8.1,"i)"))'
      DO i =1,SIZE(imatrix1,1)
        WRITE(*,formatting) TRIM(MsgPriority%SInitialMsg)//SSpaces,i, imatrix1(i,:)," |",cmatrix2(i,:)
      END DO
    END IF  
  
  END SUBROUTINE message_alongside_ic

  ! text, INTEGER, extra-text, extra-INTEGER
  SUBROUTINE messageInteger_isi ( MsgPriority, MsgTag, text1, int1, text2, int2)

    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: text1, text2
    INTEGER(IKIND), INTENT(IN) :: int1, int2

    IF ( ( my_rank==0 .OR. LPrintThisCore ) &
    .AND. (MsgPriority%LState .OR. MsgTag%LState) ) THEN
      WRITE(*,'(a,a,i4.3,a,i4.3)') TRIM(MsgPriority%SInitialMsg)//SSpaces, text1, int1, text2, int2
    END IF
  
  END SUBROUTINE messageInteger_isi

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  ! INTERFACES BELOW
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Matrix messages using corresponding vector messages
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Variable, vector, matrix interfaces for optional (priority and tag) arguments
  !---------------------------------------------------------------------

  SUBROUTINE messageRVector ( MsgPriority, MsgTag, SMainMsg, RVector )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND),DIMENSION(:),INTENT(IN) :: RVector ! REAL vector
    CALL messageRMatrix( MsgPriority, MsgTag, RESHAPE(RVector,1))
  END SUBROUTINE messageRVector

  SUBROUTINE messageCVector ( MsgPriority, MsgTag, SMainMsg, CVector )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND),DIMENSION(:),INTENT(IN) :: CVector ! complex vector
    CALL messageCMatrix( MsgPriority, MsgTag, RESHAPE(CVector,1))
  END SUBROUTINE messageCVector

  SUBROUTINE messageIVector ( MsgPriority, MsgTag, SMainMsg, IVector )  
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND),DIMENSION(:),INTENT(IN) :: IVector ! INTEGER vector
    CALL messageIMatrix( MsgPriority, MsgTag, RESHAPE(IVector,1))
  END SUBROUTINE messageIVector

  SUBROUTINE messageReal ( MsgPriority, MsgTag, SMainMsg, real_variable )
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND), INTENT(IN) :: real_variable
    CALL messageRVector ( MsgPriority, MsgTag, SMainMsg, (/ real_variable /) )  
  END SUBROUTINE messageReal

  SUBROUTINE messageReal2 ( MsgPriority, SMainMsg, real_variable )
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND), INTENT(IN) :: real_variable
    CALL messageRVector ( MsgPriority, no_tag, SMainMsg, (/ real_variable /) )  
  END SUBROUTINE messageReal2

  SUBROUTINE messageReal3 ( SMainMsg, real_variable )
    CHARACTER(*), INTENT(IN) :: SMainMsg
    REAL(RKIND), INTENT(IN) :: real_variable
    CALL messageRVector ( LS, no_tag, SMainMsg, (/ real_variable /) )  
  END SUBROUTINE messageReal3

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

  SUBROUTINE messageComplex ( MsgPriority, MsgTag, SMainMsg, c_variable )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND),INTENT(IN) :: c_variable
    CALL messageCVector ( MsgPriority, MsgTag, SMainMsg, (/ c_variable /) )
  END SUBROUTINE messageComplex

  SUBROUTINE messageComplex2 ( MsgPriority, SMainMsg, c_variable )    
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND),INTENT(IN) :: c_variable
    CALL messageCVector ( MsgPriority, no_tag, SMainMsg, (/ c_variable /) )
  END SUBROUTINE messageComplex2

  SUBROUTINE messageComplex3 ( SMainMsg, c_variable ) 
    CHARACTER(*), INTENT(IN) :: SMainMsg
    COMPLEX(CKIND),INTENT(IN) :: c_variable
    CALL messageCVector ( LS, no_tag, SMainMsg, (/ c_variable /) )
  END SUBROUTINE messageComplex3

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

  SUBROUTINE messageInteger ( MsgPriority, MsgTag, SMainMsg, int_variable )
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    TYPE (MsgTags), INTENT(IN) :: MsgTag
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND), INTENT(IN) :: int_variable
    CALL messageIVector ( MsgPriority, MsgTag, SMainMsg, (/ int_variable /) )
  END SUBROUTINE messageInteger

  SUBROUTINE messageInteger2 ( MsgPriority, SMainMsg, int_variable )
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND), INTENT(IN) :: int_variable
    CALL messageIVector ( MsgPriority, no_tag, SMainMsg, (/ int_variable /) )  
  END SUBROUTINE messageInteger2

  SUBROUTINE messageInteger3 ( SMainMsg, int_variable )
    CHARACTER(*), INTENT(IN) :: SMainMsg
    INTEGER(IKIND), INTENT(IN) :: int_variable
    CALL messageIVector ( LS, no_tag, SMainMsg, (/ int_variable /) )  
  END SUBROUTINE messageInteger3

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

  SUBROUTINE message2Strings2 ( MsgPriority, SMainMsg, str_variable )
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg, str_variable
    CALL message2Strings ( MsgPriority, no_tag, SMainMsg, str_variable )
  END SUBROUTINE message2Strings2

  SUBROUTINE message2Strings3 ( SMainMsg, str_variable )
    CHARACTER(*), INTENT(IN) :: SMainMsg, str_variable
    CALL message2Strings ( LS, no_tag, SMainMsg, str_variable )
  END SUBROUTINE message2Strings3

  SUBROUTINE message1String2 ( MsgPriority, SMainMsg )
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: SMainMsg
    CALL message1String ( MsgPriority, no_tag, SMainMsg )  
  END SUBROUTINE message1String2

  SUBROUTINE message1String3 ( SMainMsg )
    CHARACTER(*), INTENT(IN) :: SMainMsg
    CALL message1String ( LS, no_tag, SMainMsg )  
  END SUBROUTINE message1String3

  SUBROUTINE messageInteger_isi2 ( MsgPriority, text1, int1, text2, int2)
    TYPE (MsgPriorities), INTENT(IN) :: MsgPriority
    CHARACTER(*), INTENT(IN) :: text1, text2
    INTEGER(IKIND), INTENT(IN) :: int1, int2
    CALL messageInteger_isi ( MsgPriority, no_tag, text1, int1, text2, int2)
  END SUBROUTINE messageInteger_isi2

END MODULE message_mod
