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
!  along with Felix.  If not, see <http://www.gnu.org/licenses/>.
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
!>  message( Smsg_priority, msg_tag, text_to_print, variable_to_print )
!>
!>    Smsg_priority     
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
!>    msg_tag       
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
!>            a format depending upon its TYPE and dimension. E.g. a REAL matrix
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
!>        LXL messages are usually accessed via a msg_tag instead.
!>
!>    The tens component of the INTEGER relates to the msg_tag. 
!>        
!>        0   - no_tag  : no specific msg_tag's are switched on, only priority is considered
!>        30  - dbg3    : corresponds to IWriteFLAG = 3 IN THE OLD SYSTEM, many print
!>                        statements in the code previously only printed when IWriteFLAG = 3
!>
!>        For the full list of msg_tag's and their corresponding IWriteFLAG ten's value,
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
!?? standardise format top level, spaces and INTEGERs...
!?? consider all variables as 2 dimensional matrices?
!?? consider could make dbg initialise select case more concise with array of derived
!?? rearrange to have terminal_error modules, msg_output_group examples, 
!?? rename to Smsg_priority, update Smsg_priority for clarity high vs. low priority
!?? priority for top level refinement details
!?? adding new msg_tags
!?? generic -> optional arguments not implimented...


!---------------------------------------------------------------------------------------------

  use l_alert_mod             ! grants access to felix's main error handling
  use MyNumbers               ! necesary for IKIND, RKIND etc.
  use MyMPI, ONLY : my_rank   ! necesary for my_rank - used to print on core 0 only

  !?? JR this grants MyNumbers, USE MyNumbers throughout code is for clariry...

!---------------------------------------------------------------------------------------------

! HUGE INTERFACE to handle various variables TYPEs and optional arguments

  interface message

  !-------------------------------------------------------------------------------------------
  ! main REAL/complex/INTEGER vector printing - used by matrix and scalar printing
    MODULE procedure message_rvector        
    MODULE procedure message_cvector       
    MODULE procedure message_ivector
  !-------------------------------------------------------------------------------------------                                        
    MODULE procedure message_string         ! text_to_print and string_variable
    MODULE procedure message_only           ! text_to_print only
    MODULE procedure message_LOGICAL        ! text_to_print and LOGICAL_variable
  !-------------------------------------------------------------------------------------------
  ! special versions of message
    MODULE procedure message_alongside_ir   ! prints INT & REAL matrix alongside
    MODULE procedure message_alongside_ir2  ! prints INT matrix & REAL vector alongside
    MODULE procedure message_alongside_ic   ! prints INT & complex matrix alongside
    MODULE procedure message_INTEGER_isi    ! prints text, INTEGER, more text, INTEGER 
  !-------------------------------------------------------------------------------------------
  ! interfaces used for scalars, matrices and optional arguments
    MODULE procedure message_rmatrix        
    MODULE procedure message_cmatrix        ! suffix :
    MODULE procedure message_imatrix        ! 2 = without msg_tag argument
    MODULE procedure message_real           ! 3 = without msg_tag & Smsg_priority argument  
    MODULE procedure message_real2          
    MODULE procedure message_real3          
    MODULE procedure message_rvector2 
    MODULE procedure message_rvector3
    MODULE procedure message_rmatrix2
    MODULE procedure message_rmatrix3
    MODULE procedure message_complex
    MODULE procedure message_complex2
    MODULE procedure message_complex3
    MODULE procedure message_cvector2
    MODULE procedure message_cvector3
    MODULE procedure message_cmatrix2
    MODULE procedure message_cmatrix3
    MODULE procedure message_INTEGER
    MODULE procedure message_INTEGER2
    MODULE procedure message_INTEGER3
    MODULE procedure message_ivector2
    MODULE procedure message_ivector3
    MODULE procedure message_imatrix2
    MODULE procedure message_imatrix3
    MODULE procedure message_string2
    MODULE procedure message_string3
    MODULE procedure message_only2
    MODULE procedure message_only3
    MODULE procedure message_INTEGER_isi2

  END interface message

  !--------------------------------------------------------------------

  TYPE msg_priorities
    LOGICAL :: state              ! set to .true. or .false. at run-time by IWriteFLAG
    CHARACTER(50) :: initial_msg  ! used for priority indenting - tiered structure
  END TYPE

  TYPE msg_tags
    LOGICAL :: state              ! set to .true. or .false. at run-time by IWriteFLAG
    INTEGER(IKIND) :: id_number   ! comapred with IWriteFLAG to set state to .true. or .false
  END TYPE

  !--------------------------------------------------------------------

  TYPE(msg_priorities) :: LS, LM, LL, LXL

  TYPE(msg_tags) :: no_tag, dbg7, dbg3, dbg6, dbg14

  LOGICAL,private :: l_print_this_core ! used to print on all cores, usually = .false.
  
  CHARACTER(:), allocatable :: indent_spaces, spaces ! used for tiered-structure 

  INTEGER(IKIND) :: iclock_rate ! used with print_end_timer
  
CONTAINS

  !--------------------------------------------------------------------
  ! setup message (on felix start-up, before IWriteFLAG read-in)
  !--------------------------------------------------------------------
  
  SUBROUTINE init_message_LOGICALs  
 
    !--------------------------------------------------------------------
    ! msg_tags fixed ID numbers (to compare against IWriteFLAG)
    !--------------------------------------------------------------------

    ! to add another msg_tag, declare it above and add it to set_message_mod_mode below
    no_tag%id_number  = 0
    dbg7%id_number    = 70
    dbg3%id_number    = 30
    dbg6%id_number    = 60
    dbg14%id_number   = 90 !?? would use 140, but our INTEGERs cannot exceed 140

    !--------------------------------------------------------------------
    ! setup tiered-structure
    !--------------------------------------------------------------------

    ! e.g. WRITE(*,formatting) TRIM(Smsg_priority%initial_msg)//spaces, main_msg, rvector

    indent_spaces   = "" 

    LS%initial_msg  = indent_spaces
    LM%initial_msg  = indent_spaces//"   "
    LL%initial_msg  = indent_spaces//"     "
    LXL%initial_msg = indent_spaces//"       "

    spaces = " " 
    ! spaces are the global spaces after intial_msg

    !--------------------------------------------------------------------
    ! intially only set LS = .true. before IWriteFLAG read-in
    !--------------------------------------------------------------------

    LS%state = .true.; 
    LM%state = .false.;   LL%state = .false.;   LXL%state = .false.

    no_tag%state=.false.; dbg7%state = .false.; dbg3%state = .false.;
    dbg6%state = .false.; dbg14%state = .false.

    l_print_this_core = .false.   

  END SUBROUTINE

  !--------------------------------------------------------------------
  ! read-in IWriteFLAG and setup message
  !--------------------------------------------------------------------

  SUBROUTINE set_message_mod_mode ( IWriteFLAG_in, ierr )
    
    INTEGER(IKIND), INTENT(out):: ierr
    INTEGER(IKIND), INTENT(IN) :: IWriteFLAG_in
    INTEGER(IKIND) :: prio                        ! acts as IWriteFLAG
    prio = IWriteFLAG_in    
  
    !--------------------------------------------------------------------
    ! Use IWriteFLAG tens component to turn-on a msg_tag
    !--------------------------------------------------------------------

    if ( prio > 10 ) then

      if ( (dbg7%id_number <= prio) .and. (prio  <= dbg7%id_number+10) ) then
        dbg7%state = .true.
        prio = prio - dbg7%id_number
      elseif ( (dbg3%id_number <= prio) .and. (prio  <= dbg3%id_number+10) ) then
        dbg3%state = .true.
        prio = prio - dbg3%id_number
      elseif ( (dbg6%id_number <= prio) .and. (prio  <= dbg6%id_number+10) ) then
        dbg6%state = .true.
        prio = prio - dbg6%id_number
      elseif ( (dbg14%id_number <= prio) .and. (prio  <= dbg14%id_number+10) ) then
        dbg14%state = .true.
        prio = prio - dbg14%id_number
      else  ! error (minor) - tens component not recognised
        !?? set IWriteFLAG to 2, just switch LS & LM on
        IErr=1
        IF(l_alert(IErr,"set_message_mod_mode",&
              "IWriteFLAG tens component not recognised check felix.inp")) RETURN
      END if

    END if
    ! prio is just a single digit INTEGER now

    ! switch lower priorities on depending upon IWriteFLAG digit compoent
    if (prio >= 2) LM%state = .true.;
    if (prio >= 4) LL%state = .true.;
    if (prio >= 8) LXL%state = .true.;

  END SUBROUTINE

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE allow_message_on_this_core
    l_print_this_core = .true.
    CALL message_INTEGER3("MESSAGE FROM THIS CORE AS WELL, rank =",my_rank)
  END SUBROUTINE  !?? not currently used, but may be useful

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !---------------------------------------------------------------------
  ! print_end_time
  !---------------------------------------------------------------------  

  ! compare start time to current and print time-passed
  SUBROUTINE print_end_time(Smsg_priority, Istart_time, STaskName)
    USE IPARA, ONLY : IClockRate
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
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
    CALL message(Smsg_priority,SPrintString)

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
  ! main REAL/complex/INTEGER vector printing - used by matrix and scalar printing
  !---------------------------------------------------------------------

  SUBROUTINE message_rvector ( Smsg_priority, msg_tag, main_msg, rvector ) 
   
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    REAL(RKIND),dimension(:),INTENT(IN) :: rvector ! REAL vector
    CHARACTER(50) :: formatting
    
    ! check priority then print REAL vector
    if ( ( my_rank==0 .or. l_print_this_core ) & 
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      if (size(rvector) >= 2) then ! vector, so surround with brackets
        WRITE(formatting,'(a,i3.3,a)') '(a,a,"("',size(rvector),'(1x,sp,ES10.3)")")'
      else  ! scalar, so bracketless
        WRITE(formatting,'(a,i3.3,a)') '(a,a,',size(rvector),'(1x,sp,ES10.3))'
      END if
      WRITE(*,formatting) TRIM(Smsg_priority%initial_msg)//spaces, main_msg, rvector
    END if

  END SUBROUTINE message_rvector

  SUBROUTINE message_cvector ( Smsg_priority, msg_tag, main_msg, cvector )    

    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    complex(CKIND),dimension(:),INTENT(IN) :: cvector ! complex vector
    CHARACTER(50) :: formatting

    ! check priority then print complex vector
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      WRITE(formatting,'(a,i3.3,a)') '(a,a,',size(cvector),'(1x,"(",sp,ES8.1," ",ES8.1,"i)"))'
      WRITE(*,formatting) TRIM(Smsg_priority%initial_msg)//spaces, main_msg, cvector
    END if

  END SUBROUTINE message_cvector

  SUBROUTINE message_ivector ( Smsg_priority, msg_tag, main_msg, ivector )  
  
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    INTEGER(IKIND),dimension(:),INTENT(IN) :: ivector ! INTEGER vector
    CHARACTER(50) :: formatting

    ! check priority then print INTEGER vector
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      WRITE(formatting,'(a,i3.3,a)') '(a,a,',size(ivector),'(1x,i4.3))'
      WRITE(*,formatting) TRIM(Smsg_priority%initial_msg)//spaces, main_msg, ivector
    END if

  END SUBROUTINE message_ivector

  !-------------------------------------------------------------------------------------------



  SUBROUTINE message_string ( Smsg_priority, msg_tag, main_msg, str_variable )

    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg, str_variable

    ! check priority then print main_msg & String
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      WRITE(*,'(a,a,a,a,a)') TRIM(Smsg_priority%initial_msg)//spaces, main_msg,"'",str_variable,"'"
    END if
  
  END SUBROUTINE message_string

  SUBROUTINE message_only ( Smsg_priority, msg_tag, main_msg )

    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg

    ! check priority then print main_msg
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      WRITE(*,'(a,a)') TRIM(Smsg_priority%initial_msg)//spaces, main_msg
    END if
  
  END SUBROUTINE message_only

  ! For printing LOGICALs ! currently no interfaces for optional arguments, matrices, vectors
  SUBROUTINE message_LOGICAL ( Smsg_priority, msg_tag, main_msg, LOGICAL_var )  
  
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    LOGICAL, INTENT(IN) :: LOGICAL_var ! LOGICAL variable

    ! check priority then print LOGICAL variable
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      WRITE(*,'(a,a,1x,l1)') TRIM(Smsg_priority%initial_msg)//spaces, main_msg, LOGICAL_var
    END if

  END SUBROUTINE message_LOGICAL



  !---------------------------------------------------------------------
  ! Special message varieties
  !---------------------------------------------------------------------
  
  ! vertically INT matrix beside REAL matrix
  SUBROUTINE message_alongside_ir ( Smsg_priority, msg_tag, &
               main_msg, imatrix1, rmatrix2 )
    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    INTEGER(IKIND),dimension(:,:),INTENT(IN) :: imatrix1 ! 1st matrix - INTEGER
    REAL(RKIND), dimension(:,:), INTENT(IN) :: rmatrix2 ! 2nd matrix - REAL
    CHARACTER(100) :: formatting
    INTEGER(IKIND) :: i

    ! check priority then print two matrices alongside
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      WRITE(*,'(a,a)') TRIM(Smsg_priority%initial_msg)//spaces, main_msg
      WRITE(formatting,'(a,i3.3,a,i3.3,a)') '(a,i3.3,',size(imatrix1,2),'(1x,i4.3),a,',&
            size(rmatrix2,2),'(1x,sp,ES10.3))'
      do i =1,size(imatrix1,1)
        WRITE(*,formatting) TRIM(Smsg_priority%initial_msg)//spaces,i, imatrix1(i,:)," |",rmatrix2(i,:)
      END do
    END if  
  
  END SUBROUTINE message_alongside_ir

  ! vertically INT matrix and REAL vector
  SUBROUTINE message_alongside_ir2 ( Smsg_priority, msg_tag, &
               main_msg, imatrix, rvector )

    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    INTEGER(IKIND),dimension(:,:),INTENT(IN) :: imatrix ! 1st matrix - INTEGER
    REAL(RKIND), dimension(:), INTENT(IN) :: rvector ! 2nd matrix - REAL
    CHARACTER(100) :: formatting
    INTEGER(IKIND) :: i

    CALL message_alongside_ir ( Smsg_priority, msg_tag, &
               main_msg, imatrix, RESHAPE(rvector,[size(rvector),1]) ) 
  
  END SUBROUTINE message_alongside_ir2

  ! vertically INT matrix beside complex matrix
  SUBROUTINE message_alongside_ic ( Smsg_priority, msg_tag, &
               main_msg, imatrix1, cmatrix2 )
    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    INTEGER(IKIND), dimension(:,:),INTENT(IN) :: imatrix1 ! 1st matrix - INTEGER
    complex(CKIND), dimension(:,:), INTENT(IN) :: cmatrix2 ! 2nd matrix - REAL
    CHARACTER(100) :: formatting
    INTEGER(IKIND) :: i

    ! check priority then print two matrices alongside
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      WRITE(*,'(a,a)') TRIM(Smsg_priority%initial_msg)//spaces, main_msg
      WRITE(formatting,'(a,i3.3,a,i3.3,a)') '(a,i3.3,',size(imatrix1,2),'(1x,i4.3),a,',&
            size(cmatrix2,2),'(1x,"(",sp,ES8.1," ",ES8.1,"i)"))'
      do i =1,size(imatrix1,1)
        WRITE(*,formatting) TRIM(Smsg_priority%initial_msg)//spaces,i, imatrix1(i,:)," |",cmatrix2(i,:)
      END do
    END if  
  
  END SUBROUTINE message_alongside_ic

  ! text, INTEGER, extra-text, extra-INTEGER
  SUBROUTINE message_INTEGER_isi ( Smsg_priority, msg_tag, text1, int1, text2, int2)

    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: text1, text2
    INTEGER(IKIND), INTENT(IN) :: int1, int2

    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      WRITE(*,'(a,a,i4.3,a,i4.3)') TRIM(Smsg_priority%initial_msg)//spaces, text1, int1, text2, int2
    END if
  
  END SUBROUTINE message_INTEGER_isi

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  ! INTERFACES BELOW
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Matrix messages using corresponding vector messages
  !---------------------------------------------------------------------

  SUBROUTINE message_rmatrix ( Smsg_priority, msg_tag, main_msg, rmatrix )
    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    REAL(RKIND),dimension(:,:),INTENT(IN) :: rmatrix ! REAL martrix
    INTEGER(IKIND) :: i

    ! check priority then print REAL matrix using vector print
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      CALL message_only ( Smsg_priority, msg_tag, main_msg )
      do i =1,size(rmatrix,1)
        CALL message_rvector ( Smsg_priority, msg_tag, "", rmatrix(i,:) ) 
      END do
    END if

  END SUBROUTINE message_rmatrix

  SUBROUTINE message_cmatrix ( Smsg_priority, msg_tag, main_msg, cmatrix )
    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    complex(CKIND),dimension(:,:),INTENT(IN) :: cmatrix ! complex martrix
    INTEGER(IKIND) :: i

    ! check priority then print complex matrix using vector print
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      CALL message_only ( Smsg_priority, msg_tag, main_msg )
      do i =1,size(cmatrix,1)
        CALL message_cvector ( Smsg_priority, msg_tag, "", cmatrix(i,:) ) 
      END do
    END if

  END SUBROUTINE message_cmatrix

  SUBROUTINE message_imatrix ( Smsg_priority, msg_tag, main_msg, imatrix )
    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    INTEGER(IKIND),dimension(:,:),INTENT(IN) :: imatrix ! INTEGER martrix
    INTEGER(IKIND) :: i

    ! check priority then print INTEGER matrix using vector print
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      CALL message_only ( Smsg_priority, msg_tag, main_msg )
      do i =1,size(imatrix,1)
        CALL message_ivector ( Smsg_priority, msg_tag, "", imatrix(i,:) ) 
      END do
    END if

  END SUBROUTINE message_imatrix

  !---------------------------------------------------------------------
  ! Variable, vector, matrix interfaces for optional (priority and tag) arguments
  !---------------------------------------------------------------------

  SUBROUTINE message_real ( Smsg_priority, msg_tag, main_msg, real_variable )
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    REAL(RKIND), INTENT(IN) :: real_variable
    CALL message_rvector ( Smsg_priority, msg_tag, main_msg, (/ real_variable /) )  
  END SUBROUTINE message_real

  SUBROUTINE message_real2 ( Smsg_priority, main_msg, real_variable )
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    CHARACTER(*), INTENT(IN) :: main_msg
    REAL(RKIND), INTENT(IN) :: real_variable
    CALL message_rvector ( Smsg_priority, no_tag, main_msg, (/ real_variable /) )  
  END SUBROUTINE message_real2

  SUBROUTINE message_real3 ( main_msg, real_variable )
    CHARACTER(*), INTENT(IN) :: main_msg
    REAL(RKIND), INTENT(IN) :: real_variable
    CALL message_rvector ( LS, no_tag, main_msg, (/ real_variable /) )  
  END SUBROUTINE message_real3

  SUBROUTINE message_rvector2 ( Smsg_priority, main_msg, rvector )    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    CHARACTER(*), INTENT(IN) :: main_msg
    REAL(RKIND),dimension(:),INTENT(IN) :: rvector
    CALL message_rvector ( Smsg_priority, no_tag, main_msg, rvector )
  END SUBROUTINE message_rvector2

  SUBROUTINE message_rvector3 ( main_msg, rvector )    
    CHARACTER(*), INTENT(IN) :: main_msg
    REAL(RKIND),dimension(:),INTENT(IN) :: rvector
    CALL message_rvector ( LS, no_tag, main_msg, rvector )
  END SUBROUTINE message_rvector3

  SUBROUTINE message_rmatrix2 ( Smsg_priority, main_msg, rmatrix )    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    CHARACTER(*), INTENT(IN) :: main_msg
    REAL(RKIND),dimension(:,:),INTENT(IN) :: rmatrix
    CALL message_rmatrix ( Smsg_priority, no_tag, main_msg, rmatrix )
  END SUBROUTINE message_rmatrix2

  SUBROUTINE message_rmatrix3 ( main_msg, rmatrix )    
    CHARACTER(*), INTENT(IN) :: main_msg
    REAL(RKIND),dimension(:,:),INTENT(IN) :: rmatrix
    CALL message_rmatrix ( LS, no_tag, main_msg, rmatrix )
  END SUBROUTINE message_rmatrix3  

  SUBROUTINE message_complex ( Smsg_priority, msg_tag, main_msg, c_variable )    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    complex(CKIND),INTENT(IN) :: c_variable
    CALL message_cvector ( Smsg_priority, msg_tag, main_msg, (/ c_variable /) )
  END SUBROUTINE message_complex

  SUBROUTINE message_complex2 ( Smsg_priority, main_msg, c_variable )    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    CHARACTER(*), INTENT(IN) :: main_msg
    complex(CKIND),INTENT(IN) :: c_variable
    CALL message_cvector ( Smsg_priority, no_tag, main_msg, (/ c_variable /) )
  END SUBROUTINE message_complex2

  SUBROUTINE message_complex3 ( main_msg, c_variable ) 
    CHARACTER(*), INTENT(IN) :: main_msg
    complex(CKIND),INTENT(IN) :: c_variable
    CALL message_cvector ( LS, no_tag, main_msg, (/ c_variable /) )
  END SUBROUTINE message_complex3

  SUBROUTINE message_cvector2 ( Smsg_priority, main_msg, cvector )    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    CHARACTER(*), INTENT(IN) :: main_msg
    complex(CKIND), dimension(:), INTENT(IN) :: cvector
    CALL message_cvector ( Smsg_priority, no_tag, main_msg, cvector )
  END SUBROUTINE message_cvector2

  SUBROUTINE message_cvector3 ( main_msg, cvector )    
    CHARACTER(*), INTENT(IN) :: main_msg
    complex(CKIND), dimension(:), INTENT(IN) :: cvector
    CALL message_cvector( LS, no_tag, main_msg, cvector )
  END SUBROUTINE message_cvector3

  SUBROUTINE message_cmatrix2 ( Smsg_priority, main_msg, cmatrix )    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    CHARACTER(*), INTENT(IN) :: main_msg
    complex(CKIND), dimension(:,:), INTENT(IN) :: cmatrix
    CALL message_cmatrix ( Smsg_priority, no_tag, main_msg, cmatrix )
  END SUBROUTINE message_cmatrix2

  SUBROUTINE message_cmatrix3 ( main_msg, cmatrix )    
    CHARACTER(*), INTENT(IN) :: main_msg
    complex(CKIND), dimension(:,:), INTENT(IN) :: cmatrix
    CALL message_cmatrix ( LS, no_tag, main_msg, cmatrix )
  END SUBROUTINE message_cmatrix3

  SUBROUTINE message_INTEGER ( Smsg_priority, msg_tag, main_msg, int_variable )
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    TYPE (msg_tags), INTENT(IN) :: msg_tag
    CHARACTER(*), INTENT(IN) :: main_msg
    INTEGER(IKIND), INTENT(IN) :: int_variable
    CALL message_ivector ( Smsg_priority, msg_tag, main_msg, (/ int_variable /) )
  END SUBROUTINE message_INTEGER

  SUBROUTINE message_INTEGER2 ( Smsg_priority, main_msg, int_variable )
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    CHARACTER(*), INTENT(IN) :: main_msg
    INTEGER(IKIND), INTENT(IN) :: int_variable
    CALL message_ivector ( Smsg_priority, no_tag, main_msg, (/ int_variable /) )  
  END SUBROUTINE message_INTEGER2

  SUBROUTINE message_INTEGER3 ( main_msg, int_variable )
    CHARACTER(*), INTENT(IN) :: main_msg
    INTEGER(IKIND), INTENT(IN) :: int_variable
    CALL message_ivector ( LS, no_tag, main_msg, (/ int_variable /) )  
  END SUBROUTINE message_INTEGER3

  SUBROUTINE message_ivector2 ( Smsg_priority, main_msg, ivector )    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    CHARACTER(*), INTENT(IN) :: main_msg
    INTEGER(IKIND),dimension(:),INTENT(IN) :: ivector
    CALL message_ivector ( Smsg_priority, no_tag, main_msg, ivector )
  END SUBROUTINE message_ivector2

  SUBROUTINE message_ivector3 ( main_msg, ivector )    
    CHARACTER(*), INTENT(IN) :: main_msg
    INTEGER(IKIND),dimension(:),INTENT(IN) :: ivector
    CALL message_ivector ( LS, no_tag, main_msg, ivector )
  END SUBROUTINE message_ivector3

  SUBROUTINE message_imatrix2 ( Smsg_priority, main_msg, imatrix )    
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    CHARACTER(*), INTENT(IN) :: main_msg
    INTEGER(IKIND),dimension(:,:),INTENT(IN) :: imatrix
    CALL message_imatrix ( Smsg_priority, no_tag, main_msg, imatrix )
  END SUBROUTINE message_imatrix2

  SUBROUTINE message_imatrix3 ( main_msg, imatrix )    
    CHARACTER(*), INTENT(IN) :: main_msg
    INTEGER(IKIND),dimension(:,:),INTENT(IN) :: imatrix
    CALL message_imatrix ( LS, no_tag, main_msg, imatrix )
  END SUBROUTINE message_imatrix3

  SUBROUTINE message_string2 ( Smsg_priority, main_msg, str_variable )
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    CHARACTER(*), INTENT(IN) :: main_msg, str_variable
    CALL message_string ( Smsg_priority, no_tag, main_msg, str_variable )
  END SUBROUTINE message_string2

  SUBROUTINE message_string3 ( main_msg, str_variable )
    CHARACTER(*), INTENT(IN) :: main_msg, str_variable
    CALL message_string ( LS, no_tag, main_msg, str_variable )
  END SUBROUTINE message_string3

  SUBROUTINE message_only2 ( Smsg_priority, main_msg )
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    CHARACTER(*), INTENT(IN) :: main_msg
    CALL message_only ( Smsg_priority, no_tag, main_msg )  
  END SUBROUTINE message_only2

  SUBROUTINE message_only3 ( main_msg )
    CHARACTER(*), INTENT(IN) :: main_msg
    CALL message_only ( LS, no_tag, main_msg )  
  END SUBROUTINE message_only3

  SUBROUTINE message_INTEGER_isi2 ( Smsg_priority, text1, int1, text2, int2)
    TYPE (msg_priorities), INTENT(IN) :: Smsg_priority
    CHARACTER(*), INTENT(IN) :: text1, text2
    INTEGER(IKIND), INTENT(IN) :: int1, int2
    CALL message_INTEGER_isi ( Smsg_priority, no_tag, text1, int1, text2, int2)
  END SUBROUTINE message_INTEGER_isi2

END MODULE message_mod
