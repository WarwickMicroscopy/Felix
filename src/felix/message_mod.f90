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

module message_mod
!>--------------------------------------------------------------------------------------------
!>  major-authors: Jacob Richardson (2017)
!>--------------------------------------------------------------------------------------------
!>
!>  MODULE OVERVIEW:
!>
!>    This module directly contains everything for message(),
!>      which is a major subroutine used throughout felix. It replaces the FORTRAN print
!>      statement and provides formatted output chosen at run-time.
!>
!>    message_mod also accesses l_alert_mod:
!>      l_alert() is a simple but important function, used extensively for error handling.
!>
!>    message_mod also conatins print_end_time() used to neatly print the time elapsed.
!>      
!>
!>--------------------------------------------------------------------------------------------
!>
!>  message( ... ) FEATURES:
!>
!>    IWriteFLAG in felix.inp is the run-time input used to select the output mode of Felix.
!>    message() uses IWriteFLAG to decide whether a message should be printed.
!>
!>    - message() prints various different variable types with consistent formatting
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
!>            DERIVED TYPE ( msg_priority_type ) = LS, LM, LL or LXL
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
!>            DERIVED TYPE ( msg_tag_type )
!>              
!>            This allows for messages to be grouped together with tags. At run-time you 
!>            can select that only key messages and messages with a certain tag are printed.
!>
!>    text_to_print
!>     
!>            TYPE: string (1D character array)
!>
!>            This is some text to print to the terminal.
!>            Any variable values will be printed after this text. 
!>                        
!>    variable_to_print 
!>    (optional)
!> 
!>            TYPE: complex, integer, real, logical, string 
!>            DIMENSION: scalar, vector or matrix
!>
!>            The value of this variable will be printed to the terminal with
!>            a format depending upon its type and dimension. E.g. a real matrix
!>            will be printed on multiple lines in columns with scientific notation.                
!>                       
!>--------------------------------------------------------------------------------------------
!>
!>  Using IWriteFLAG to choose the terminal output mode:
!>
!>    In felix.inp, IWriteFLAG can be assigned different integer values to affect which 
!>    messages are printed. 
!>
!>    The digit of the integer affects the minimum priority of the message printed:
!>      
!>        0 = LS,   2 = LM,   4 = LL,   8 = LXL   
!>
!>        sidenote: 8 rarely used, it prints every message and creates a very verbose output.
!>        LXL messages are usually accessed via a msg_tag instead.
!>
!>    The tens component of the integer relates to the msg_tag. 
!>        
!>        0   - no_tag  : no specific msg_tag's are switched on, only priority is considered
!>        30  - dbg3    : corresponds to IWriteFLAG = 3 IN THE OLD SYSTEM, many print
!>                        statements in the code previously only printed when IWriteFLAG = 3
!>
!>        For the full list of msg_tag's and their corresponding IWriteFLAG ten's value,
!>        refer to init_message_logicals() below.
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
!?? currently using derived type to force particular priority/dbg use
!?? standardise format top level, spaces and integers...
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

! HUGE INTERFACE to handle various variables types and optional arguments

  interface message

  !-------------------------------------------------------------------------------------------
  ! main real/complex/integer vector printing - used by matrix and scalar printing
    module procedure message_rvector        
    module procedure message_cvector       
    module procedure message_ivector
  !-------------------------------------------------------------------------------------------                                        
    module procedure message_string         ! text_to_print and string_variable
    module procedure message_only           ! text_to_print only
    module procedure message_logical        ! text_to_print and logical_variable
  !-------------------------------------------------------------------------------------------
  ! special versions of message
    module procedure message_alongside_ir   ! prints int & real matrix alongside
    module procedure message_alongside_ir2  ! prints int matrix & real vector alongside
    module procedure message_alongside_ic   ! prints int & complex matrix alongside
    module procedure message_integer_isi    ! prints text, integer, more text, integer 
  !-------------------------------------------------------------------------------------------
  ! interfaces used for scalars, matrices and optional arguments
    module procedure message_rmatrix        
    module procedure message_cmatrix        ! suffix :
    module procedure message_imatrix        ! 2 = without msg_tag argument
    module procedure message_real           ! 3 = without msg_tag & Smsg_priority argument  
    module procedure message_real2          
    module procedure message_real3          
    module procedure message_rvector2 
    module procedure message_rvector3
    module procedure message_rmatrix2
    module procedure message_rmatrix3
    module procedure message_complex
    module procedure message_complex2
    module procedure message_complex3
    module procedure message_cvector2
    module procedure message_cvector3
    module procedure message_cmatrix2
    module procedure message_cmatrix3
    module procedure message_integer
    module procedure message_integer2
    module procedure message_integer3
    module procedure message_ivector2
    module procedure message_ivector3
    module procedure message_imatrix2
    module procedure message_imatrix3
    module procedure message_string2
    module procedure message_string3
    module procedure message_only2
    module procedure message_only3
    module procedure message_integer_isi2

  end interface message

  !--------------------------------------------------------------------

  type msg_priorities
    logical :: state              ! set to .true. or .false. at run-time by IWriteFLAG
    character(50) :: initial_msg  ! used for priority indenting - tiered structure
  end type

  type msg_tags
    logical :: state              ! set to .true. or .false. at run-time by IWriteFLAG
    integer(IKIND) :: id_number   ! comapred with IWriteFLAG to set state to .true. or .false
  end type

  !--------------------------------------------------------------------

  type(msg_priorities) :: LS, LM, LL, LXL

  type(msg_tags) :: no_tag, dbg7, dbg3, dbg6, dbg14

  logical,private :: l_print_this_core ! used to print on all cores, usually = .false.
  
  character(:), allocatable :: indent_spaces, spaces ! used for tiered-structure 

  integer(IKIND) :: iclock_rate ! used with print_end_timer()
  
contains

  !--------------------------------------------------------------------
  ! setup message() (on felix start-up, before IWriteFLAG read-in)
  !--------------------------------------------------------------------
  
  subroutine init_message_logicals()  
 
    !--------------------------------------------------------------------
    ! msg_tags fixed ID numbers (to compare against IWriteFLAG)
    !--------------------------------------------------------------------

    ! to add another msg_tag, declare it above and add it to set_message_mod_mode() below
    no_tag%id_number  = 0
    dbg7%id_number    = 70
    dbg3%id_number    = 30
    dbg6%id_number    = 60
    dbg14%id_number   = 90 !?? would use 140, but our integers cannot exceed 140

    !--------------------------------------------------------------------
    ! setup tiered-structure
    !--------------------------------------------------------------------

    ! e.g. write(*,formatting) trim(Smsg_priority%initial_msg)//spaces, main_msg, rvector

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

  end subroutine

  !--------------------------------------------------------------------
  ! read-in IWriteFLAG and setup message()
  !--------------------------------------------------------------------

  subroutine set_message_mod_mode ( IWriteFLAG_in, ierr )
    
    integer(IKIND), intent(out):: ierr
    integer(IKIND), intent(in) :: IWriteFLAG_in
    integer(IKIND) :: prio                        ! acts as IWriteFLAG
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
      end if

    end if
    ! prio is just a single digit integer now

    ! switch lower priorities on depending upon IWriteFLAG digit compoent
    if (prio >= 2) LM%state = .true.;
    if (prio >= 4) LL%state = .true.;
    if (prio >= 8) LXL%state = .true.;

  end subroutine

  subroutine allow_message_on_this_core()
    l_print_this_core = .true.
    CALL message_integer3("MESSAGE FROM THIS CORE AS WELL, rank =",my_rank)
  end subroutine  !?? not currently used, but may be useful

  !---------------------------------------------------------------------
  ! print_end_time()
  !---------------------------------------------------------------------  

  ! compare start time to current and print time-passed
  SUBROUTINE print_end_time( msg_priority, istart_time, completed_task_name )
    type (msg_priorities), intent(in) :: msg_priority
    character(*), intent(in) :: completed_task_name
    integer(IKIND), intent(in) :: istart_time
    integer(IKIND) :: ihours,iminutes,iseconds,icurrent_time
    real(RKIND) :: duration
    character(100) :: string

    call system_clock(icurrent_time)
    ! converts ticks from system clock into seconds
    duration = real(icurrent_time-istart_time)/real(iclock_rate)
    ihours = floor(duration/3600.0d0)
    iminutes = floor(mod(duration,3600.0d0)/60.0d0)
    iseconds = int(mod(duration,3600.0d0)-iminutes*60)
    write(string,fmt='(a,1x,a,i3,a5,i2,a6,i2,a4)') completed_task_name,&
          "completed in ",ihours," hrs ",iminutes," mins ",iseconds," sec"
    call message_only2(msg_priority,trim(string))
    !?? currently message can't print the combined 3 integers and strings directly
  END SUBROUTINE print_end_time

  !   EXAMPLE USAGE - use intrinisic system_clock(), 'IStartTime2' local variable 
  !
  !   CALL SYSTEM_CLOCK( IStartTime2 )
  !   CALL Absorption (IErr)
  !   IF(l_alert(IErr,"felixrefine","Absorption()")) CALL abort()
  !   CALL print_end_time( LM, IStartTime2, "Absorption" )
  !
  !   EXAMPLE TERMINAL OUTPUT:
  !   @ ---- Absorption completed in   0 hrs  0 mins  2 sec

  !---------------------------------------------------------------------
  ! main real/complex/integer vector printing - used by matrix and scalar printing
  !---------------------------------------------------------------------

  subroutine message_rvector ( Smsg_priority, msg_tag, main_msg, rvector ) 
   
    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    real(RKIND),dimension(:),intent(in) :: rvector ! real vector
    character(50) :: formatting
    
    ! check priority then print real vector
    if ( ( my_rank==0 .or. l_print_this_core ) & 
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      if (size(rvector) >= 2) then ! vector, so surround with brackets
        write(formatting,'(a,i3.3,a)') '(a,a,"("',size(rvector),'(1x,sp,ES10.3)")")'
      else  ! scalar, so bracketless
        write(formatting,'(a,i3.3,a)') '(a,a,',size(rvector),'(1x,sp,ES10.3))'
      end if
      write(*,formatting) trim(Smsg_priority%initial_msg)//spaces, main_msg, rvector
    end if

  end subroutine message_rvector

  subroutine message_cvector ( Smsg_priority, msg_tag, main_msg, cvector )    

    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    complex(CKIND),dimension(:),intent(in) :: cvector ! complex vector
    character(50) :: formatting

    ! check priority then print complex vector
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      write(formatting,'(a,i3.3,a)') '(a,a,',size(cvector),'(1x,"(",sp,ES8.1," ",ES8.1,"i)"))'
      write(*,formatting) trim(Smsg_priority%initial_msg)//spaces, main_msg, cvector
    end if

  end subroutine message_cvector

  subroutine message_ivector ( Smsg_priority, msg_tag, main_msg, ivector )  
  
    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:),intent(in) :: ivector ! integer vector
    character(50) :: formatting

    ! check priority then print integer vector
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      write(formatting,'(a,i3.3,a)') '(a,a,',size(ivector),'(1x,i4.3))'
      write(*,formatting) trim(Smsg_priority%initial_msg)//spaces, main_msg, ivector
    end if

  end subroutine message_ivector

  !-------------------------------------------------------------------------------------------



  subroutine message_string ( Smsg_priority, msg_tag, main_msg, str_variable )

    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg, str_variable

    ! check priority then print main_msg & string
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      write(*,'(a,a,a,a,a)') trim(Smsg_priority%initial_msg)//spaces, main_msg,"'",str_variable,"'"
    end if
  
  end subroutine message_string

  subroutine message_only ( Smsg_priority, msg_tag, main_msg )

    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg

    ! check priority then print main_msg
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      write(*,'(a,a)') trim(Smsg_priority%initial_msg)//spaces, main_msg
    end if
  
  end subroutine message_only

  ! For printing logicals ! currently no interfaces for optional arguments, matrices, vectors
  subroutine message_logical ( Smsg_priority, msg_tag, main_msg, logical_var )  
  
    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    logical, intent(in) :: logical_var ! logical variable

    ! check priority then print logical variable
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      write(*,'(a,a,1x,l1)') trim(Smsg_priority%initial_msg)//spaces, main_msg, logical_var
    end if

  end subroutine message_logical



  !---------------------------------------------------------------------
  ! Special message varieties
  !---------------------------------------------------------------------
  
  ! vertically int matrix beside real matrix
  subroutine message_alongside_ir ( Smsg_priority, msg_tag, &
               main_msg, imatrix1, rmatrix2 )
    
    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:,:),intent(in) :: imatrix1 ! 1st matrix - integer
    real(RKIND), dimension(:,:), intent(in) :: rmatrix2 ! 2nd matrix - real
    character(100) :: formatting
    integer(IKIND) :: i

    ! check priority then print two matrices alongside
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      write(*,'(a,a)') trim(Smsg_priority%initial_msg)//spaces, main_msg
      write(formatting,'(a,i3.3,a,i3.3,a)') '(a,i3.3,',size(imatrix1,2),'(1x,i4.3),a,',&
            size(rmatrix2,2),'(1x,sp,ES10.3))'
      do i =1,size(imatrix1,1)
        write(*,formatting) trim(Smsg_priority%initial_msg)//spaces,i, imatrix1(i,:)," |",rmatrix2(i,:)
      end do
    end if  
  
  end subroutine message_alongside_ir

  ! vertically int matrix and real vector
  subroutine message_alongside_ir2 ( Smsg_priority, msg_tag, &
               main_msg, imatrix, rvector )

    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:,:),intent(in) :: imatrix ! 1st matrix - integer
    real(RKIND), dimension(:), intent(in) :: rvector ! 2nd matrix - real
    character(100) :: formatting
    integer(IKIND) :: i

    call message_alongside_ir ( Smsg_priority, msg_tag, &
               main_msg, imatrix, RESHAPE(rvector,[size(rvector),1]) ) 
  
  end subroutine message_alongside_ir2

  ! vertically int matrix beside complex matrix
  subroutine message_alongside_ic ( Smsg_priority, msg_tag, &
               main_msg, imatrix1, cmatrix2 )
    
    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    integer(IKIND), dimension(:,:),intent(in) :: imatrix1 ! 1st matrix - integer
    complex(CKIND), dimension(:,:), intent(in) :: cmatrix2 ! 2nd matrix - real
    character(100) :: formatting
    integer(IKIND) :: i

    ! check priority then print two matrices alongside
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      write(*,'(a,a)') trim(Smsg_priority%initial_msg)//spaces, main_msg
      write(formatting,'(a,i3.3,a,i3.3,a)') '(a,i3.3,',size(imatrix1,2),'(1x,i4.3),a,',&
            size(cmatrix2,2),'(1x,"(",sp,ES8.1," ",ES8.1,"i)"))'
      do i =1,size(imatrix1,1)
        write(*,formatting) trim(Smsg_priority%initial_msg)//spaces,i, imatrix1(i,:)," |",cmatrix2(i,:)
      end do
    end if  
  
  end subroutine message_alongside_ic

  ! text, integer, extra-text, extra-integer
  subroutine message_integer_isi ( Smsg_priority, msg_tag, text1, int1, text2, int2)

    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: text1, text2
    integer(IKIND), intent(in) :: int1, int2

    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      write(*,'(a,a,i4.3,a,i4.3)') trim(Smsg_priority%initial_msg)//spaces, text1, int1, text2, int2
    end if
  
  end subroutine message_integer_isi

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  ! INTERFACES BELOW
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Matrix messages using corresponding vector messages
  !---------------------------------------------------------------------

  subroutine message_rmatrix ( Smsg_priority, msg_tag, main_msg, rmatrix )
    
    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    real(RKIND),dimension(:,:),intent(in) :: rmatrix ! real martrix
    integer(IKIND) :: i

    ! check priority then print real matrix using vector print
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      call message_only ( Smsg_priority, msg_tag, main_msg )
      do i =1,size(rmatrix,1)
        call message_rvector ( Smsg_priority, msg_tag, "", rmatrix(i,:) ) 
      end do
    end if

  end subroutine message_rmatrix

  subroutine message_cmatrix ( Smsg_priority, msg_tag, main_msg, cmatrix )
    
    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    complex(CKIND),dimension(:,:),intent(in) :: cmatrix ! complex martrix
    integer(IKIND) :: i

    ! check priority then print complex matrix using vector print
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      call message_only ( Smsg_priority, msg_tag, main_msg )
      do i =1,size(cmatrix,1)
        call message_cvector ( Smsg_priority, msg_tag, "", cmatrix(i,:) ) 
      end do
    end if

  end subroutine message_cmatrix

  subroutine message_imatrix ( Smsg_priority, msg_tag, main_msg, imatrix )
    
    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:,:),intent(in) :: imatrix ! integer martrix
    integer(IKIND) :: i

    ! check priority then print integer matrix using vector print
    if ( ( my_rank==0 .or. l_print_this_core ) &
    .and. (Smsg_priority%state .or. msg_tag%state) ) then
      call message_only ( Smsg_priority, msg_tag, main_msg )
      do i =1,size(imatrix,1)
        call message_ivector ( Smsg_priority, msg_tag, "", imatrix(i,:) ) 
      end do
    end if

  end subroutine message_imatrix

  !---------------------------------------------------------------------
  ! Variable, vector, matrix interfaces for optional (priority and tag) arguments
  !---------------------------------------------------------------------

  subroutine message_real ( Smsg_priority, msg_tag, main_msg, real_variable )
    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    real(RKIND), intent(in) :: real_variable
    call message_rvector ( Smsg_priority, msg_tag, main_msg, (/ real_variable /) )  
  end subroutine message_real

  subroutine message_real2 ( Smsg_priority, main_msg, real_variable )
    type (msg_priorities), intent(in) :: Smsg_priority
    character(*), intent(in) :: main_msg
    real(RKIND), intent(in) :: real_variable
    call message_rvector ( Smsg_priority, no_tag, main_msg, (/ real_variable /) )  
  end subroutine message_real2

  subroutine message_real3 ( main_msg, real_variable )
    character(*), intent(in) :: main_msg
    real(RKIND), intent(in) :: real_variable
    call message_rvector ( LS, no_tag, main_msg, (/ real_variable /) )  
  end subroutine message_real3

  subroutine message_rvector2 ( Smsg_priority, main_msg, rvector )    
    type (msg_priorities), intent(in) :: Smsg_priority
    character(*), intent(in) :: main_msg
    real(RKIND),dimension(:),intent(in) :: rvector
    call message_rvector ( Smsg_priority, no_tag, main_msg, rvector )
  end subroutine message_rvector2

  subroutine message_rvector3 ( main_msg, rvector )    
    character(*), intent(in) :: main_msg
    real(RKIND),dimension(:),intent(in) :: rvector
    call message_rvector ( LS, no_tag, main_msg, rvector )
  end subroutine message_rvector3

  subroutine message_rmatrix2 ( Smsg_priority, main_msg, rmatrix )    
    type (msg_priorities), intent(in) :: Smsg_priority
    character(*), intent(in) :: main_msg
    real(RKIND),dimension(:,:),intent(in) :: rmatrix
    call message_rmatrix ( Smsg_priority, no_tag, main_msg, rmatrix )
  end subroutine message_rmatrix2

  subroutine message_rmatrix3 ( main_msg, rmatrix )    
    character(*), intent(in) :: main_msg
    real(RKIND),dimension(:,:),intent(in) :: rmatrix
    call message_rmatrix ( LS, no_tag, main_msg, rmatrix )
  end subroutine message_rmatrix3  

  subroutine message_complex ( Smsg_priority, msg_tag, main_msg, c_variable )    
    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    complex(CKIND),intent(in) :: c_variable
    call message_cvector ( Smsg_priority, msg_tag, main_msg, (/ c_variable /) )
  end subroutine message_complex

  subroutine message_complex2 ( Smsg_priority, main_msg, c_variable )    
    type (msg_priorities), intent(in) :: Smsg_priority
    character(*), intent(in) :: main_msg
    complex(CKIND),intent(in) :: c_variable
    call message_cvector ( Smsg_priority, no_tag, main_msg, (/ c_variable /) )
  end subroutine message_complex2

  subroutine message_complex3 ( main_msg, c_variable ) 
    character(*), intent(in) :: main_msg
    complex(CKIND),intent(in) :: c_variable
    call message_cvector ( LS, no_tag, main_msg, (/ c_variable /) )
  end subroutine message_complex3

  subroutine message_cvector2 ( Smsg_priority, main_msg, cvector )    
    type (msg_priorities), intent(in) :: Smsg_priority
    character(*), intent(in) :: main_msg
    complex(CKIND), dimension(:), intent(in) :: cvector
    call message_cvector ( Smsg_priority, no_tag, main_msg, cvector )
  end subroutine message_cvector2

  subroutine message_cvector3 ( main_msg, cvector )    
    character(*), intent(in) :: main_msg
    complex(CKIND), dimension(:), intent(in) :: cvector
    call message_cvector( LS, no_tag, main_msg, cvector )
  end subroutine message_cvector3

  subroutine message_cmatrix2 ( Smsg_priority, main_msg, cmatrix )    
    type (msg_priorities), intent(in) :: Smsg_priority
    character(*), intent(in) :: main_msg
    complex(CKIND), dimension(:,:), intent(in) :: cmatrix
    call message_cmatrix ( Smsg_priority, no_tag, main_msg, cmatrix )
  end subroutine message_cmatrix2

  subroutine message_cmatrix3 ( main_msg, cmatrix )    
    character(*), intent(in) :: main_msg
    complex(CKIND), dimension(:,:), intent(in) :: cmatrix
    call message_cmatrix ( LS, no_tag, main_msg, cmatrix )
  end subroutine message_cmatrix3

  subroutine message_integer ( Smsg_priority, msg_tag, main_msg, int_variable )
    type (msg_priorities), intent(in) :: Smsg_priority
    type (msg_tags), intent(in) :: msg_tag
    character(*), intent(in) :: main_msg
    integer(IKIND), intent(in) :: int_variable
    call message_ivector ( Smsg_priority, msg_tag, main_msg, (/ int_variable /) )
  end subroutine message_integer

  subroutine message_integer2 ( Smsg_priority, main_msg, int_variable )
    type (msg_priorities), intent(in) :: Smsg_priority
    character(*), intent(in) :: main_msg
    integer(IKIND), intent(in) :: int_variable
    call message_ivector ( Smsg_priority, no_tag, main_msg, (/ int_variable /) )  
  end subroutine message_integer2

  subroutine message_integer3 ( main_msg, int_variable )
    character(*), intent(in) :: main_msg
    integer(IKIND), intent(in) :: int_variable
    call message_ivector ( LS, no_tag, main_msg, (/ int_variable /) )  
  end subroutine message_integer3

  subroutine message_ivector2 ( Smsg_priority, main_msg, ivector )    
    type (msg_priorities), intent(in) :: Smsg_priority
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:),intent(in) :: ivector
    call message_ivector ( Smsg_priority, no_tag, main_msg, ivector )
  end subroutine message_ivector2

  subroutine message_ivector3 ( main_msg, ivector )    
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:),intent(in) :: ivector
    call message_ivector ( LS, no_tag, main_msg, ivector )
  end subroutine message_ivector3

  subroutine message_imatrix2 ( Smsg_priority, main_msg, imatrix )    
    type (msg_priorities), intent(in) :: Smsg_priority
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:,:),intent(in) :: imatrix
    call message_imatrix ( Smsg_priority, no_tag, main_msg, imatrix )
  end subroutine message_imatrix2

  subroutine message_imatrix3 ( main_msg, imatrix )    
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:,:),intent(in) :: imatrix
    call message_imatrix ( LS, no_tag, main_msg, imatrix )
  end subroutine message_imatrix3

  subroutine message_string2 ( Smsg_priority, main_msg, str_variable )
    type (msg_priorities), intent(in) :: Smsg_priority
    character(*), intent(in) :: main_msg, str_variable
    call message_string ( Smsg_priority, no_tag, main_msg, str_variable )
  end subroutine message_string2

  subroutine message_string3 ( main_msg, str_variable )
    character(*), intent(in) :: main_msg, str_variable
    call message_string ( LS, no_tag, main_msg, str_variable )
  end subroutine message_string3

  subroutine message_only2 ( Smsg_priority, main_msg )
    type (msg_priorities), intent(in) :: Smsg_priority
    character(*), intent(in) :: main_msg
    call message_only ( Smsg_priority, no_tag, main_msg )  
  end subroutine message_only2

  subroutine message_only3 ( main_msg )
    character(*), intent(in) :: main_msg
    call message_only ( LS, no_tag, main_msg )  
  end subroutine message_only3

  subroutine message_integer_isi2 ( Smsg_priority, text1, int1, text2, int2)
    type (msg_priorities), intent(in) :: Smsg_priority
    character(*), intent(in) :: text1, text2
    integer(IKIND), intent(in) :: int1, int2
    call message_integer_isi ( Smsg_priority, no_tag, text1, int1, text2, int2)
  end subroutine message_integer_isi2

end module message_mod
