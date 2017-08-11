!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>
!! Procedure-description: Printing to the terminal via this command allows
!! different messages to be printed depeding upon the debug mode selected at 
!! run-time. It is a generic subroutine and can handle a variety of variable
!! types as input.
!!
!! Major-Authors: Jacob Richardson (2017)
!!
module message_mod

  !A few examples below of accepted inputs

  !CALL message( priority_logical, debug_mode_logical, main_msg, int_variable )
  !CALL message( priority_logical, main_msg, int_variable )
  !CALL message( priority_logical, debug_mode_logical, main_msg)
  !CALL message( main_msg )
  !CALL message( priority_logical, debug_mode_logical, main_msg, imatrix1, rmatrix2 )

  ! a special case - if you pass two matrix variables, message assumes and prints
  ! them alongisde one another  

  !todo - implicit none & private
  !todo - add iwriteflag initialise
  !todo - consider debug modes for other corresponding input variables
  !todo - add current procedure state
  !todo - currently using derived type to force particular priority/dbg use
  !todo - possible standardise format...
  !todo - error message subroutine, however problematic goto...
  !rodo - add space after main_msg
  !todo - consider all variables as 2 dimensional matrices...
  !todo - consider could make dbg initialise select case more concise with array of derived

  use MyNumbers
  use IPARA !todo - IWriteFLAG?
  use SPARA !todo - CurrentProcedure?
  use IChannels
  use MPI
  use MYMPI

  interface message

    module procedure message_rvector
    module procedure message_cvector
    module procedure message_ivector
    module procedure message_string
    module procedure message_only

    module procedure message_logical

    module procedure message_alongside_ir ! int & real matrix alongside
    module procedure message_alongside_ic ! int & complex matrix alongside
    module procedure message_integer_isi ! logicals, main_msg, integer, string, integer

    module procedure message_rmatrix
    module procedure message_cmatrix
    module procedure message_imatrix

    module procedure message_real
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

  end interface message

  type priority_logicals
    logical :: state
    character(:),allocatable :: initial_msg
  end type

  type debug_mode_logicals
    logical :: state
    integer(IKIND) :: id_number
  end type

  ! priority logicals
  type(priority_logicals) :: LS, LM, LL, LXL
  ! specific debug modes logicals 
  type(debug_mode_logicals) :: dbg_default, dbg7, dbg3, dbg6, dbg14
  
  character(:), allocatable :: set_initial_msg

contains
  
  ! initialising immediately before reading-in felix.inp
  subroutine init_message_logicals()   
    
    dbg_default%id_number = 0
    dbg7%id_number = 70
    dbg3%id_number = 30
    dbg6%id_number = 60
    dbg14%id_number = 90


    LS%state = .true.; LM%state = .false.; LL%state = .false.; LXL%state = .false.
    dbg_default%state=.false.; dbg7%state = .false.; dbg3%state = .false.;
    dbg6%state = .false.; dbg14%state = .false.   
    set_initial_msg = "@ "  
    LS%initial_msg  = set_initial_msg
    LM%initial_msg  = "----"//set_initial_msg
    LL%initial_msg  = "--------"//set_initial_msg
    LXL%initial_msg = "-----------"//set_initial_msg

  end subroutine

  subroutine set_message_logicals( prio )
    
    integer(IKIND) :: prio

    if (prio < 10) then ! normal output mode (not debug)
      set_initial_msg = "@ "
    else

      set_initial_msg = "@ msg:"
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
      else ! error (minor) - incorrect input
        WRITE(*,*) "Error (minor) No debug mode matching that IWriteFLAG input from felix.inp"
        prio = 0
      end if

    end if
      
    LS%state = .true.;
    if (prio >= 2) LM%state = .true.;
    if (prio >= 4) LL%state = .true.;
    if (prio >= 8) LXL%state = .true.;
    
    ! indent messages depending upon priority / depth-of-detail
    LS%initial_msg  = set_initial_msg
    LM%initial_msg  = "----"//set_initial_msg
    LL%initial_msg  = "--------"//set_initial_msg
    LXL%initial_msg = "-----------"//set_initial_msg

  end subroutine

  !---------------------------------------------------------------------
  ! Key variable messaging - vectors used for single element and matrix messaging 
  subroutine message_rvector ( priority_logical, debug_mode_logical, main_msg, rvector ) 
   
    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg
    real(RKIND),dimension(:),intent(in) :: rvector ! real vector
    character(50) :: formatting
    
    ! check priority then print real vector
    if ( my_rank==0 .and. (priority_logical%state .or. debug_mode_logical%state) ) then
      if (size(rvector) >= 2) then ! vector, so surround with brackets
        write(formatting,'(a,i3.3,a)') '(1x,a,a,"("',size(rvector),'(1x,sp,ES10.3)")")'
      else  ! scalar, so bracketless
        write(formatting,'(a,i3.3,a)') '(1x,a,a,',size(rvector),'(1x,sp,ES10.3))'
      end if
      write(*,formatting) priority_logical%initial_msg, main_msg, rvector
    end if

  end subroutine message_rvector

  subroutine message_cvector ( priority_logical, debug_mode_logical, main_msg, cvector )    

    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg
    complex(CKIND),dimension(:),intent(in) :: cvector ! complex vector
    character(50) :: formatting

    ! check priority then print complex vector
    if ( my_rank==0 .and. (priority_logical%state .or. debug_mode_logical%state) ) then
      write(formatting,'(a,i3.3,a)') '(1x,a,a,',size(cvector),'(1x,"(",sp,ES8.1," ",ES8.1,"i)"))'
      write(*,formatting) priority_logical%initial_msg, main_msg, cvector
    end if

  end subroutine message_cvector

  subroutine message_ivector ( priority_logical, debug_mode_logical, main_msg, ivector )  
  
    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:),intent(in) :: ivector ! integer vector
    character(50) :: formatting

    ! check priority then print integer vector
    if ( my_rank==0 .and. (priority_logical%state .or. debug_mode_logical%state) ) then
      write(formatting,'(a,i3.3,a)') '(1x,a,a,',size(ivector),'(1x,i4.3))'
      write(*,formatting) priority_logical%initial_msg, main_msg, ivector
    end if

  end subroutine message_ivector

  subroutine message_string ( priority_logical, debug_mode_logical, main_msg, str_variable )

    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg, str_variable

    ! check priority then print main_msg & string
    if ( my_rank==0 .and. (priority_logical%state .or. debug_mode_logical%state) ) then
      write(*,'(1x,a,a,a)') priority_logical%initial_msg, main_msg, str_variable
    end if
  
  end subroutine message_string

  subroutine message_only ( priority_logical, debug_mode_logical, main_msg )

    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg

    ! check priority then print main_msg
    if ( my_rank==0 .and. (priority_logical%state .or. debug_mode_logical%state) ) then
      write(*,'(1x,a,a)') priority_logical%initial_msg, main_msg
    end if
  
  end subroutine message_only



  !---------------------------------------------------------------------
  ! For printing logicals, currently without interfaces, vectors, matrices
  subroutine message_logical ( priority_logical, debug_mode_logical, main_msg, logical_var )  
  
    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg
    logical, intent(in) :: logical_var ! logical variable

    ! check priority then print logical variable
    if ( my_rank==0 .and. (priority_logical%state .or. debug_mode_logical%state) ) then
      write(*,'(1x,a,a,1x,l1)') priority_logical%initial_msg, main_msg, logical_var
    end if

  end subroutine message_logical



  !---------------------------------------------------------------------
  ! Special cases like printing matrices alongside or multiple integers amongst text
  
  ! vertically int matrix beside real matrix
  subroutine message_alongside_ir ( priority_logical, debug_mode_logical, &
               main_msg, imatrix1, rmatrix2 )
    
    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:,:),intent(in) :: imatrix1 ! 1st matrix - integer
    real(RKIND), dimension(:,:), intent(in) :: rmatrix2 ! 2nd matrix - real
    character(100) :: formatting
    integer(IKIND) :: i

    ! check priority then print two matrices alongside
    if ( my_rank==0 .and. (priority_logical%state .or. debug_mode_logical%state) ) then
      write(*,'(1x,a,a)') priority_logical%initial_msg, main_msg
      write(formatting,'(a,i3.3,a,i3.3,a)') '(1x,a,i3.3,a,',size(imatrix1,2),'(1x,i4.3),a,',&
            size(rmatrix2,2),'(1x,sp,ES10.3))'
      do i =1,size(imatrix1,1)
        write(*,formatting) "@ ",i,":", imatrix1(i,:)," |",rmatrix2(i,:)
      end do
    end if  
  
  end subroutine message_alongside_ir

  ! vertically int matrix beside complex matrix
  subroutine message_alongside_ic ( priority_logical, debug_mode_logical, &
               main_msg, imatrix1, cmatrix2 )
    
    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg
    integer(IKIND), dimension(:,:),intent(in) :: imatrix1 ! 1st matrix - integer
    complex(CKIND), dimension(:,:), intent(in) :: cmatrix2 ! 2nd matrix - real
    character(100) :: formatting
    integer(IKIND) :: i

    ! check priority then print two matrices alongside
    if ( my_rank==0 .and. (priority_logical%state .or. debug_mode_logical%state) ) then
      write(*,'(1x,a,a)') priority_logical%initial_msg, main_msg
      write(formatting,'(a,i3.3,a,i3.3,a)') '(1x,a,i3.3,a,',size(imatrix1,2),'(1x,i4.3),a,',&
            size(cmatrix2,2),'(1x,"(",sp,ES8.1," ",ES8.1,"i)"))'
      do i =1,size(imatrix1,1)
        write(*,formatting) "@ ",i,":", imatrix1(i,:)," |",cmatrix2(i,:)
      end do
    end if  
  
  end subroutine message_alongside_ic

  ! text, integer, extra-text, extra-integer
  subroutine message_integer_isi ( priority_logical, debug_mode_logical, text1, int1, text2, int2)

    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: text1, text2
    integer(IKIND), intent(in) :: int1, int2

    if ( my_rank==0 .and. (priority_logical%state .or. debug_mode_logical%state) ) then
      write(*,'(1x,a,a,i4.3,a,i4.3)') priority_logical%initial_msg, text1, int1, text2, int2
    end if
  
  end subroutine message_integer_isi


  !---------------------------------------------------------------------
  ! Matrix messages using corresponding vector messages
  subroutine message_rmatrix ( priority_logical, debug_mode_logical, main_msg, rmatrix )
    
    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg
    real(RKIND),dimension(:,:),intent(in) :: rmatrix ! real martrix
    integer(IKIND) :: i

    ! check priority then print real matrix using vector print
    if ( my_rank==0 .and. (priority_logical%state .or. debug_mode_logical%state) ) then
      call message_only ( priority_logical, debug_mode_logical, main_msg )
      do i =1,size(rmatrix,1)
        call message_rvector ( priority_logical, debug_mode_logical, "", rmatrix(i,:) ) 
      end do
    end if

  end subroutine message_rmatrix

  subroutine message_cmatrix ( priority_logical, debug_mode_logical, main_msg, cmatrix )
    
    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg
    complex(CKIND),dimension(:,:),intent(in) :: cmatrix ! complex martrix
    integer(IKIND) :: i

    ! check priority then print complex matrix using vector print
    if ( my_rank==0 .and. (priority_logical%state .or. debug_mode_logical%state) ) then
      call message_only ( priority_logical, debug_mode_logical, main_msg )
      do i =1,size(cmatrix,1)
        call message_cvector ( priority_logical, debug_mode_logical, "", cmatrix(i,:) ) 
      end do
    end if

  end subroutine message_cmatrix

  subroutine message_imatrix ( priority_logical, debug_mode_logical, main_msg, imatrix )
    
    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:,:),intent(in) :: imatrix ! integer martrix
    integer(IKIND) :: i

    ! check priority then print integer matrix using vector print
    if ( my_rank==0 .and. (priority_logical%state .or. debug_mode_logical%state) ) then
      call message_only ( priority_logical, debug_mode_logical, main_msg )
      do i =1,size(imatrix,1)
        call message_ivector ( priority_logical, debug_mode_logical, "", imatrix(i,:) ) 
      end do
    end if

  end subroutine message_imatrix


  !---------------------------------------------------------------------
  ! Variable, vector, matrix interfaces for optional logical arguments
  subroutine message_real ( priority_logical, debug_mode_logical, main_msg, real_variable )
    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg
    real(RKIND), intent(in) :: real_variable
    call message_rvector ( priority_logical, debug_mode_logical, main_msg, (/ real_variable /) )  
  end subroutine message_real

  subroutine message_real2 ( priority_logical, main_msg, real_variable )
    type (priority_logicals), intent(in) :: priority_logical
    character(*), intent(in) :: main_msg
    real(RKIND), intent(in) :: real_variable
    call message_rvector ( priority_logical, dbg_default, main_msg, (/ real_variable /) )  
  end subroutine message_real2

  subroutine message_real3 ( main_msg, real_variable )
    character(*), intent(in) :: main_msg
    real(RKIND), intent(in) :: real_variable
    call message_rvector ( LS, dbg_default, main_msg, (/ real_variable /) )  
  end subroutine message_real3

  subroutine message_rvector2 ( priority_logical, main_msg, rvector )    
    type (priority_logicals), intent(in) :: priority_logical
    character(*), intent(in) :: main_msg
    real(RKIND),dimension(:),intent(in) :: rvector
    call message_rvector ( priority_logical, dbg_default, main_msg, rvector )
  end subroutine message_rvector2

  subroutine message_rvector3 ( main_msg, rvector )    
    character(*), intent(in) :: main_msg
    real(RKIND),dimension(:),intent(in) :: rvector
    call message_rvector ( LS, dbg_default, main_msg, rvector )
  end subroutine message_rvector3

  subroutine message_rmatrix2 ( priority_logical, main_msg, rmatrix )    
    type (priority_logicals), intent(in) :: priority_logical
    character(*), intent(in) :: main_msg
    real(RKIND),dimension(:,:),intent(in) :: rmatrix
    call message_rmatrix ( priority_logical, dbg_default, main_msg, rmatrix )
  end subroutine message_rmatrix2

  subroutine message_rmatrix3 ( main_msg, rmatrix )    
    character(*), intent(in) :: main_msg
    real(RKIND),dimension(:,:),intent(in) :: rmatrix
    call message_rmatrix ( LS, dbg_default, main_msg, rmatrix )
  end subroutine message_rmatrix3  

  subroutine message_complex ( priority_logical, debug_mode_logical, main_msg, c_variable )    
    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg
    complex(CKIND),intent(in) :: c_variable
    call message_cvector ( priority_logical, debug_mode_logical, main_msg, (/ c_variable /) )
  end subroutine message_complex

  subroutine message_complex2 ( priority_logical, main_msg, c_variable )    
    type (priority_logicals), intent(in) :: priority_logical
    character(*), intent(in) :: main_msg
    complex(CKIND),intent(in) :: c_variable
    call message_cvector ( priority_logical, dbg_default, main_msg, (/ c_variable /) )
  end subroutine message_complex2

  subroutine message_complex3 ( main_msg, c_variable ) 
    character(*), intent(in) :: main_msg
    complex(CKIND),intent(in) :: c_variable
    call message_cvector ( LS, dbg_default, main_msg, (/ c_variable /) )
  end subroutine message_complex3

  subroutine message_cvector2 ( priority_logical, main_msg, cvector )    
    type (priority_logicals), intent(in) :: priority_logical
    character(*), intent(in) :: main_msg
    complex(CKIND), dimension(:), intent(in) :: cvector
    call message_cvector ( priority_logical, dbg_default, main_msg, cvector )
  end subroutine message_cvector2

  subroutine message_cvector3 ( main_msg, cvector )    
    character(*), intent(in) :: main_msg
    complex(CKIND), dimension(:), intent(in) :: cvector
    call message_cvector( LS, dbg_default, main_msg, cvector )
  end subroutine message_cvector3

  subroutine message_cmatrix2 ( priority_logical, main_msg, cmatrix )    
    type (priority_logicals), intent(in) :: priority_logical
    character(*), intent(in) :: main_msg
    complex(CKIND), dimension(:,:), intent(in) :: cmatrix
    call message_cmatrix ( priority_logical, dbg_default, main_msg, cmatrix )
  end subroutine message_cmatrix2

  subroutine message_cmatrix3 ( main_msg, cmatrix )    
    character(*), intent(in) :: main_msg
    complex(CKIND), dimension(:,:), intent(in) :: cmatrix
    call message_cmatrix ( LS, dbg_default, main_msg, cmatrix )
  end subroutine message_cmatrix3

  subroutine message_integer ( priority_logical, debug_mode_logical, main_msg, int_variable )
    type (priority_logicals), intent(in) :: priority_logical
    type (debug_mode_logicals), intent(in) :: debug_mode_logical
    character(*), intent(in) :: main_msg
    integer(IKIND), intent(in) :: int_variable
    call message_ivector ( priority_logical, debug_mode_logical, main_msg, (/ int_variable /) )
  end subroutine message_integer

  subroutine message_integer2 ( priority_logical, main_msg, int_variable )
    type (priority_logicals), intent(in) :: priority_logical
    character(*), intent(in) :: main_msg
    integer(IKIND), intent(in) :: int_variable
    call message_ivector ( priority_logical, dbg_default, main_msg, (/ int_variable /) )  
  end subroutine message_integer2

  subroutine message_integer3 ( main_msg, int_variable )
    character(*), intent(in) :: main_msg
    integer(IKIND), intent(in) :: int_variable
    call message_ivector ( LS, dbg_default, main_msg, (/ int_variable /) )  
  end subroutine message_integer3

  subroutine message_ivector2 ( priority_logical, main_msg, ivector )    
    type (priority_logicals), intent(in) :: priority_logical
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:),intent(in) :: ivector
    call message_ivector ( priority_logical, dbg_default, main_msg, ivector )
  end subroutine message_ivector2

  subroutine message_ivector3 ( main_msg, ivector )    
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:),intent(in) :: ivector
    call message_ivector ( LS, dbg_default, main_msg, ivector )
  end subroutine message_ivector3

  subroutine message_imatrix2 ( priority_logical, main_msg, imatrix )    
    type (priority_logicals), intent(in) :: priority_logical
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:,:),intent(in) :: imatrix
    call message_imatrix ( priority_logical, dbg_default, main_msg, imatrix )
  end subroutine message_imatrix2

  subroutine message_imatrix3 ( main_msg, imatrix )    
    character(*), intent(in) :: main_msg
    integer(IKIND),dimension(:,:),intent(in) :: imatrix
    call message_imatrix ( LS, dbg_default, main_msg, imatrix )
  end subroutine message_imatrix3

  subroutine message_string2 ( priority_logical, main_msg, str_variable )
    type (priority_logicals), intent(in) :: priority_logical
    character(*), intent(in) :: main_msg, str_variable
    call message_string ( priority_logical, dbg_default, main_msg, str_variable )
  end subroutine message_string2

  subroutine message_string3 ( main_msg, str_variable )
    character(*), intent(in) :: main_msg, str_variable
    call message_string ( LS, dbg_default, main_msg, str_variable )
  end subroutine message_string3

  subroutine message_only2 ( priority_logical, main_msg )
    type (priority_logicals), intent(in) :: priority_logical
    character(*), intent(in) :: main_msg
    call message_only ( priority_logical, dbg_default, main_msg )  
  end subroutine message_only2

  subroutine message_only3 ( main_msg )
    character(*), intent(in) :: main_msg
    call message_only ( LS, dbg_default, main_msg )  
  end subroutine message_only3

end module message_mod

!CALL message(LL,dbg6,"For Pattern ",ind,", thickness ",jnd)
!CALL message(LL,dbg6,"  the FoM = ",RImageCorrelation)
