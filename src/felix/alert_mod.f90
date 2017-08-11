!>
!! Module-description: 
!!
MODULE alert_mod

  IMPLICIT NONE

  CONTAINS

  ! returns false if ierr is non-zero and also prints errors using input info 
  logical function l_alert(ierr,SCurrentProcedure,SAlertedActivity)
    use mynumbers, only : ikind
    character(*),intent(in) :: SCurrentProcedure,SAlertedActivity
    integer(ikind),intent(in) :: ierr
    l_alert = .false.
    if ( ierr /= 0 ) then
      l_alert = .true.
      call alert_message(SCurrentProcedure,SAlertedActivity)
    end if
  end function

  ! error message - to be used instead of print & my_rank in code
  subroutine alert_message(SCurrentProcedure,SAlertedActivity)
    use mympi, only : my_rank
    character(*),intent(in) :: SCurrentProcedure,SAlertedActivity
    write(*,'(2x,i1,a,a,a,a)') my_rank," = rank, error in ",SCurrentProcedure,&
          " after attempting to ",SAlertedActivity
  end subroutine

  ! error message - to be used instead of print & my_rank in code
  subroutine error_message(SCurrentProcedure,error_msg)
    use mympi, only : my_rank
    character(*),intent(in) :: SCurrentProcedure,error_msg

    write(*,'(2x,i1,a,a,a,a)') my_rank," = rank, error in ", SCurrentProcedure, &
          ": ",error_msg
  end subroutine

end module alert_mod
