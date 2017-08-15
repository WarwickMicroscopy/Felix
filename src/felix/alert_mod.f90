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
    write(*,'(2x,i1,a,a,a,a,a)') my_rank," = rank, error in ", SCurrentProcedure,&
          "(", SAlertedActivity, ")"
  end subroutine

  ! error message - to be used instead of print & my_rank in code
  subroutine error_message(SCurrentProcedure,error_msg)
    use mympi, only : my_rank
    character(*),intent(in) :: SCurrentProcedure,error_msg

    write(*,'(2x,i1,a,a,a,a)') my_rank," = rank, error in ", SCurrentProcedure, &
          ": ",error_msg
  end subroutine

end module alert_mod
