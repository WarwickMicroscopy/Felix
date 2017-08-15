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
module alert_mod

  implicit none

  contains

  ! returns false if ierr is non-zero and also prints errors using input info 
  logical function l_alert(ierr,SCurrentProcedure,SAlertedActivity)
    use mynumbers, only : ikind
    use MyMPI, only : my_rank
    character(*),intent(in) :: SCurrentProcedure,SAlertedActivity
    integer(ikind),intent(in) :: ierr
    l_alert = .false.
    if ( ierr /= 0 ) then
      l_alert = .true.
      write(*,'(2x,i1,a,a,a,a,a)') my_rank," = rank, error in ", SCurrentProcedure,&
            "(", SAlertedActivity, ")"
    end if
  end function

end module alert_mod


!   Example output - (error handling on 4 cores from problem on felix.inp on line 27)
!
!   """"
!
!      @ -----------------------------------------------------------------
!      1 = rank, error in ReadInpFile(READ() felix.inp line number =          27)
!      1 = rank, error in felixrefine(ReadInpFile())
!      1 = rank, error in felixrefine(ABORTING)
!      @ felixrefine: 'Version: :VERSION: / :BUILD: / :AUTHOR:           '
!      @              'Date: :DATE:                                      '
!      @              '(:RLSTATUS:) multipole atom test & debug          '
!      @ -----------------------------------------------------------------
!      @ total number of MPI ranks is 004, screen messages via rank= 000
!      @ -----------------------------------------------------------------
!      3 = rank, error in ReadInpFile(READ() felix.inp line number =          27)
!      3 = rank, error in felixrefine(ReadInpFile())
!      3 = rank, error in felixrefine(ABORTING)
!      0 = rank, error in ReadInpFile(READ() felix.inp line number =          27)
!      0 = rank, error in felixrefine(ReadInpFile())
!      0 = rank, error in felixrefine(ABORTING)
!    --------------------------------------------------------------------------
!    MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD 
!    with errorcode 1.

!    NOTE: invoking MPI_ABORT causes Open MPI to kill all MPI processes.
!    You may or may not see output from other processes, depending on
!    exactly when Open MPI kills them.
!    --------------------------------------------------------------------------
!      2 = rank, error in ReadInpFile(READ() felix.inp line number =          27)
!      2 = rank, error in felixrefine(ReadInpFile())
!      2 = rank, error in felixrefine(ABORTING)
!
!   """"
!
!   Example usage:
!
!     CALL MPI_Comm_size(MPI_COMM_WORLD,p,IErr) ! get size of the current communicator
!     IF(l_alert(IErr,"felixrefine","MPI_Comm_size()")) CALL abort()
!
!     ALLOCATE(RAtomCoordinate(IMaxPossibleNAtomsUnitCell,ITHREE),STAT=IErr)
!     IF(l_alert(IErr,"felixrefine","allocate RAtomCoordinate")) CALL abort()
!
!     CALL BlochCoefficientCalculation(ind,jnd,knd,ILocalPixelCountMin,IErr)
!     IF(l_alert(IErr,"FelixFunction","BlochCoefficientCalculation")) RETURN
!
!
!
!   CALL abort() used in felixrefine program
!
!   RETURN used in any subroutine, 
!   expecting to be picked up by another IF(l_alert(IErr outside
!
!
!    
!   Used following style when printing simple variables with error messages:
!
!     IErr=1
!     WRITE(SPrintString,*) ILine
!     IF(l_alert( IErr,"ReadInpFile",&
!           "READ() felix.inp, premature end of file, line number =" // &
!           TRIM(SPrintString) )) RETURN
!
!     IF(l_alert(IErr,"ReadExperimentalImages",&
!           "OPEN() an experimental image, filename ="//TRIM(ADJUSTL(filename)))) RETURN

