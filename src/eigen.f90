!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! BlochSim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all right reserved
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
!  This file is part of BlochSim.
!
!  BlochSim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  BlochSim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with BlochSim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: eigen.f90,v 1.10 2014/03/25 15:35:34 phsht Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE EigenSpectrum(isize, matU, evals, evecs, IErr)

  USE MyNumbers

  IMPLICIT NONE

  INTEGER(KIND=IKIND) isize, IErr

  COMPLEX(KIND=RKIND) matU(isize,isize), evals(isize), evecs(isize,isize)

  INTEGER(KIND=IKIND) LWORK, LRWORK, LIWORK
  COMPLEX(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: WORK
  REAL(KIND=RKIND),    DIMENSION(:), ALLOCATABLE :: RWORK
  !REAL(KIND=IKIND),    DIMENSION(:), ALLOCATABLE :: IWORK
  EXTERNAL ZGEEV

  INTEGER(IKIND) ind,jnd
  REAL(RKIND) norm

  !PRINT*,"DBG: EigenSpectrum()"

  ! ------------------------------------------------
  ! find optimum size of arrays
  ! ------------------------------------------------
  LWORK=1
  ALLOCATE(WORK(LWORK), STAT = IErr)
  ALLOCATE(RWORK(2*isize), STAT = IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ALLOCATE() for work arrays (query stage)"
     RETURN
  ENDIF

  LWORK=-1
  CALL ZGEEV('N','V', isize, matU, isize, evals, 0,1, evecs,isize, &
       WORK, LWORK, RWORK, IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ZGEEV determining work arrays"
     RETURN
  ENDIF

  LWORK = INT(WORK(1))
  !PRINT*,"DBG: LWORK=", LWORK

  ! ------------------------------------------------
  ! ALLOCATE necessary memory
  ! ------------------------------------------------
  DEALLOCATE(WORK)
  ALLOCATE(WORK(LWORK), STAT = IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ALLOCATE() for work arrays (final stage)"
     RETURN
  ENDIF

  ! ------------------------------------------------
  ! do the actual call to get the spectrum
  ! ------------------------------------------------
  CALL ZGEEV('N','V', isize, matU, isize, evals, 0,1, evecs,isize, &
       WORK, LWORK, RWORK, IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error ", IErr, " in ZGEEV"
     RETURN
  ENDIF
  DEALLOCATE(WORK,RWORK)

  !PRINT*,"DBG: evals=", evals

  RETURN

  ! this part is not needed
  ! normalize the eigenvectors
  DO ind=1,isize
     norm= SQRT(DOT_PRODUCT(evecs(ind,:),evecs(ind,:)))
     !PRINT*, "DBG: before norm:", ind,norm

     evecs(ind,:)= evecs(ind,:)/norm

     norm= SQRT(DOT_PRODUCT(evecs(ind,:),evecs(ind,:)))
     !PRINT*, "DBG: after norm:", ind,norm
  ENDDO
  
  ! orthogonality test, comment out of not needed
  DO ind=1,isize
     DO jnd=1,isize
        
        norm= SQRT(DOT_PRODUCT( &
             evecs(ind,:),evecs(jnd,:)))
        
        !PRINT*,"DBG: (ind,jnd)", ind,jnd,norm
     ENDDO
  ENDDO ! orthogonality
  
  RETURN

END SUBROUTINE EigenSpectrum
