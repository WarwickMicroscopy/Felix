!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! BlochSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
!
! (C) 2013/14, all right reserved
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Log: eigen.f90,v $
! Revision 1.10  2014/03/25 15:35:34  phsht
! included consistent start of each source file and GPL statements
!
! Revision 1.9  2013/11/27 12:30:21  phsht
! more consistency in MPI output routines and their error handling
!
! Revision 1.8  2013/10/02 20:43:27  phsht
! replaced old input names with new names;
! structure of input file still needs changing
!
! Revision 1.7  2013/10/02 08:26:08  phsht
! removed two DBG lines
!
! Revision 1.6  2013/10/01 20:38:01  phsht
! trying to write as a single loop over all pixels
!
! Revision 1.5  2013/09/23 16:52:29  phsht
! OPEN, write and CLOSE of data files now implemented
!
! Revision 1.4  2013/09/19 11:19:33  phsht
! eigenvector calculation now working
!
! Revision 1.3  2013/09/18 16:08:42  phsht
! work with Keiths to include changes done with Richard yesterday and
! also to continue until intensities; reached these, but don't quite agree
! yet. Nevertheless, up to the diagonalization, things do work! And eigenvalues
! are also correct, just eigenvectors not yet.
!
! Revision 1.2  2013/09/17 17:00:21  phsht
! diagonalization by LAPACK routines now included
!
! Revision 1.1  2013/09/11 14:26:48  phsht
! interface to an eigensolver routine
!
! Revision 1.1  2011/05/06 08:13:09  phsht
! 1st installement
!
! Revision 1.3  2010/10/26 14:28:07  phrkaj
! Replaced !! with !, fixed the OpenOutputRGamma of negative energy in inout.f90
!
! Revision 1.2  2010/10/26 09:43:39  phrkaj
! Deleted naive debugging statements, got rid of ILevelflag and IConvflag, deleted old logs
!
! Revision 1.1.1.1  2010/10/22 12:23:38  phsht
! ParaTMM
!
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
