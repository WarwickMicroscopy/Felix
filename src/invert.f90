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
! $Id: invert.f90,v 1.3 2014/03/25 15:35:34 phsht Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Log: invert.f90,v $
! Revision 1.3  2014/03/25 15:35:34  phsht
! included consistent start of each source file and GPL statements
!
! Revision 1.2  2013/09/18 16:08:42  phsht
! work with Keiths to include changes done with Richard yesterday and
! also to continue until intensities; reached these, but don't quite agree
! yet. Nevertheless, up to the diagonalization, things do work! And eigenvalues
! are also correct, just eigenvectors not yet.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE INVERT(M,A,B)  
!
!   Invert an M*M Complex Matrix
!   A: the Matrix (Destroyed)
!   B: the Inverse
  USE MyNUMBERS
  IMPLICIT NONE
  
  INTEGER :: M, LWORK, INFO, I, IErr
  COMPLEX(KIND=RKIND), DIMENSION(1:M,1:M) :: A
  COMPLEX(KIND=RKIND), DIMENSION(1:M,1:M) :: B
  
  INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
  COMPLEX(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: WORK
  
  PRINT*,"DBG: Invert()"

  B = CZERO 
  DO I=1,M
     B(I,I) = CONE
  END DO
  INFO=0
  
  ALLOCATE(IPIV(M),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Invert(): ERR in ALLOCATE(IPIV(M)) statement, M=", M
     STOP
  ENDIF
  
  CALL ZGETRF(M,M,A,M,IPIV,INFO)
  LWORK = M*M
  
  ALLOCATE(WORK(LWORK),STAT=IErr)   
  IF( IErr.NE.0 ) THEN
     PRINT*,"Invert(): ERR in ALLOCATE(WORK(LWORK)) statement, LWORK=", LWORK
     STOP
  ENDIF
  
  CALL ZGETRI(M,A,M,IPIV,WORK,LWORK,INFO)
  IF ( INFO.NE.0 ) THEN
     PRINT *,'Inversion Error: IFAIL=',INFO
     STOP
  END IF
  DEALLOCATE(IPIV,WORK)
  B = A  
  RETURN
END SUBROUTINE INVERT


