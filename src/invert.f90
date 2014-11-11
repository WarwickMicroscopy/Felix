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
! $Id: invert.f90,v 1.3 2014/03/25 15:35:34 phsht Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE INVERT(MatrixSize,Matrix,InvertedMatrix,IErr)  
!
!   Invert an M*M Complex Matrix
!   Matrix: the Matrix (Destroyed)
!   InvertedMatrix: the Inverse
  USE MyNUMBERS
  
  
  IMPLICIT NONE
  
  INTEGER :: MatrixSize, LWORK, INFO, I, IErr
  COMPLEX(KIND=CKIND), DIMENSION(1:MatrixSize,1:MatrixSize) :: Matrix
  COMPLEX(KIND=CKIND), DIMENSION(1:MatrixSize,1:MatrixSize) :: InvertedMatrix
  
  INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
  COMPLEX(KIND=CKIND), DIMENSION(:), ALLOCATABLE :: WORK
  
  !PRINT*,"DBG: Invert()"
!!$
!!$  B = CZERO 
!!$  DO I=1,M
!!$     B(I,I) = CONE
!!$  END DO
  !INFO=0
  
  ALLOCATE(IPIV(MatrixSize),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Invert(): ERR in ALLOCATE(IPIV(MatrixSize)) statement, MatrixSize=", MatrixSize
     RETURN
  ENDIF
  
  CALL ZGETRF(MatrixSize,MatrixSize,Matrix,MatrixSize,IPIV,IErr)
  LWORK = MatrixSize*MatrixSize
  !LWORK = 0
  IF ( IErr.NE.0 ) THEN
     PRINT *,'Invert() : Datatype Error: IFAIL=',INFO
     RETURN
  END IF
  ALLOCATE(WORK(LWORK),STAT=IErr)   
  IF( IErr.NE.0 ) THEN
     PRINT*,"Invert(): ERR in ALLOCATE(WORK(LWORK)) statement, LWORK=", LWORK
     RETURN
  ENDIF
  
  CALL ZGETRI(MatrixSize,Matrix,MatrixSize,IPIV,WORK,LWORK,IErr)
  IF ( IErr.NE.0 ) THEN
     PRINT *,'Inversion Error: IFAIL=',INFO
     RETURN
  END IF
  DEALLOCATE(IPIV,WORK,STAT=IErr)
  IF ( IErr.NE.0 ) THEN
     PRINT *,'Invert : Deallocation Error',INFO
     RETURN
  END IF
  !DEALLOCATE(IPIV,WORK)
  InvertedMatrix = Matrix  
  RETURN
END SUBROUTINE INVERT


