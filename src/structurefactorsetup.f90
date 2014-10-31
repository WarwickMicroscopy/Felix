!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
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
!  This file is part of felixsim.
!
!  felixsim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  felixsim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with felixsim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-------------------------------------------------------------------------------
!setup structure factors for each material - to be used in Bloch Calculation
!-------------------------------------------------------------------------------

SUBROUTINE StructureFactorSetup(IErr)

  USE MyNumbers
  
  USE IPara; USE RPara ; USE CPara
  USE BlochPara

  USE MPI
  USE MyMPI


  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,jnd

  COMPLEX(CKIND),DIMENSION(:,:), ALLOCATABLE :: &
       CZeroMat

  !--------------------------------------------------------------------
  ! Calculate Reflection Matrix
  !--------------------------------------------------------------------
  ALLOCATE( &  
       RgMatMat(nReflections,nReflections,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     !call error function
     RETURN
  ENDIF
       
  ALLOCATE( &  
       RgMatMag(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     !call error function
     RETURN
  ENDIF

  CALL GMatrixInitialisation (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
          " in GMatrixInitialisation"
     !error function call
     RETURN
  ENDIF

  !--------------------------------------------------------------------
  ! calculating Ug matrix
  !--------------------------------------------------------------------

  !Allocate memory for Ug Matrix

  ALLOCATE( & 
       CUgMat(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     !call error function
     RETURN
  ENDIF  

  ALLOCATE( & 
       CUgMatPrime(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     !call error function
     RETURN
  ENDIF       
 
  ALLOCATE( & 
       CZeroMat(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     !refinemain was here
     PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CZeroMat"
     !call error function
     RETURN
  ENDIF

  CALL StructureFactorInitialisation (IErr, CZeroMat)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
          " in StructureFactorInitialisation"
     !call error function
     RETURN
  ENDIF


END SUBROUTINE StructureFactorSetup
