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

SUBROUTINE StructureFactorSetup(IErr)

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$%
!!$%     Calculate g-vector matrix (All inter g vectors) and from them
!!$%     the Structure factors which will be available
!!$%
!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  USE WriteToScreen
  USE MyNumbers
  USE IConst

  USE IPara; USE RPara ; USE CPara
  USE BlochPara

  USE MPI
  USE MyMPI


  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr

  CALL Message("StructureFactorSetup",IMust,IErr)

  !--------------------------------------------------------------------
  ! Calculate Reflection Matrix
  !--------------------------------------------------------------------

  ALLOCATE(RgMatMat(nReflections,nReflections,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RgMatMat"
     !call error function
     RETURN
  END IF
       
  ALLOCATE(RgMatMag(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RgMatMag"
     !call error function
     RETURN
  END IF
!RB  NB Also deallocated in felixfunction!!!
  !RB Matrix for sum of indices - for symmetry equivalence  
  ALLOCATE(RgSumMat(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(",my_rank,") error ",IErr, &
          "in ALLOCATE() of DYNAMIC variables RgSumMat"
     !call error function
     RETURN
  END IF  

  CALL GMatrixInitialisation (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
          " in GMatrixInitialisation"
     !error function call
     RETURN
  END IF

  !--------------------------------------------------------------------
  ! calculating Ug matrix
  !--------------------------------------------------------------------
!RB  NB Also deallocated in felixfunction!!!
  !Allocate memory for Ug Matrix
  !RB Matrix that is sum of real+abs
  !PRINT*,"Allocating CUgMat,CUgMatNoAbs,CUgMatPrime in structurefactorsetup"
  ALLOCATE(CUgMat(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(",my_rank,") error", IErr, &
          "in ALLOCATE() of DYNAMIC variables CUgMat"
     !call error function
     RETURN
  END IF
  
  !RB Matrix without absorption
  ALLOCATE(CUgMatNoAbs(nReflections,nReflections), &!RB
       STAT=IErr)!RB
  IF( IErr.NE.0 ) THEN!RB
     PRINT*,"StructureFactorSetup(",my_rank,") error",IErr, &!RB
          "in ALLOCATE() of DYNAMIC variables CUgMatNoAbs"!RB
     !call error function
     RETURN!RB
  END IF  !RB

  !RB Matrix for absorption  
  ALLOCATE(CUgMatPrime(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(",my_rank,") error ",IErr, &
          "in ALLOCATE() of DYNAMIC variables CUgMatPrime"
     !call error function
     RETURN
  END IF  

  CALL StructureFactorInitialisation (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(",my_rank,") error", IErr, &
          "in StructureFactorInitialisation"
     !call error function
     RETURN
  END IF

  DEALLOCATE(RgMatMat,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation of RgMatMat"
     RETURN
  END IF
  
  DEALLOCATE(RgMatMag,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation of RgMatMag"
     RETURN
  END IF

  DEALLOCATE(RrVecMat,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in DEALLOCATE() "
     RETURN
  ENDIF

END SUBROUTINE StructureFactorSetup
