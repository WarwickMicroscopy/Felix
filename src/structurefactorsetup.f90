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

SUBROUTINE StructureFactorSetup

USE MyNumbers

USE IPara; USE RPara ; USE CPara;
USE BlochPara

IMPLICIT NONE

INTEGER(IKIND) :: IErr,ind,jnd

COMPLEX(CKIND),DIMENSION(:,:), ALLOCATABLE :: &
       CZeroMat
INTEGER(IKIND),DIMENSION(2) :: ILoc

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

  CALL StructureFactorCalculation (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
          " in StructureFactorCalculation"
     !call error function
     RETURN
  ENDIF
  
  !--------------------------------------------------------------------
  ! high-energy approximation (not HOLZ compatible)
  !--------------------------------------------------------------------
  
  RBigK= SQRT(RElectronWaveVectorMagnitude**2 + RMeanInnerCrystalPotential)

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"StructureFactorSetup(", my_rank, ") BigK=", RBigK
  END IF
    
  IF(IAbsorbFLAG.GT.1) THEN
     CZeroMAT = CZERO
     
     DO ind = 1,nReflections
        DO jnd = 1,ind
           CZeroMAT(ind,jnd) = CONE
        END DO
     END DO
     
     ALLOCATE( &  
          ISymmetryRelations(nReflections,nReflections), &
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        !refinemain was here
        PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
        !call error function
        RETURN
     ENDIF
     
     ISymmetryRelations = ISymmetryRelations*CZeroMat  
     
     CALL DetermineSymmetryRelatedUgs (IErr)
     IF( IErr.NE.0 ) THEN
        !refinemain was here
        PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
             " in DetermineSymmetryRelatedUgs"
        !call error function
        RETURN
     ENDIF
     
     ALLOCATE( &  
          RUniqueUgPrimeValues((SIZE(ISymmetryStrengthKey,DIM=1))), &
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        !refinemain was here
        PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
        !call error function
        RETURN
     ENDIF
     
     DO ind = 1,(SIZE(ISymmetryStrengthKey,DIM=1))
        ILoc = MINLOC(ABS(ISymmetryRelations-ind))
        ISymmetryStrengthKey(ind,:) = ILoc
     END DO
  END IF
  
  IF(IAbsorbFLAG.GE.1) THEN
     
     CALL UgAddAbsorption(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
             " in UgAddAbsorption"
        !call error function
        RETURN
     ENDIF
     IF(IAbsorbFLAG.GE.2) THEN
        DO ind = 2,(SIZE(ISymmetryStrengthKey,DIM=1))
           WHERE (ISymmetryRelations.EQ.ind)
              CUgMatPrime = RUniqueUgPrimeValues(ind)*CIMAGONE
           END WHERE
        END DO
        DO ind = 1,nReflections
           CUgMatPrime(ind,ind) = RUniqueUgPrimeValues(1)*CIMAGONE
        END DO
     END IF
     CUgMat =  CUgMat+CUgMatPrime
  
  ENDIF       

END SUBROUTINE StructureFactorSetup
