!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all rights reserved
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

SUBROUTINE GMatrixInitialisation (IErr)

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$%
!!$%   Creates a matrix of every inter G vector (i.e. g1-g2) and
!!$%   their magnitudes 
!!$%
!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara

  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,jnd,IErr

  CALL Message("GMatrixInitialisation",IMust,IErr)
!Ug RgPool is a list of g-vectors in the microscope ref frame, units of 1/A, multiplied by 2 pi
  DO ind=1,nReflections
     DO jnd=1,nReflections
        RgMatMat(ind,jnd,:)= RgPoolT(ind,:)-RgPoolT(jnd,:)
        RgMatMag(ind,jnd)= SQRT(DOT_PRODUCT(RgMatMat(ind,jnd,:),RgMatMat(ind,jnd,:)))
     ENDDO
  ENDDO
  !Ug take the 2 pi back out of the magnitude...   
  RgMatMag = RgMatMag/TWOPI
  !Ug for symmetry determination
  RgSumMat = SUM(ABS(RgMatMat),3)
  
END SUBROUTINE GMatrixInitialisation

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SymmetryRelatedStructureFactorDetermination (IErr)

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$%
!!$%    Determines which structure factors are related by symmetry, by assuming 
!!$%    that two structure factors with identical absolute values are related
!!$%    (allowing for the hermiticity)
!!$%
!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara

  USE IChannels

  USE MPI
  USE MyMPI
    
  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,jnd,ierr,knd,Iuid
  CHARACTER*200 :: SPrintString

  CALL Message("SymmetryRelatedStructureFactorDetermination",IMust,IErr)

  RgSumMat = RgSumMat+ABS(REAL(CUgMatNoAbs))+ABS(AIMAG(CUgMatNoAbs))

  ISymmetryRelations = 0_IKIND 
  Iuid = 0_IKIND
  
  DO ind = 1,nReflections
     DO jnd = 1,ind
        IF(ISymmetryRelations(ind,jnd).NE.0) THEN
           CYCLE
        ELSE
           Iuid = Iuid + 1_IKIND
           !Ug Fill the symmetry relation matrix with incrementing numbers that have the sign of the imaginary part
		   WHERE (ABS(RgSumMat-RgSumMat(ind,jnd)).LE.RTolerance)
              ISymmetryRelations = Iuid*SIGN(1_IKIND,NINT(AIMAG(CUgMatNoAbs)/TINY**2))
           END WHERE
        END IF
     END DO
  END DO

!yy DO ind = 1,4!yy
!yy  WRITE(SPrintString,FMT='(I1,A1,4(1X,I2))') ind,":",ISymmetryRelations(ind,1:4)
!yy  PRINT*,TRIM(ADJUSTL(SPrintString))
!yy END DO

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     WRITE(SPrintString,FMT='(I5,A25)') Iuid," unique structure factors"
     PRINT*,TRIM(ADJUSTL(SPrintString))
!     PRINT*,"Unique Ugs = ",Iuid
  END IF
 
  ALLOCATE(IEquivalentUgKey(Iuid),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SymmetryRelatedStructureFactorDetermination(", my_rank, ") error ", IErr, &
          " in ALLOCATE() IEquivalentUgKey"
     RETURN
  END IF
 
  ALLOCATE(CUgToRefine(Iuid),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SymmetryRelatedStructureFactorDetermination(", my_rank, ") error ", IErr, " in ALLOCATE() CUgToRefine"
     RETURN
  END IF
  
END SUBROUTINE SymmetryRelatedStructureFactorDetermination

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE StructureFactorInitialisation (IErr)

  USE MyNumbers
  USE WriteToScreen

  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: ind,jnd,knd,lnd,oddindlorentz,evenindlorentz,oddindgauss, &
       evenindgauss,imaxj,IFound,ICount,currentatom,IErr
  INTEGER(IKIND),DIMENSION(2) :: IPos,ILoc
  COMPLEX(CKIND) :: CVgij
  REAL(RKIND) :: RMeanInnerPotentialVolts,RAtomicFormFactor, Lorentzian,Gaussian

  CALL Message("StructureFactorInitialisation",IMust,IErr)

  CUgMatNoAbs = CZERO

  DO ind=1,nReflections
     DO jnd=1,ind 
        CVgij= 0.0D0
        DO lnd=1,INAtomsUnitCell
           ICurrentAtom = IAtoms(lnd)!Atomic number

           SELECT CASE (IScatterFactorMethodFLAG)! calculate f_e(q) as in Eq. C.15 of Kirkland, "Advanced Computing in EM"


           CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians
              RAtomicFormFactor = ZERO
              DO knd = 1,3
                 !odd and even indicies for Lorentzian function
                 evenindlorentz = knd*2
                 oddindlorentz = knd*2 -1
                 !odd and even indicies for Gaussian function
                 evenindgauss = evenindlorentz + 6
                 oddindgauss = oddindlorentz + 6
                 !Kirkland Method uses summation of 3 Gaussians and 3 Lorentzians (summed in loop)
                 RAtomicFormFactor = RAtomicFormFactor + &
                                !3 Lorentzians
                      LORENTZIAN(RScattFactors(ICurrentAtom,oddindlorentz), RgMatMag(ind,jnd),ZERO,&
                      RScattFactors(ICurrentAtom,evenindlorentz))+ &
                                !3 Gaussians
                      GAUSSIAN(RScattFactors(ICurrentAtom,oddindgauss),RgMatMag(ind,jnd),ZERO, & 
                      1/(SQRT(2*RScattFactors(ICurrentAtom,evenindgauss))),ZERO)
              END DO

           CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
              RAtomicFormFactor = ZERO
              DO knd = 1, 4
                 !Peng Method uses summation of 4 Gaussians
                 RAtomicFormFactor = RAtomicFormFactor + &
                                !4 Gaussians
                      GAUSSIAN(RScattFactors(ICurrentAtom,knd),RgMatMag(ind,jnd),ZERO, & 
                      SQRT(2/RScattFactors(ICurrentAtom,knd+4)),ZERO)
              END DO
			  
           CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)
              RAtomicFormFactor = ZERO
              DO knd = 1, 4
                 evenindgauss = knd*2
                 oddindgauss = knd*2 -1
                 !Doyle &Turner uses summation of 4 Gaussians
                 RAtomicFormFactor = RAtomicFormFactor + &
                                !4 Gaussians
                      GAUSSIAN(RScattFactors(ICurrentAtom,oddindgauss),RgMatMag(ind,jnd),ZERO, & 
                      SQRT(2/RScattFactors(ICurrentAtom,evenindgauss)),ZERO)
              END DO

           CASE(3) ! 10 Parameter method with Scattering Parameters from Lobato et al. 2014
              RAtomicFormFactor = ZERO
              DO knd = 1,5
                 evenindlorentz=knd+5
                 RAtomicFormFactor = RAtomicFormFactor + &
                      LORENTZIAN(RScattFactors(ICurrentAtom,knd)* &
                      (TWO+RScattFactors(ICurrentAtom,evenindlorentz)*(RgMatMag(ind,jnd)**TWO)), &
                      ONE, &
                      RScattFactors(ICurrentAtom,evenindlorentz)*(RgMatMag(ind,jnd)**TWO),ZERO)
              END DO

           END SELECT

           ! initialize potential as in Eq. (6.10) of Kirkland

           RAtomicFormFactor = RAtomicFormFactor*ROcc(lnd)
           IF (IAnisoDebyeWallerFactorFlag.EQ.0) THEN
              IF(RDWF(lnd).GT.10.OR.RDWF(lnd).LT.0) THEN
                 RDWF(lnd) = RDebyeWallerConstant
              END IF
              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-((RgMatMag(ind,jnd)/2.D0)**2)*RDWF(lnd))
           ELSE
              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-TWOPI*DOT_PRODUCT(RgMatMat(ind,jnd,:), &
                   MATMUL( RAnisotropicDebyeWallerFactorTensor( &
                   IAnisoDWFT(lnd),:,:), &
                   RgMatMat(ind,jnd,:))))
           END IF
           CVgij = CVgij + RAtomicFormFactor * &
                EXP(-CIMAGONE* &
                DOT_PRODUCT(RgMatMat(ind,jnd,:), RrVecMat(lnd,:)) &
                )
        ENDDO
        CUgMatNoAbs(ind,jnd)=((((TWOPI**2)* RRelativisticCorrection) / &!Ug
             (PI * RVolume)) * CVgij)
     ENDDO
  ENDDO

  RMeanInnerCrystalPotential= REAL(CUgMatNoAbs(1,1))!Ug

  !NB Only the lower half of the Ug matrix was calculated, this completes the upper half
  !and also doubles the values on the diagonal
  CUgMatNoAbs = CUgMatNoAbs + CONJG(TRANSPOSE(CUgMatNoAbs))!Ug

  DO ind=1,nReflections!Ug now halve the diagonal again
     CUgMatNoAbs(ind,ind)=CUgMatNoAbs(ind,ind)-RMeanInnerCrystalPotential!Ug
  ENDDO
  
  RMeanInnerPotentialVolts = RMeanInnerCrystalPotential*(((RPlanckConstant**2)/ &
       (TWO*RElectronMass*RElectronCharge*TWOPI**2))*&
       RAngstromConversion*RAngstromConversion)

  CALL Message("StructureFactorInitialisation",IMoreInfo,IErr, &
       MessageVariable = "RMeanInnerCrystalPotential", &
       RVariable = RMeanInnerCrystalPotential)

  CALL Message("StructureFactorInitialisation",IMoreInfo,IErr, &
       MessageVariable = "RMeanInnerPotentialVolts", &
       RVariable = RMeanInnerPotentialVolts)

  !Now initialisation calls the Ug calculation subroutines
  !Used to be in Setup

  !--------------------------------------------------------------------
  ! high-energy approximation (not HOLZ compatible)
  !--------------------------------------------------------------------

  RBigK= SQRT(RElectronWaveVectorMagnitude**2 + RMeanInnerCrystalPotential)

  CALL Message("StructureFactorInitialisation",IInfo,IErr, &
       MessageVariable = "RBigK", RVariable = RBigK)

END SUBROUTINE StructureFactorInitialisation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE StructureFactorsWithAbsorption(IErr)         
!RB this is a lot of bumf for 3 lines of code

  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 
  
  INTEGER(IKIND) :: IErr,ind
  CHARACTER*200 :: SPrintString

   CALL Message("StructureFactorsWithAbsorption",IMust,IErr)

  CUgMatPrime = CZERO
    
  SELECT CASE (IAbsorbFLAG)

  CASE(1)

!!$     THE PROPORTIONAL MODEL OF ABSORPTION
     
     CUgMatPrime = CUgMatNoAbs*EXP(CIMAGONE*PI/2)*(RAbsorptionPercentage/100_RKIND)!Ug
     CUgMat =  CUgMatNoAbs+CUgMatPrime!Ug
!  WRITE(SPrintString,FMT='(A2,8(1X,F5.2))') "1:",CUgMatNoAbs(1,1),CUgMatNoAbs(1,2),CUgMatNoAbs(1,3),CUgMatNoAbs(1,4)
!  PRINT*,TRIM(ADJUSTL(SPrintString))
!  WRITE(SPrintString,FMT='(A2,8(1X,F5.2))') "2:",CUgMatNoAbs(2,1),CUgMatNoAbs(2,2),CUgMatNoAbs(2,3),CUgMatNoAbs(2,4)
!  PRINT*,TRIM(ADJUSTL(SPrintString))
!  WRITE(SPrintString,FMT='(A2,8(1X,F5.2))') "3:",CUgMatNoAbs(3,1),CUgMatNoAbs(3,2),CUgMatNoAbs(3,3),CUgMatNoAbs(3,4)
!  PRINT*,TRIM(ADJUSTL(SPrintString))
!  WRITE(SPrintString,FMT='(A2,8(1X,F5.2))') "4:",CUgMatNoAbs(4,1),CUgMatNoAbs(4,2),CUgMatNoAbs(4,3),CUgMatNoAbs(4,4)
!  PRINT*,TRIM(ADJUSTL(SPrintString))
  
  CASE Default
 
  END SELECT
  
END SUBROUTINE StructureFactorsWithAbsorption
  
