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
!RB RgVecMat is a list of g-vectors in the microscope ref frame, units of 1/A, multiplied by 2 pi
  DO ind=1,nReflections
     DO jnd=1,nReflections
        RgMatMat(ind,jnd,:)= RgVecMatT(ind,:)-RgVecMatT(jnd,:)
        RgMatMag(ind,jnd)= SQRT(DOT_PRODUCT(RgMatMat(ind,jnd,:),RgMatMat(ind,jnd,:)))
     ENDDO
  ENDDO
  !RB take the 2 pi back out of the magnitude...   
  RgMatMag = RgMatMag/TWOPI
  !RB for symmetry determination
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
           !RB Fill the symmetry relation matrix with incrementing numbers that have the sign of the imaginary part
		   WHERE (ABS(RgSumMat-RgSumMat(ind,jnd)).LE.RTolerance)
              ISymmetryRelations = Iuid*SIGN(1_IKIND,NINT(AIMAG(CUgMatNoAbs)/TINY**2))
           END WHERE
        END IF
     END DO
  END DO

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     WRITE(SPrintString,FMT='(I4,A24)') Iuid," unique structure factors"
     PRINT*,TRIM(ADJUSTL(SPrintString))
!     PRINT*,"Unique Ugs = ",Iuid
  END IF
 
  ALLOCATE(ISymmetryStrengthKey(Iuid),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SymmetryRelatedStructureFactorDetermination(", my_rank, ") error ", IErr, &
          " in ALLOCATE() ISymmetryStrengthKey"
     RETURN
  END IF
 
  ALLOCATE(CSymmetryStrengthKey(Iuid),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SymmetryRelatedStructureFactorDetermination(", my_rank, ") error ", IErr, " in ALLOCATE() CSymmetryStrengthKey"
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

  INTEGER(IKIND) :: &
       ind, jnd, knd, oddindlorentz, evenindlorentz, oddindgauss, &
       evenindgauss,imaxj, IFound, ICount, currentatom,IErr
  INTEGER(IKIND),DIMENSION(2) :: &
       IPos, ILoc
  COMPLEX(CKIND) :: &
       CVgij
  REAL(RKIND) :: &
       RMeanInnerPotentialVolts,RAtomicFormFactor, Lorentzian,Gaussian

  CALL Message("StructureFactorInitialisation",IMust,IErr)

  CUgMatNoAbs = CZERO

  DO ind=1,nReflections
     imaxj = ind
     DO jnd=1,imaxj 

        CVgij= 0.0D0

        DO iAtom=1, INAtomsUnitCell
           ICurrentAtom = IAtoms(iAtom)
           ! calculate f_e(q) as in Eq. (C.15) of Kirkland, "Advanced Computing in EM"

           SELECT CASE (IScatterFactorMethodFLAG)

           CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians

              RAtomicFormFactor = ZERO
              DO knd = 1, 3

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

           RAtomicFormFactor = RAtomicFormFactor*ROcc(iAtom)

           IF (IAnisoDebyeWallerFactorFlag.EQ.0) THEN

              IF(RDWF(iAtom).GT.10.OR.RDWF(iAtom).LT.0) THEN
                 RDWF(iAtom) = RDebyeWallerConstant
              END IF

              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-((RgMatMag(ind,jnd)/2.D0)**2)*RDWF(iAtom))

           ELSE

              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-TWOPI*DOT_PRODUCT(RgMatMat(ind,jnd,:), &
                   MATMUL( RAnisotropicDebyeWallerFactorTensor( &
                   IAnisoDWFT(iAtom),:,:), &
                   RgMatMat(ind,jnd,:))))

           END IF

           CVgij = CVgij + &
                RAtomicFormFactor * &
                EXP(-CIMAGONE* &
                DOT_PRODUCT(RgMatMat(ind,jnd,:), RrVecMat(iAtom,:)) &
                )
        ENDDO

        CUgMatNoAbs(ind,jnd)=((((TWOPI**2)* RRelativisticCorrection) / &!RB
             (PI * RVolume)) * CVgij)

     ENDDO
  ENDDO

  RMeanInnerCrystalPotential= REAL(CUgMatNoAbs(1,1))!RB

  !RB Only the lower half of the Ug matrix was calculated, this completes the upper half
  !and also doubles the values on the diagonal
  CUgMatNoAbs = CUgMatNoAbs + CONJG(TRANSPOSE(CUgMatNoAbs))!RB

  DO ind=1,nReflections!RB now halve the diagonal again
     CUgMatNoAbs(ind,ind)=CUgMatNoAbs(ind,ind)-RMeanInnerCrystalPotential!RB
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
       MessageVariable = "RBigK", &
       RVariable = RBigK)

END SUBROUTINE StructureFactorInitialisation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE StructureFactorsWithAbsorptionDetermination(IErr)         


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

   CALL Message("StructureFactorsWithAbsorptionDetermination",IMust,IErr)

  CUgMatPrime = CZERO
    
  SELECT CASE (IAbsorbFLAG)

  CASE(1)

!!$     THE PROPORTIONAL MODEL OF ABSORPTION
     
     CUgMatPrime = CUgMatNoAbs*EXP(CIMAGONE*PI/2)*(RAbsorptionPercentage/100_RKIND)!RB changed
     CUgMat =  CUgMatNoAbs+CUgMatPrime!RB
!RB   PRINT*,"CUgMat(1,1)= ",CUgMatNoAbs(1,1),CUgMatPrime(1,1)!,CUgMat(1,1)   
!RB   PRINT*,"CUgMat(2,1)= ",CUgMatNoAbs(2,1),CUgMatPrime(2,1)!,CUgMat(2,1)   
!RB   PRINT*,"CUgMat(1,2)= ",CUgMatNoAbs(1,2),CUgMatPrime(1,2)!,CUgMat(1,2)   
!RB   PRINT*,"CUgMat(3,1)= ",CUgMatNoAbs(3,1),CUgMatPrime(3,1)!,CUgMat(3,1)   
!RB   PRINT*,"CUgMat(1,3)= ",CUgMatNoAbs(1,3),CUgMatPrime(1,3)!,CUgMat(1,3)   
!RB   PRINT*,"CUgMat(3,2)= ",CUgMatNoAbs(3,2),CUgMatPrime(3,2)!,CUgMat(3,2)   
!RB   PRINT*,"CUgMat(2,3)= ",CUgMatNoAbs(2,3),CUgMatPrime(2,3)!,CUgMat(2,3)   
!RB   PRINT*,"CUgMat(4,1)= ",CUgMatNoAbs(4,1),CUgMatPrime(4,1)!,CUgMat(2,2)   
!RB   PRINT*,"CUgMat(1,4)= ",CUgMatNoAbs(1,4),CUgMatPrime(1,4)!,CUgMat(3,3) 
  
  CASE Default
 
  END SELECT
  
END SUBROUTINE StructureFactorsWithAbsorptionDetermination
  
