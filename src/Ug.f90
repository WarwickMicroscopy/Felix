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
 ! This is now redundant, moved to StructureFactorSetup
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
  !Ug RgPool is a list of g-vectors in the microscope ref frame, units of 1/A
  ! Note that reciprocal lattice vectors dot not have two pi included, we are using the optical convention exp(2*pi*i*g.r)
  DO ind=1,nReflections
     DO jnd=1,nReflections
        RgMatrix(ind,jnd,:)= RgPool(ind,:)-RgPool(jnd,:)
        RgMatrixMagnitude(ind,jnd)= SQRT(DOT_PRODUCT(RgMatrix(ind,jnd,:),RgMatrix(ind,jnd,:)))
     ENDDO
  ENDDO
  !For symmetry determination, only in Ug refinement
  IF (IRefineMode(1).EQ.1 .OR. IRefineMode(12).EQ.1) THEN
    RgSumMat = SUM(ABS(RgMatrix),3)
  END IF
  
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
  !!This is now redundant, moved to SetupUgsToRefine
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

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     WRITE(SPrintString,FMT='(I5,A25)') Iuid," unique structure factors"
     PRINT*,TRIM(ADJUSTL(SPrintString))
!     PRINT*,"Unique Ugs = ",Iuid
  END IF
 
  ALLOCATE(IEquivalentUgKey(Iuid),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SymmetryRelatedStructureFactorDetermination(",my_rank,")error allocating IEquivalentUgKey"
     RETURN
  END IF
  ALLOCATE(CUgToRefine(Iuid),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SymmetryRelatedStructureFactorDetermination(",my_rank,")error allocating CUgToRefine"
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
       evenindgauss,currentatom,IErr
  INTEGER(IKIND),DIMENSION(2) :: IPos,ILoc
  COMPLEX(CKIND) :: CVgij
  REAL(RKIND) :: RMeanInnerPotentialVolts,RScatteringFactor,Lorentzian,Gaussian,Kirkland,RScattFacToVolts
  CHARACTER*200 :: SPrintString
  
  CALL Message("StructureFactorInitialisation",IMust,IErr)

  !Conversion factor from scattering factors to volts. h^2/(2pi*m0*e*CellVolume), see e.g. Kirkland eqn. C.5
  RScattFacToVolts=(RPlanckConstant**2)*(RAngstromConversion**3)/(TWOPI*RElectronMass*RElectronCharge*RVolume)
  !Calculate Ug matrix
  CUgMatNoAbs = CZERO
  DO ind=1,nReflections
    DO jnd=1,ind
      !The Fourier component of the potential Vg goes in location (i,j)
      CVgij= 0.0D0!this is in Volts
      DO lnd=1,INAtomsUnitCell
        ICurrentAtom = IAtomicNumber(lnd)!Atomic number
        SELECT CASE (IScatterFactorMethodFLAG)! calculate f_e(q)

        CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians, NB Kirkland scattering factor is in Angstrom units
          RScatteringFactor = Kirkland(IAtomicNumber(lnd),RgMatrixMagnitude(ind,jnd))
  
        CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
          RScatteringFactor = ZERO
          DO knd = 1, 4
            !Peng Method uses summation of 4 Gaussians
            RScatteringFactor = RScatteringFactor + &
              GAUSSIAN(RScattFactors(ICurrentAtom,knd),RgMatrixMagnitude(ind,jnd),ZERO, & 
              SQRT(2/RScattFactors(ICurrentAtom,knd+4)),ZERO)
          END DO
			  
        CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)
          RScatteringFactor = ZERO
          DO knd = 1, 4
            evenindgauss = knd*2
            oddindgauss = knd*2 -1
            !Doyle &Turner uses summation of 4 Gaussians
            RScatteringFactor = RScatteringFactor + &
              GAUSSIAN(RScattFactors(ICurrentAtom,oddindgauss),RgMatrixMagnitude(ind,jnd),ZERO, & 
              SQRT(2/RScattFactors(ICurrentAtom,evenindgauss)),ZERO)
          END DO

        CASE(3) ! 10 Parameter method with Scattering Parameters from Lobato et al. 2014
          RScatteringFactor = ZERO
          DO knd = 1,5
            evenindlorentz=knd+5
            RScatteringFactor = RScatteringFactor + &
              LORENTZIAN(RScattFactors(ICurrentAtom,knd)* &
             (TWO+RScattFactors(ICurrentAtom,evenindlorentz)*(RgMatrixMagnitude(ind,jnd)**TWO)),ONE, &
              RScattFactors(ICurrentAtom,evenindlorentz)*(RgMatrixMagnitude(ind,jnd)**TWO),ZERO)
          END DO

        END SELECT
        ! Occupancy
        RScatteringFactor = RScatteringFactor*ROccupancy(lnd)
        IF (IAnisoDebyeWallerFactorFlag.EQ.0) THEN
          IF(RIsoDW(lnd).GT.10.OR.RIsoDW(lnd).LT.0) THEN
            RIsoDW(lnd) = RDebyeWallerConstant
          END IF
          !Isotropic D-W factor exp(-B sin(theta)^2/lamda^2) = exp(-Bs^2)=exp(-Bg^2/16pi^2), see e.g. Bird&King
          RScatteringFactor = RScatteringFactor*EXP(-RIsoDW(lnd)*(RgMatrixMagnitude(ind,jnd)**2)/(4*TWOPI**2) )
        ELSE!this will need sorting out, not sure if it works
          RScatteringFactor = RScatteringFactor * &
            EXP(-DOT_PRODUCT(RgMatrix(ind,jnd,:), &
            MATMUL( RAnisotropicDebyeWallerFactorTensor( &
            RAnisoDW(lnd),:,:),RgMatrix(ind,jnd,:))))
        END IF
		!The structure factor equation, CVgij in Volts
        CVgij = CVgij + RScattFacToVolts*RScatteringFactor * EXP(-CIMAGONE* &
        DOT_PRODUCT(RgMatrix(ind,jnd,:), RAtomCoordinate(lnd,:)) )/RAngstromConversion
      ENDDO
	  !This is actually still the Vg matrix, converted at the end to Ug
      CUgMatNoAbs(ind,jnd)=CVgij
      !CUgMatNoAbs(ind,jnd)=((((TWOPI**2)*RRelativisticCorrection)/(PI*RVolume))*CVgij)
    ENDDO
  ENDDO
  RMeanInnerPotential= REAL(CUgMatNoAbs(1,1))
  IF(my_rank.EQ.0) THEN
    WRITE(SPrintString,FMT='(A20,F5.2,1X,A6)') "MeanInnerPotential= ",RMeanInnerPotential," Volts"
	PRINT*,TRIM(ADJUSTL(SPrintString))
  END IF
  DO ind=1,nReflections!Take the Mean Inner Potential off the diagonal 
     CUgMatNoAbs(ind,ind)=CUgMatNoAbs(ind,ind)-RMeanInnerPotential
  ENDDO
  !NB Only the lower half of the Vg matrix was calculated, this completes the upper half
  CUgMatNoAbs = CUgMatNoAbs + CONJG(TRANSPOSE(CUgMatNoAbs))
  !Now convert to Ug=Vg*(2*m*e/h^2)
  CUgMatNoAbs=CUgMatNoAbs*TWO*RElectronMass*RRelativisticCorrection*RElectronCharge/(RPlanckConstant**2)
  !Divide U0 by 10^20 to convert Planck constant to A 
  CUgMatNoAbs=CUgMatNoAbs/(RAngstromConversion**2)
  
  !Alternative way of calculating the mean inner potential as the sum of scattering factors at g=0 multiplied by h^2/(2pi*m0*e*CellVolume)
  !RMeanInnerPotential=ZERO
  !DO ind=1,INAtomsUnitCell
  !  RMeanInnerPotential = RMeanInnerPotential+Kirkland(IAtomicNumber(ind),ZERO)/RAngstromConversion
  !END DO
  !RMeanInnerPotential = RMeanInnerPotential*RScattFacToVolts
  !PRINT*,"MeanInnerPotential2=",RMeanInnerPotential
  CALL Message("StructureFactorInitialisation",IMoreInfo,IErr, &
       MessageVariable = "RMeanInnerPotential", &
       RVariable = RMeanInnerPotential)

  !--------------------------------------------------------------------
  ! high-energy approximation (not HOLZ compatible)
  !Wave vector in crystal
  !K^2=k^2+U0
  RBigK= SQRT(RElectronWaveVectorMagnitude**2 + REAL(CUgMatNoAbs(1,1)))
  CALL Message("StructureFactorInitialisation",IInfo,IErr, &
       MessageVariable = "RBigK", RVariable = RBigK)
  IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
    PRINT*,"RBigK=",RBigK
  END IF
  
  !Absorption
  CUgMatPrime = CZERO
    
  SELECT CASE (IAbsorbFLAG)

  CASE(1)
!!$ Proportional
    CUgMatPrime = CUgMatNoAbs*EXP(CIMAGONE*PI/2)*(RAbsorptionPercentage/100_RKIND)!Ug
    CUgMat =  CUgMatNoAbs+CUgMatPrime!Ug

  CASE(2)
  !!$ Bird & King

  
  CASE Default
 
  END SELECT
  
  IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
   PRINT*,"Ug matrix, including absorption (nm^-2)"
	DO ind =1,20
     WRITE(SPrintString,FMT='(3(1X,I3),A1,8(1X,F6.2,F6.2))') NINT(Rhkl(ind,:)),":",100*CUgMat(ind,1:8)
     PRINT*,TRIM(SPrintString)
    END DO
  END IF	   
	   
END SUBROUTINE StructureFactorInitialisation