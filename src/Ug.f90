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
  ALLOCATE(CUniqueUg(Iuid),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SymmetryRelatedStructureFactorDetermination(",my_rank,")error allocating CUniqueUg"
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

  INTEGER(IKIND) :: ind,jnd,knd,lnd,mnd,oddindlorentz,evenindlorentz,oddindgauss, &
       evenindgauss,currentatom,IErr,Iuid
  INTEGER(IKIND),DIMENSION(2) :: IPos,ILoc
  COMPLEX(CKIND) :: CVgij
  REAL(RKIND) :: RMeanInnerPotentialVolts,RScatteringFactor,Lorentzian,Gaussian,Kirkland,&
        RScattFacToVolts
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
        ICurrentZ = IAtomicNumber(lnd)!Atomic number
        RCurrentG = RgMatrixMagnitude(ind,jnd)!g-vector magnitude, global variable
        SELECT CASE (IScatterFactorMethodFLAG)! calculate f_e(q)

        CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians, NB Kirkland scattering factor is in Angstrom units
		  !NB atomic number and g-vector passed as global variables
          RScatteringFactor = Kirkland(RgMatrixMagnitude(ind,jnd))
  
        CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
          RScatteringFactor = ZERO
          DO knd = 1, 4
            !Peng Method uses summation of 4 Gaussians
            RScatteringFactor = RScatteringFactor + &
              GAUSSIAN(RScattFactors(ICurrentZ,knd),RgMatrixMagnitude(ind,jnd),ZERO, & 
              SQRT(2/RScattFactors(ICurrentZ,knd+4)),ZERO)
          END DO
			  
        CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)
          RScatteringFactor = ZERO
          DO knd = 1, 4
            evenindgauss = knd*2
            oddindgauss = knd*2 -1
            !Doyle &Turner uses summation of 4 Gaussians
            RScatteringFactor = RScatteringFactor + &
              GAUSSIAN(RScattFactors(ICurrentZ,oddindgauss),RgMatrixMagnitude(ind,jnd),ZERO, & 
              SQRT(2/RScattFactors(ICurrentZ,evenindgauss)),ZERO)
          END DO

        CASE(3) ! 10 Parameter method with Scattering Parameters from Lobato et al. 2014
          RScatteringFactor = ZERO
          DO knd = 1,5
            evenindlorentz=knd+5
            RScatteringFactor = RScatteringFactor + &
              LORENTZIAN(RScattFactors(ICurrentZ,knd)* &
             (TWO+RScattFactors(ICurrentZ,evenindlorentz)*(RgMatrixMagnitude(ind,jnd)**TWO)),ONE, &
              RScattFactors(ICurrentZ,evenindlorentz)*(RgMatrixMagnitude(ind,jnd)**TWO),ZERO)
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

  ! high-energy approximation (not HOLZ compatible)
  !Wave vector in crystal
  !K^2=k^2+U0
  RBigK= SQRT(RElectronWaveVectorMagnitude**2 + REAL(CUgMatNoAbs(1,1)))
  CALL Message("StructureFactorInitialisation",IInfo,IErr, &
       MessageVariable = "RBigK", RVariable = RBigK)
  IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
    WRITE(SPrintString,FMT='(A4,F5.1,A10)') "K = ",RBigK," Angstroms"
	PRINT*,TRIM(ADJUSTL(SPrintString))
  END IF

  !--------------------------------------------------------------------
  !Count equivalent Ugs
  !Equivalent Ug's are identified by the sum of their abs(indices)plus the sum of abs(Ug)'s with no absorption
  RgSumMat = SUM(ABS(RgMatrix),3)+RgMatrixMagnitude+ABS(REAL(CUgMatNoAbs))+ABS(AIMAG(CUgMatNoAbs))
  ISymmetryRelations = 0_IKIND 
  Iuid = 0_IKIND 
  DO ind = 1,nReflections
    DO jnd = 1,ind
      IF(ISymmetryRelations(ind,jnd).NE.0) THEN
        CYCLE
      ELSE
        Iuid = Iuid + 1_IKIND
        !Ug Fill the symmetry relation matrix with incrementing numbers that have the sign of the imaginary part
	    WHERE (ABS(RgSumMat-ABS(RgSumMat(ind,jnd))).LE.RTolerance)
          ISymmetryRelations = Iuid*SIGN(1_IKIND,NINT(AIMAG(CUgMatNoAbs)/TINY**2))
        END WHERE
      END IF
    END DO
  END DO

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
    WRITE(SPrintString,FMT='(I5,A25)') Iuid," unique structure factors"
    PRINT*,TRIM(ADJUSTL(SPrintString))
  END IF
  IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
	PRINT*,"hkl: symmetry matrix"
    DO ind =1,20
      WRITE(SPrintString,FMT='(3(1X,I3),A1,12(2X,I3))') NINT(Rhkl(ind,:)),":",ISymmetryRelations(ind,1:12)
      PRINT*,TRIM(SPrintString)
    END DO
  END IF

!Link each key with its Ug, from 1 to the number of unique Ug's Iuid
  ALLOCATE(IEquivalentUgKey(Iuid),STAT=IErr)
  ALLOCATE(CUniqueUg(Iuid),STAT=IErr)
  IF( IErr.NE.0 ) THEN
    PRINT*,"SetupUgsToRefine(",my_rank,")error allocating IEquivalentUgKey or CUniqueUg"
    RETURN
  END IF
  DO ind = 1,Iuid
    ILoc = MINLOC(ABS(ISymmetryRelations-ind))
    IEquivalentUgKey(ind) = ind
    CUniqueUg(ind) = CUgMatNoAbs(ILoc(1),ILoc(2))
  END DO
  
!Put them in descending order of magnitude  
  CALL ReSortUgs(IEquivalentUgKey,CUniqueUg,Iuid)  
  
END SUBROUTINE StructureFactorInitialisation

SUBROUTINE Absorption (IErr)  

  USE MyNumbers
  USE WriteToScreen

  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: ind,jnd,knd,lnd,mnd,IErr,Ieval
  REAL(RKIND) :: Rintegral,RfPrime
  COMPLEX(CKIND) :: CfPrime
  CHARACTER*200 :: SPrintString
  
  CUgMatPrime = CZERO
  SELECT CASE (IAbsorbFLAG)
    CASE(1)
    !!$ Proportional
    CUgMatPrime = CUgMatNoAbs*EXP(CIMAGONE*PI/2)*(RAbsorptionPercentage/100_RKIND)
	  
    CASE(2)
    !!$ Bird & King
    IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
      PRINT*,"Starting absorptive form factor calculation..."
    END IF
    !Uses numerical integration of a function BirdKing to calculate absorptive form factor
	DO ind=2,2!nReflections!work down the first column of the Ug matrix
      IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
        PRINT*,ind,"of",nReflections
      END IF
	  RfPrime=0
      DO knd=1,1!INAtomsUnitCell
	    ICurrentZ = IAtomicNumber(knd)!Atomic number, global variable
        RCurrentB = RIsoDW(knd)!Debye-Waller constant, global variable
        RCurrentG = RgMatrixMagnitude(ind,1)!g-vector magnitude, global variable
		CALL DoubleIntegrate(Rintegral,IErr)
        IF( IErr.NE.0 ) THEN
          PRINT*,"Absorption(",my_rank,") error in Integrate"
          RETURN
        ENDIF
		RfPrime=RfPrime+Rintegral
      END DO
	  CfPrime=Rfprime*CIMAGONE
      IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
        PRINT*,"U'g=",CfPrime!CUgMatPrime(ind,1)
      END IF
   END DO

    !NB Only the lower half of the U'g matrix was calculated, this completes the upper half
    CUgMatPrime = CUgMatPrime + CONJG(TRANSPOSE(CUgMatPrime))!Need to think about this
    !Keep with proportional model while debugging
	CUgMatPrime = CUgMatNoAbs*EXP(CIMAGONE*PI/2)*(RAbsorptionPercentage/100_RKIND)
	
    CASE Default
	
  END SELECT
  CUgMat = CUgMatNoAbs+CUgMatPrime
  
  IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
   PRINT*,"Ug matrix, including absorption (nm^-2)"
	DO ind =1,20
     WRITE(SPrintString,FMT='(3(1X,I3),A1,8(1X,F6.2,F6.2))') NINT(Rhkl(ind,:)),":",100*CUgMat(ind,1:8)
     PRINT*,TRIM(SPrintString)
    END DO
  END IF	   
	   
END SUBROUTINE Absorption

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DoubleIntegrate(RResult,IErr) 

  USE MyNumbers
  USE WriteToScreen

  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr
  REAL(RKIND) :: RResult

  !Will become a 2d integral, nothing here yet until it works
  RSprimeY=0!second dimension, global variable

  CALL Integrate(RResult,IErr)

  
END SUBROUTINE DoubleIntegrate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE Integrate(RResult,IErr) 

  USE MyNumbers
  USE WriteToScreen

  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,Ieval
  INTEGER(IKIND), PARAMETER :: inf=1
  INTEGER(IKIND), PARAMETER :: limit=500
  INTEGER(IKIND), PARAMETER :: lenw= limit*4
!  REAL(RKIND), EXTERNAL :: BirdKing
  REAL(RKIND) :: BirdKing
!  REAL(RKIND), EXTERNAL :: debug
  REAL(RKIND) :: RAccuracy,RError,RResult,dd,ee
  COMPLEX(CKIND) :: CfPrime

  INTEGER(IKIND) last, iwork(limit)
  REAL(RKIND) work(lenw)

  ee=0.0
  dd=BirdKing(ee)
  PRINT*,ee,":", dd
  ee=1.0
  dd=BirdKing(ee)
  PRINT*,ee,":", dd
  ee=10.0
  dd=BirdKing(ee) 
  PRINT*,ee,":", dd
  ee=1000.0
  dd=BirdKing(ee) 
  PRINT*,ee,":", dd
 
  RAccuracy=0.1D0!accuracy of integration
  !CALL dqagi(BirdKing,ZERO,inf,0,RAccuracy,RResult,RError,Ieval,IErr,&
  !CALL dqagi(debug,ZERO,inf,ZERO,RAccuracy,RResult,RError,Ieval,IErr,&
!       limit, lenw, last, iwork, work )
!  PRINT*,RResult,RError
  
END SUBROUTINE Integrate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Defines a Bird & King integrand to calculate an absorptive scattering factor 
FUNCTION BirdKing(RSprimeX)
  !From Bird and King, Acta Cryst A46, 202 (1990)
  !ICurrentZ is atomic number, global variable
  !RCurrentB is Debye-Waller constant b=8*pi*<u^2>, where u is mean square thermal vibration amplitude in Angstroms, global variable
  !RCurrentG is magnitude of scattering vector in 1/A (NB exp(i*g.r), physics convention, global variable
  !RSprime is dummy parameter for integration [s'x s'y]
  !RSprimeY is passed as a global variable so we just have an integral in 1D from 0 to inf
  USE MyNumbers
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind
  REAL(RKIND):: BirdKing,Rs,Rg1,Rg2,RsEff,Kirkland
  REAL(RKIND), INTENT(IN) :: RSprimeX
  REAL(RKIND),DIMENSION(2) :: RGprime
  RGprime=2*TWOPI*(/RSprimeX,RSprimeY/)
  !Since [s'x s'y]  is a dummy parameter for integration I can assign s'x //g
  Rg1=SQRT( (RCurrentG/2+RGprime(1))**2 + RGprime(2)**2 )
  Rg2=SQRT( (RCurrentG/2-RGprime(1))**2 - RGprime(2)**2 )
  RsEff=RSprimeX**2+RSprimeY**2-RCurrentG**2/(16*TWOPI**2)
  BirdKing=Kirkland(Rg1)*Kirkland(Rg2)*(1-EXP(-2*RCurrentB*RsEff ) )
  PRINT*,"g'=",RGprime,":",Kirkland(Rg1),BirdKing
  
END FUNCTION BirdKing

FUNCTION debug(x)
  USE MyNumbers
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  
  IMPLICIT NONE
  REAL(RKIND):: debug,x
  debug=x/(1+x**4)

END FUNCTION debug
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Defines a Kirkland scattering factor 
FUNCTION Kirkland(Rg)
  !From Appendix C of Kirkland, "Advanced Computing in Electron Microscopy", 2nd ed.
  !ICurrentZ is atomic number, global variable
  !RCurrentG is magnitude of scattering vector in 1/A (NB exp(i*g.r), physics convention), global variable
  !Kirkland scattering factor is in Angstrom units
  !Rg is a dummy variable, just to give the function an argument, and is not used
  USE MyNumbers
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,IErr
  REAL(RKIND):: Kirkland,Ra,Rb,Rc,Rd,Rq,Rg

  !NB Kirkland scattering factors are calculated in the optics convention exp(2*pi*i*q.r)
  Rq=RCurrentG/TWOPI
  Kirkland=ZERO;
  !Equation C.15
  DO ind = 1,3
    Ra=RScattFactors(ICurrentZ,ind*2-1);
    Rb=RScattFactors(ICurrentZ,ind*2);
    Rc=RScattFactors(ICurrentZ,ind*2+5);
    Rd=RScattFactors(ICurrentZ,ind*2+6);
    Kirkland = Kirkland + Ra/((Rq**2)+Rb)+Rc*EXP(-(Rd*Rq**2));
  END DO
  
END FUNCTION Kirkland