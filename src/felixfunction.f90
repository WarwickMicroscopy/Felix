!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixrefine
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
!  This file is part of felixrefine.
!
!  felixrefine is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  felixrefine is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with felixrefine.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: Felixrefine.f90,v 1.89 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SimulateAndFit(RFigureofMerit,RIndependentVariable,Iiter,IExitFLAG,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,IExitFLAG,IThickness,ind,jnd
  REAL(RKIND),DIMENSION(INoOfVariables) :: RIndependentVariable
  REAL(RKIND) :: RFigureofMerit
  INTEGER(IKIND),INTENT(IN) :: Iiter
  COMPLEX(CKIND),DIMENSION(nReflections,nReflections) :: CUgMatDummy
  
  IF(IWriteFLAG.GE.10.AND.my_rank.EQ.0) THEN
     PRINT*,"SimulateAndFit(",my_rank,")"
  END IF

  IF (IRefineMode(1).EQ.1 .OR. IRefineMode(12).EQ.1) THEN  !Ug refinement; update structure factors 
  !Dummy Matrix to contain new iterative values
    CUgMatDummy = CZERO    !NB these are Ug's without absorption
    jnd=1
    DO ind = 1+IUgOffset,INoofUgs+IUgOffset!=== temp changes so real part only***
      IF ( (ABS(REAL(CUniqueUg(ind),RKIND)).GE.RTolerance).AND.&!===
           (ABS(AIMAG(CUniqueUg(ind))).GE.RTolerance)) THEN!===use both real and imag parts
        CUniqueUg(ind)=CMPLX(RIndependentVariable(jnd),RIndependentVariable(jnd+1))!===
        jnd=jnd+2!===
      ELSEIF ( ABS(AIMAG(CUniqueUg(ind))).LT.RTolerance ) THEN!===use only real part
        CUniqueUg(ind)=CMPLX(RIndependentVariable(jnd),ZERO)!===
!===        CUniqueUg(ind)=CMPLX(RIndependentVariable(jnd),AIMAG(CUniqueUg(ind)))!replacement line, remove to revert
        jnd=jnd+1
      ELSEIF ( ABS(REAL(CUniqueUg(ind),RKIND)).LT.RTolerance ) THEN!===use only imag part
        CUniqueUg(ind)=CMPLX(ZERO,RIndependentVariable(jnd))!===
        jnd=jnd+1!===
      ELSE!===should never happen
        PRINT*,"Warning - zero structure factor!",ind,":",CUniqueUg(IEquivalentUgKey(ind))!===
		IErr=1
      END IF!===
      WHERE(ISymmetryRelations.EQ.IEquivalentUgKey(ind))
        CUgMatDummy = CUniqueUg(ind)
      END WHERE
      WHERE(ISymmetryRelations.EQ.-IEquivalentUgKey(ind))
        CUgMatDummy = CONJG(CUniqueUg(ind))
      END WHERE
    END DO
    WHERE(ABS(CUgMatDummy).GT.TINY)
      CUgMatNoAbs = CUgMatDummy
    END WHERE
	CALL Absorption(IErr)
    IF( IErr.NE.0 ) THEN
      PRINT*,"SimulateAndFit(",my_rank,")error in Absorption"
      RETURN
    END IF	
    RAbsorptionPercentage = RIndependentVariable(jnd)!===![[[

  ELSE !everything else

	!Change variables
    CALL UpdateVariables(RIndependentVariable,IErr)
    IF( IErr.NE.0 ) THEN
      PRINT*,"SimulateAndFit(",my_rank,")error in UpdateVariables"
      RETURN
    END IF
	!recalculate unit cell
    CALL UniqueAtomPositions(IErr)
    IF( IErr.NE.0 ) THEN
      PRINT*,"SimulateAndFit(",my_rank,")error in UniqueAtomPositions"
      RETURN
    END IF

  END IF

  IF (my_rank.EQ.0) THEN
    CALL PrintVariables(IErr)
    IF( IErr.NE.0 ) THEN
      PRINT*,"SimulateAndFit(",my_rank,")error in PrintVariables"
      RETURN
    END IF
  END IF

  CALL FelixFunction(IErr) ! Simulate !!  
  IF( IErr.NE.0 ) THEN
    PRINT*,"SimulateAndFit(",my_rank,")error in FelixFunction"
    RETURN
  END IF

  IF(my_rank.EQ.0) THEN   
    CALL CreateImagesAndWriteOutput(Iiter,IExitFLAG,IErr) 
    IF( IErr.NE.0 ) THEN
      PRINT*,"SimulateAndFit(",my_rank,")error in CreateImagesAndWriteOutput"
      RETURN
    ENDIF
    !This is the key parameter!!!****     
    RFigureofMerit = RCrossCorrelation
  END IF
  CALL MPI_BCAST(RFigureofMerit,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)

END SUBROUTINE SimulateAndFit

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE FelixFunction(IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  !--------------------------------------------------------------------
  ! local variable definitions
  !--------------------------------------------------------------------
  
  INTEGER(IKIND) :: IErr,ind,jnd,knd,pnd,IThicknessIndex,IIterationFLAG
  INTEGER(IKIND) :: IAbsorbTag = 0
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RFinalMontageImageRoot

  IF(IWriteFLAG.GE.10.AND.my_rank.EQ.0) THEN
     PRINT*,"Felixfunction(", my_rank, "): starting the eigenvalue problem"
     PRINT*,"Felixfunction(",my_rank,")pixels",ILocalPixelCountMin," to ",ILocalPixelCountMax
  END IF

  !Reset simuation--------------------------------------------------------------------  
  RIndividualReflections = ZERO
  !Update scattering matrix, if it's not a Ug refinement and it's not the baseline simulation-------------------------------------  
  IF (IRefineMode(1).NE.1 .AND. IRefineMode(12).NE.1 .AND. IInitialSimulationFLAG.NE.1) THEN
    CALL StructureFactorInitialisation(IErr)
    CALL Absorption (IErr)
    IF( IErr.NE.0 ) THEN
      PRINT*,"felixfunction(",my_rank,")error in StructureFactorInitialisation"
      RETURN
    END IF
  END IF
  IMAXCBuffer = 200000!RB what are these?
  IPixelComputed= 0

!Simulation (different local pixels for each core)--------------------------------------------------------------------  
  IF(IWriteFLAG.GE.0.AND.my_rank.EQ.0) THEN
    PRINT*,"Bloch wave calculation..."
  END IF
  DO knd = ILocalPixelCountMin,ILocalPixelCountMax,1
    jnd = IPixelLocations(knd,1)
    ind = IPixelLocations(knd,2)
    CALL BlochCoefficientCalculation(ind,jnd,knd,ILocalPixelCountMin,IErr)
    IF( IErr.NE.0 ) THEN
      PRINT*,"Felixfunction(",my_rank,") error in BlochCofficientCalculation"
      RETURN
    END IF
  END DO
 
  !MPI gatherv into RSimulatedPatterns--------------------------------------------------------------------  
  CALL MPI_GATHERV(RIndividualReflections,SIZE(RIndividualReflections),MPI_DOUBLE_PRECISION,&
                   RSimulatedPatterns,ICount,IDisplacements,MPI_DOUBLE_PRECISION,&
				   root,MPI_COMM_WORLD,IErr)
  IF( IErr.NE.0 ) THEN
    PRINT*,"Felixfunction(",my_rank,")error",IErr,"in MPI_GATHERV"
    RETURN
  END IF
  !We have done at least one simulation now
  IInitialSimulationFLAG=0

END SUBROUTINE FelixFunction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CalculateFigureofMeritandDetermineThickness(IThicknessCountFinal,IErr)
  !NB core 0 only
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: ind,jnd,knd,IErr,ICountedPixels,IThickness,hnd
  INTEGER(IKIND),DIMENSION(INoOfLacbedPatterns) :: IThicknessByReflection
  INTEGER(IKIND),INTENT(OUT) :: IThicknessCountFinal
  REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: RSimulatedImage,RExperimentalImage
  REAL(RKIND) :: RCrossCorrelationOld,RIndependentCrossCorrelation,RThickness,&
       PhaseCorrelate,Normalised2DCrossCorrelation,ResidualSumofSquares,RThicknessRange
  REAL(RKIND),DIMENSION(INoOfLacbedPatterns) :: RReflectionCrossCorrelations,RReflectionThickness
  CHARACTER*200 :: SPrintString
       
  
  IF (IWriteFLAG.GE.10) THEN
     PRINT*,"CalculateFigureofMeritandDetermineThickness(",my_rank,")"
  END IF

  RReflectionCrossCorrelations = ZERO

  DO hnd = 1,INoOfLacbedPatterns
    RCrossCorrelationOld = 1.0E15 !A large Number
    RThickness = ZERO
    IThicknessByReflection(hnd) = 1!default value if no correlation
    DO ind = 1,IThicknessCount
      ICountedPixels = 0
      RSimulatedImage = ZERO
	  !put 1D array RSimulatedPatterns into 2D image RSimulatedImage
      !remember dimensions of RSimulatedPatterns(INoOfLacbedPatterns,IThicknessCount,IPixelTotal)
      DO jnd = 1,2*IPixelCount
        DO knd = 1,2*IPixelCount
          ICountedPixels = ICountedPixels+1
          RSimulatedImage(jnd,knd) = RSimulatedPatterns(hnd,ind,ICountedPixels)
        END DO
      END DO
               
      SELECT CASE (IImageProcessingFLAG)

      CASE(0)!no processing
        RExperimentalImage = RImageExpi(:,:,hnd)
           
      CASE(1)!square root before perfoming corration
        RSimulatedImage=SQRT(RSimulatedImage)
        RExperimentalImage =  SQRT(RImageExpi(:,:,hnd))
           
      CASE(2)!log before performing correlation
        WHERE (RSimulatedImage.GT.TINY**2)
          RSimulatedImage=LOG(RSimulatedImage)
        ELSEWHERE
          RSimulatedImage = TINY**2
        END WHERE
        WHERE (RExperimentalImage.GT.TINY**2)
          RExperimentalImage = LOG(RImageExpi(:,:,hnd))
        ELSEWHERE
          RExperimentalImage =  TINY**2
        END WHERE
              
      CASE(4)!Apply gaussian blur to simulated image
        RExperimentalImage = RImageExpi(:,:,hnd)
       ! IF(my_rank.EQ.0) THEN
       !   PRINT*,"Gaussian blur radius =",RBlurRadius
       ! END IF
        CALL BlurG(RSimulatedImage,IErr)
		
      END SELECT

      RIndependentCrossCorrelation = ZERO   
      SELECT CASE (ICorrelationFLAG)
         
      CASE(0) ! Phase Correlation
        RIndependentCrossCorrelation=ONE-& ! So Perfect Correlation = 0 not 1
           PhaseCorrelate(RSimulatedImage,RExperimentalImage,&
           IErr,2*IPixelCount,2*IPixelCount)
           
      CASE(1) ! Residual Sum of Squares (Non functional) 
        RIndependentCrossCorrelation = ResidualSumofSquares(&
                RSimulatedImage,RImageExpi(:,:,hnd),IErr)
           
      CASE(2) ! Normalised Cross Correlation
 		   RIndependentCrossCorrelation = ONE-& ! So Perfect Correlation = 0 not 1
           Normalised2DCrossCorrelation(RSimulatedImage,RExperimentalImage,IErr)

      END SELECT
                
      IF(ABS(RIndependentCrossCorrelation).LT.RCrossCorrelationOld) THEN 
        RCrossCorrelationOld = RIndependentCrossCorrelation
        IThicknessByReflection(hnd) = ind
        RReflectionThickness(hnd) = RInitialThickness +&
        IThicknessByReflection(hnd)*RDeltaThickness
      END IF
    END DO
    RReflectionCrossCorrelations(hnd) = RCrossCorrelationOld
  END DO

  !RB assume that the best thickness is given by the mean of individual thicknesses  
  IF (ISimFLAG.EQ.0) THEN !Refine Mode (Output Single Image)
     IThicknessCountFinal = SUM(IThicknessByReflection)/INoOfLacbedPatterns
     RThickness = RInitialThickness + (IThicknessCountFinal-1)*RDeltaThickness
     RThicknessRange=( MAXVAL(IThicknessByReflection)-MINVAL(IThicknessByReflection) )*&
          RDeltaThickness
     RCrossCorrelation = SUM(RReflectionCrossCorrelations*RWeightingCoefficients)/&
          REAL(INoOfLacbedPatterns,RKIND)
   
  ELSE !Sim Mode (Output all thicknesses)
     IThicknessCountFinal = INoOfLacbedPatterns
     RThickness = RInitialThickness + (IThicknessCountFinal-1)*RDeltaThickness
     RThicknessRange=( MAXVAL(IThicknessByReflection)-MINVAL(IThicknessByReflection) )*&
          RDeltaThickness
  END IF

  IF(my_rank.eq.0) THEN
    WRITE(SPrintString,FMT='(A18,I4,A10)') "Specimen thickness ",NINT(RThickness)," Angstroms"
    PRINT*,TRIM(ADJUSTL(SPrintString))
    WRITE(SPrintString,FMT='(A15,I4,A10)') "Thickness range",NINT(RThicknessRange)," Angstroms"
    PRINT*,TRIM(ADJUSTL(SPrintString))
  END IF

END SUBROUTINE CalculateFigureofMeritandDetermineThickness

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CreateImagesAndWriteOutput(Iiter,IExitFLAG,IErr)

!NB core 0 only
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,IThicknessIndex,IThicknessCountFinal,Iiter,IExitFLAG

  
  CALL CalculateFigureofMeritandDetermineThickness(IThicknessCountFinal,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateImagesAndWriteOutput(", my_rank, ") error ", IErr, &
          "Calling function CalculateFigureofMeritandDetermineThickness"
     RETURN
  ENDIF

!!$     OUTPUT -------------------------------------  
  IF (ISimFLAG.EQ.0) THEN !Refine Mode, Only one Thickness output
     IThicknessIndex=IThicknessCountFinal
     CALL WriteIterationOutput(Iiter,IThicknessIndex,IExitFLAG,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"CreateImagesAndWriteOutput(",my_rank,")error in WriteIterationOutput"
        RETURN
     ENDIF
          
  ELSE !Sim mode - All Thicknesses Output, IThicknessCountFinal equals NoofLacbed patterns
     DO IThicknessIndex = 1,IThicknessCountFinal
        CALL WriteIterationOutput(Iiter,IThicknessIndex,IExitFLAG,IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"CreateImagesAndWriteOutput(",my_rank,")error in WriteIterationOutput"
           RETURN
        ENDIF
     END DO
  END IF
!!$           ----------------------------------------
END SUBROUTINE CreateImagesAndWriteOutput
   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE UpdateVariables(RIndependentVariable,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IVariableType,IVectorID,IAtomID,IErr,ind
  REAL(RKIND),DIMENSION(INoOfVariables),INTENT(IN) :: RIndependentVariable

  !!$  Fill the Independent Value array with values
  IF(IRefineMode(2).EQ.1) THEN     
     RBasisAtomPosition = RInitialAtomPosition!RB is this redundant? 
  END IF

  DO ind = 1,INoOfVariables
    IVariableType = IIterativeVariableUniqueIDs(ind,2)
    SELECT CASE (IVariableType)
    CASE(1)
      !RB structure factor refinement, do in UpdateStructureFactors
    CASE(2)
      !CALL ConvertVectorMovementsIntoAtomicCoordinates(ind,RIndependentVariable,IErr)
	  !The vector being used
	  IVectorID = IIterativeVariableUniqueIDs(ind,3)
	  !The atom being moved
      IAtomID = IAllowedVectorIDs(IVectorID)
	  !Change in position
      RBasisAtomPosition(IAtomID,:) = RBasisAtomPosition(IAtomID,:) + &
         RIndependentVariable(ind)*RAllowedVectors(IVectorID,:)
     CASE(3)
        RBasisOccupancy(IIterativeVariableUniqueIDs(ind,3)) = &
             RIndependentVariable(ind)
     CASE(4)
        RBasisIsoDW(IIterativeVariableUniqueIDs(ind,3)) = &
             RIndependentVariable(ind)
     CASE(5)
        RAnisotropicDebyeWallerFactorTensor(&
             IIterativeVariableUniqueIDs(ind,3),&
             IIterativeVariableUniqueIDs(ind,4),&
             IIterativeVariableUniqueIDs(ind,5)) = & 
             RIndependentVariable(ind)
     CASE(6)
        SELECT CASE(IIterativeVariableUniqueIDs(ind,3))
        CASE(1)
           RLengthX = RIndependentVariable(ind)
        CASE(2)
           RLengthY = RIndependentVariable(ind)
        CASE(3)
           RLengthZ = RIndependentVariable(ind)
        END SELECT
     CASE(7)
        SELECT CASE(IIterativeVariableUniqueIDs(ind,3))
        CASE(1)
           RAlpha = RIndependentVariable(ind)
        CASE(2)
           RBeta = RIndependentVariable(ind)
        CASE(3)
           RGamma = RIndependentVariable(ind)
        END SELECT
     CASE(8)
        RConvergenceAngle = RIndependentVariable(ind)
     CASE(9)
        RAbsorptionPercentage = RIndependentVariable(ind)
     CASE(10)
        RAcceleratingVoltage = RIndependentVariable(ind)
     CASE(11)
        RRSoSScalingFactor = RIndependentVariable(ind)
     END SELECT
  END DO

 END SUBROUTINE UpdateVariables
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE PrintVariables(IErr)

  USE WriteToScreen
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,IVariableType,jnd,knd
  REAL(RKIND),DIMENSION(3) :: RCrystalVector
  CHARACTER*200 :: SPrintString

  RCrystalVector = [RLengthX,RLengthY,RLengthZ]

  DO ind = 1,IRefinementVariableTypes
     IF (IRefineMode(ind).EQ.1) THEN
        SELECT CASE(ind)
        CASE(1)
           WRITE(SPrintString,FMT='(A18,1X,F5.2)') "Current Absorption",RAbsorptionPercentage
           PRINT*,TRIM(ADJUSTL(SPrintString))
           PRINT*,"Current Structure Factors nm^-2: amplitude, phase (deg)"!RB should also put in hkl here
           DO jnd = 1+IUgOffset,INoofUgs+IUgOffset
              WRITE(SPrintString,FMT='(2(1X,F7.3),2X,A1,1X,F6.3,1X,F6.2)') 100*CUniqueUg(jnd),":",&
              ABS(CUniqueUg(jnd)),180*ATAN2(AIMAG(CUniqueUg(jnd)),REAL(CUniqueUg(jnd)))/PI
              PRINT*,TRIM(ADJUSTL(SPrintString))
           END DO           

           CASE(2)
           PRINT*,"Current Atomic Coordinates"
           DO jnd = 1,SIZE(RBasisAtomPosition,DIM=1)
              WRITE(SPrintString,FMT='(A2,3(1X,F9.4))') SBasisAtomName(jnd),RBasisAtomPosition(jnd,:)
              PRINT*,TRIM(ADJUSTL(SPrintString))              
           END DO

           CASE(3)
           PRINT*,"Current Atomic Occupancy"
           DO jnd = 1,SIZE(RBasisOccupancy,DIM=1)
              WRITE(SPrintString,FMT='(A2,1X,F9.6)') SBasisAtomName(jnd),RBasisOccupancy(jnd)
              PRINT*,TRIM(ADJUSTL(SPrintString))
           END DO

           CASE(4)
           PRINT*,"Current Isotropic Debye Waller Factors"
           DO jnd = 1,SIZE(RBasisIsoDW,DIM=1)
              WRITE(SPrintString,FMT='(A2,1X,F9.6)') SBasisAtomName(jnd),RBasisIsoDW(jnd)
              PRINT*,TRIM(ADJUSTL(SPrintString))
           END DO

           CASE(5)
           PRINT*,"Current Anisotropic Debye Waller Factors"
           DO jnd = 1,SIZE(RAnisotropicDebyeWallerFactorTensor,DIM=1)
              DO knd = 1,3
                 WRITE(SPrintString,FMT='(A2,3(1X,F9.4))') SBasisAtomName(jnd),RAnisotropicDebyeWallerFactorTensor(jnd,knd,:)
              PRINT*,TRIM(ADJUSTL(SPrintString))
              END DO
           END DO

           CASE(6)
           PRINT*,"Current Unit Cell Parameters"
           WRITE(SPrintString,FMT='(3(1X,F9.6))') RLengthX,RLengthY,RLengthZ
              PRINT*,TRIM(ADJUSTL(SPrintString))

          CASE(7)
           PRINT*,"Current Unit Cell Angles"
           WRITE(SPrintString,FMT='(3(F9.6,1X))') RAlpha,RBeta,RGamma
              PRINT*,TRIM(ADJUSTL(SPrintString))

          CASE(8)
           PRINT*,"Current Convergence Angle"
           WRITE(SPrintString,FMT='((F9.6,1X))') RConvergenceAngle
              PRINT*,TRIM(ADJUSTL(SPrintString))

          CASE(9)
           PRINT*,"Current Absorption Percentage"
           WRITE(SPrintString,FMT='((F9.6,1X))') RAbsorptionPercentage
              PRINT*,TRIM(ADJUSTL(SPrintString))

          CASE(10)
           PRINT*,"Current Accelerating Voltage"
           WRITE(SPrintString,FMT='((F9.6,1X))') RAcceleratingVoltage
              PRINT*,TRIM(ADJUSTL(SPrintString))

          CASE(11)
           PRINT*,"Current Residual Sum of Squares Scaling Factor"
           WRITE(SPrintString,FMT='((F9.6,1X))') RRSoSScalingFactor
              PRINT*,TRIM(ADJUSTL(SPrintString))
        END SELECT
     END IF
  END DO

END SUBROUTINE PrintVariables

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE UpdateStructureFactors(RIndependentVariable,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,ind,jnd
  REAL(RKIND),DIMENSION(INoOfVariables),INTENT(IN) :: RIndependentVariable
  CHARACTER*200 :: SPrintString
    
!NB these are Ug's without absorption
  jnd=1
  DO ind = 1+IUgOffset,INoofUgs+IUgOffset!=== temp changes so real part only***
    IF ( (ABS(REAL(CUniqueUg(ind),RKIND)).GE.RTolerance).AND.&!===
       (ABS(AIMAG(CUniqueUg(ind))).GE.RTolerance)) THEN!use both real and imag parts!===
      CUniqueUg(ind)=CMPLX(RIndependentVariable(jnd),RIndependentVariable(jnd+1))!===
      jnd=jnd+2!===
    ELSEIF ( ABS(AIMAG(CUniqueUg(ind))).LT.RTolerance ) THEN!use only real part!===
      CUniqueUg(ind)=CMPLX(RIndependentVariable(jnd),ZERO)!===
!===      CUniqueUg(ind)=CMPLX(RIndependentVariable(jnd),AIMAG(CUniqueUg(ind)))!===replacement line, delete to revert
      jnd=jnd+1
    ELSEIF ( ABS(REAL(CUniqueUg(ind),RKIND)).LT.RTolerance ) THEN!===use only imag part
      CUniqueUg(ind)=CMPLX(ZERO,RIndependentVariable(jnd))!===
      jnd=jnd+1!===
    ELSE!should never happen!===
      PRINT*,"Warning - zero structure factor!",ind,":",CUniqueUg(IEquivalentUgKey(ind))!===
    END IF!===
  END DO
  RAbsorptionPercentage = RIndependentVariable(jnd)!===
  
END SUBROUTINE UpdateStructureFactors

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ConvertVectorMovementsIntoAtomicCoordinates(IVariableID,RIndependentVariable,IErr)
!RB this is now redundant, moved up to Update Variables
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,jnd,IVariableID,IVectorID,IAtomID
  REAL(RKIND),DIMENSION(INoOfVariables),INTENT(IN) :: RIndependentVariable

!!$  Use IVariableID to determine which vector is being applied (IVectorID)
  IVectorID = IIterativeVariableUniqueIDs(IVariableID,3)
!!$  Use IVectorID to determine which atomic coordinate the vector is to be applied to (IAtomID)
  IAtomID = IAllowedVectorIDs(IVectorID)
!!$  Use IAtomID to applied the IVectodID Vector to the IAtomID atomic coordinate
  RBasisAtomPosition(IAtomID,:) = RBasisAtomPosition(IAtomID,:) + &
       RIndependentVariable(IVariableID)*RAllowedVectors(IVectorID,:)
  
END SUBROUTINE ConvertVectorMovementsIntoAtomicCoordinates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE BlurG(RImageToBlur,IErr)
  !performs a 2D Gaussian blur on the input image
  !renormalises the output image to have the same min and max as the input image
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,jnd,IKernelRadius,IKernelSize
  REAL(RKIND),DIMENSION(:), ALLOCATABLE :: RGauss1D
  REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: RImageToBlur,RTempImage,RShiftImage
  REAL(RKIND) :: Rind,Rsum,Rmin,Rmax
  
  !get min and max of input image
  Rmin=MINVAL(RImageToBlur)
  Rmax=MAXVAL(RImageToBlur)

  !set up a 1D kernel of appropriate size  
  IKernelRadius=NINT(3*RBlurRadius)
  ALLOCATE(RGauss1D(2*IKernelRadius+1),STAT=IErr)!ffs
  Rsum=0
  DO ind=-IKernelRadius,IKernelRadius
    Rind=REAL(ind)
    RGauss1D(ind+IKernelRadius+1)=EXP(-(Rind**2)/((2*RBlurRadius)**2))
    Rsum=Rsum+RGauss1D(ind+IKernelRadius+1)
  END DO
  RGauss1D=RGauss1D/Rsum!normalise
  RTempImage=RImageToBlur*0_RKIND;!reset the temp image
  
  !apply the kernel in direction 1
  DO ind = -IKernelRadius,IKernelRadius
    IF (ind.LT.0) THEN
      RShiftImage(1:2*IPixelCount+ind,:)=RImageToBlur(1-ind:2*IPixelCount,:)
      DO jnd = 1,1-ind!edge fill on right
        RShiftImage(2*IPixelCount-jnd+1,:)=RImageToBlur(2*IPixelCount,:)
      END DO
    ELSE
      RShiftImage(1+ind:2*IPixelCount,:)=RImageToBlur(1:2*IPixelCount-ind,:)
      DO jnd = 1,1+ind!edge fill on left
        RShiftImage(jnd,:)=RImageToBlur(1,:)
      END DO
    END IF
    RTempImage=RTempImage+RShiftImage*RGauss1D(ind+IKernelRadius+1)
  END DO
  
  !make the 1D blurred image the input for the next direction
  RImageToBlur=RTempImage;
  RTempImage=RImageToBlur*0_RKIND;!reset the temp image

  !apply the kernel in direction 2  
  DO ind = -IKernelRadius,IKernelRadius
    IF (ind.LT.0) THEN
      RShiftImage(:,1:2*IPixelCount+ind)=RImageToBlur(:,1-ind:2*IPixelCount)
      DO jnd = 1,1-ind!edge fill on bottom
        RShiftImage(:,2*IPixelCount-jnd+1)=RImageToBlur(:,2*IPixelCount)
      END DO
    ELSE
      RShiftImage(:,1+ind:2*IPixelCount)=RImageToBlur(:,1:2*IPixelCount-ind)
      DO jnd = 1,1+ind!edge fill on top
        RShiftImage(:,jnd)=RImageToBlur(:,1)
      END DO
    END IF
    RTempImage=RTempImage+RShiftImage*RGauss1D(ind+IKernelRadius+1)
  END DO
  DEALLOCATE(RGauss1D,STAT=IErr)

  !set intensity range of outpt image to match that of the input image
  RTempImage=RTempImage-MINVAL(RTempImage)
  RTempImage=RTempImage*(Rmax-Rmin)/MAXVAL(RTempImage)+Rmin
  !return the blurred image
  RImageToBlur=RTempImage;
  
END SUBROUTINE BlurG

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REAL(RKIND) FUNCTION RStandardError(RStandardDeviation,RMean,RFigureofMerit,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr
  REAL(RKIND),INTENT(INOUT) :: RStandardDeviation,RMean
  REAL(RKIND),INTENT(IN) :: RFigureofMerit
  
  IF (IStandardDeviationCalls.GT.1) THEN
     RMean = (RMean*REAL(IStandardDeviationCalls,RKIND) + &
          RFigureofMerit)/REAL(IStandardDeviationCalls+1,RKIND)
     RStandardDeviation = SQRT(&
          ((REAL(IStandardDeviationCalls,RKIND)*RStandardDeviation**2)+&
          (RFigureofMerit-RMean)**2)/ &
          REAL(IStandardDeviationCalls+1,RKIND))   
  ELSE
     RMean = RFigureofMerit
     RStandardDeviation = ZERO
  END IF
  RStandardError = RStandardDeviation/SQRT(REAL(IStandardDeviationCalls+1,RKIND))
  IStandardDeviationCalls = IStandardDeviationCalls + 1
     
END FUNCTION  RStandardError
