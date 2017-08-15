!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-17, all rights reserved
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
!  Felix is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  Felix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with Felix.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: Felixrefine.f90,v 1.89 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! All procedures conatained in this file:
! SimulateAndFit( )                               ---
! FelixFunction( )                                ---
! CalculateFigureofMeritandDetermineThickness( )  ---
! UpdateVariables( )                              ---
! PrintVariables( )                               ---
! BlurG( )                                        ---


!?? uses UniqueAtomPositions from crystallography.f90
!?? uses WriteOutVariables from RefineWriteOut.f90

!>
!! Module-description: Holds main top level simulating subroutines considering
!! different thicknesses
!!
MODULE felixfunction_mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: SimulateAndFit, FelixFunction, CalculateFigureofMeritandDetermineThickness
  !?? JR CalculateFigureofMeritandDetermineThickness could be made private with some tweaks

  CONTAINS

  !>
  !! Procedure-description: Simulate and fit
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!  
  SUBROUTINE SimulateAndFit(RIndependentVariable,Iter,IExitFLAG,IErr)

    ! JR non-Ug case: RIndependentVariable changeble ones, UpdateVariables() updates
    ! corresponding global variables (non changeable stay as input from .inp)
    ! UniqueAtomPositions() - recalculate unit cell ?
    ! CUgMat = CUgMatNoAbs + CUgMatPrime
    ! simulate, felixfunction - CUgMat + others ---> RImageSimi
    ! CalculateFigureofMeritandDetermineThickness() - im process, compares sim/expi images
    ! ---> RFigureofMerit
    ! RImageExpi(x,y,LACBED_ID), RImageSimi(x,y,LACBED_ID, thickness_ID)
    ! MPI_BCAST(RFigureofMerit) - send to all cores ? parallel    

    !?? uses UniqueAtomPositions from crystallography.f90 every (non-Ug) iteration
    !?? uses WriteOutVariables from RefineWriteOut.f90

    USE MyNumbers
    USE message_mod; USE alert_mod

    USE MyMPI
    USE Ug_mod
    USE crystallography_mod

    !?? global inputs and outputs of subroutines using global
    ! global inputs
    USE IPARA, ONLY : INoOfVariables, nReflections, IAbsorbFLAG, IMethodFLAG, INoofUgs, &
          IPixelCount, IPrint, ISimFLAG, ISymmetryRelations, IUgOffset, IRefineMode, &
          IEquivalentUgKey
    USE RPARA, ONLY : RConvergenceAngle, RMinimumGMag, RTolerance, RRelativisticCorrection, &
          RVolume, RgMatrix, RgMatrixMagnitude, RCurrentGMagnitude

    ! gloabl outputs
    USE IPARA, ONLY : IPreviousPrintedIteration 
    USE CPARA, ONLY : CUgMatNoAbs, CUniqueUg
    USE RPARA, ONLY : RAbsorptionPercentage, RDeltaK, RFigureofMerit, RSimulatedPatterns

    IMPLICIT NONE

    REAL(RKIND),INTENT(INOUT) :: RIndependentVariable(INoOfVariables)
    INTEGER(IKIND),INTENT(INOUT) :: Iter
    INTEGER(IKIND),INTENT(OUT) :: IErr
    INTEGER(IKIND) :: ind,jnd, IExitFLAG,IThicknessIndex, ILoc(2), IUniqueUgs
    REAL(RKIND) :: RCurrentG(3), RScatteringFactor
    COMPLEX(CKIND) :: CUgMatDummy(nReflections,nReflections),CVgij
    CHARACTER*200 :: SFormat,SPrintString

    CALL message(LM,"Iteration ",Iter)
    
    !\/----------------------------------------------------------------------
    IF (IRefineMode(1).EQ.1) THEN  ! Ug refinement; update structure factors 
      ! Dummy Matrix to contain new iterative values
      CUgMatDummy = CZERO    ! NB these are Ug's without absorption
      jnd=1
      ! work through the Ug's to update
      DO ind = 1+IUgOffset,INoofUgs+IUgOffset
      ! Don't update components smaller than RTolerance:
      ! 3 possible types of Ug, complex, real and imaginary
      IF ( (ABS(REAL(CUniqueUg(ind),RKIND)).GE.RTolerance).AND.&
            (ABS(AIMAG(CUniqueUg(ind))).GE.RTolerance)) THEN ! use both real and imag parts
        CUniqueUg(ind)=CMPLX(RIndependentVariable(jnd),RIndependentVariable(jnd+1))
        jnd=jnd+2
      ELSEIF ( ABS(AIMAG(CUniqueUg(ind))).LT.RTolerance ) THEN ! use only real part
        CUniqueUg(ind)=CMPLX(RIndependentVariable(jnd),ZERO)
        jnd=jnd+1
      ELSEIF ( ABS(REAL(CUniqueUg(ind),RKIND)).LT.RTolerance ) THEN ! use only imag part
        CUniqueUg(ind)=CMPLX(ZERO,RIndependentVariable(jnd))
        jnd=jnd+1
      ELSE ! should never happen
        CALL message(LS,"Warning - zero structure factor!")
        CALL message(LS,no_tag,"A CUniqueUg element has value = ",&
              CUniqueUg(IEquivalentUgKey(ind)))
        CALL message(LS,"element number = ",ind)
        IErr=1
        IF(l_alert(IErr,"SimulateAndFit"," to Check Ug type")) RETURN
      END IF

      ! Update the Ug matrix for this Ug
      WHERE(ISymmetryRelations.EQ.IEquivalentUgKey(ind))
         CUgMatDummy = CUniqueUg(ind)
      END WHERE
      WHERE(ISymmetryRelations.EQ.-IEquivalentUgKey(ind))
         CUgMatDummy = CONJG(CUniqueUg(ind))
      END WHERE
      END DO
      ! put the changes into CUgMatNoAbs
      WHERE(ABS(CUgMatDummy).GT.TINY)
        CUgMatNoAbs = CUgMatDummy
      END WHERE
      CALL Absorption(IErr)
      IF(l_alert(IErr,"SimulateAndFit()","Absorption()")) RETURN
      IF (IAbsorbFLAG.EQ.1) THEN ! proportional absorption
        RAbsorptionPercentage = RIndependentVariable(jnd)
      END IF

    ELSE ! not Ug refinement
      ! Update variables
      CALL UpdateVariables(RIndependentVariable,IErr)
      IF(l_alert(IErr,"SimulateAndFit()","UpdateVariables()")) RETURN
      IF (IRefineMode(8).EQ.1) THEN ! convergence angle
        ! recalculate resolution in k space
        RDeltaK = RMinimumGMag*RConvergenceAngle/REAL(IPixelCount,RKIND) 
        !?? in the past wrote following to iterationlog.txt
        !?? Iter,RFigureofMerit,RConvergenceAngle
      ELSE
        ! basis has changed in some way, recalculate unit cell
        CALL UniqueAtomPositions(IErr)
        IF(l_alert(IErr,"SimulateAndFit()","UniqueAtomPositions()")) RETURN
      END IF

      !--------------------------------------------------------------------
      ! update scattering matrix Ug
      !--------------------------------------------------------------------

      ! calculate CUgMatNoAbs
      CUgMatNoAbs = CZERO
      ! Work through unique Ug's
      IUniqueUgs = SIZE(IEquivalentUgKey)
      DO ind = 1, IUniqueUgs
        ! number of this Ug
        jnd = IEquivalentUgKey(ind)
        ! find the position of this Ug in the matrix
        ILoc = MINLOC(ABS(ISymmetryRelations-jnd))
        RCurrentG = RgMatrix(ILoc(1),ILoc(2),:) ! g-vector, local variable
        RCurrentGMagnitude = RgMatrixMagnitude(ILoc(1),ILoc(2)) ! g-vector magnitude

        CALL GetVgContributionij(RScatteringFactor,ILoc(1),ILoc(2),CVgij,IErr)
        IF(l_alert(IErr,"SimulateAndFit()","GetVgContributionij()")) RETURN

        ! now replace the values in the Ug matrix
        WHERE(ISymmetryRelations.EQ.jnd)
          CUgMatNoAbs=CVgij
        END WHERE
        ! NB for imaginary potential U(g)=U(-g)*
        WHERE(ISymmetryRelations.EQ.-jnd)
          CUgMatNoAbs = CONJG(CVgij)
        END WHERE
      END DO
      CUgMatNoAbs = CUgMatNoAbs*RRelativisticCorrection/(PI*RVolume)
      DO ind=1,nReflections ! zero diagonal
         CUgMatNoAbs(ind,ind) = ZERO
      ENDDO

      ! calculate CUgMat
      ! calculate CUgMatPrime, then CUgMat = CUgMatNoAbs + CUgMatPrime
      CALL Absorption(IErr)
      IF(l_alert(IErr,"SimulateAndFit()","Absorption()")) RETURN

    END IF
    !/\----------------------------------------------------------------------

    IF (my_rank.EQ.0) THEN ! send current values to screen
      CALL PrintVariables(IErr)
      IF(l_alert(IErr,"SimulateAndFit()","PrintVariables()")) RETURN
    END IF

    ! simulate

    RSimulatedPatterns = ZERO ! Reset simulation
    CALL FelixFunction(IErr) ! simulate 
    IF(l_alert(IErr,"SimulateAndFit()","FelixFunction()")) RETURN
    IF (IMethodFLAG.NE.1) Iter=Iter+1 ! iterate for methods other than simplex

    IF(my_rank.EQ.0) THEN
      ! Only calculate figure of merit if we are refining
      IF (ISimFLAG.EQ.0) THEN
        CALL CalculateFigureofMeritandDetermineThickness(Iter,IThicknessIndex,IErr)
        IF(l_alert(IErr,"SimulateAndFit()",&
              "CalculateFigureofMeritandDetermineThickness()")) RETURN
      END IF
      ! Write current variable list and fit to IterationLog.txt
      CALL WriteOutVariables(Iter,IErr)
      IF(l_alert(IErr,"SimulateAndFit()","WriteOutVariables()")) RETURN
      ! write images to disk every IPrint iterations, or when finished
      IF(IExitFLAG.EQ.1.OR.(Iter.GE.(IPreviousPrintedIteration+IPrint))) THEN
        CALL WriteIterationOutput(Iter,IThicknessIndex,IExitFLAG,IErr)
        IF(l_alert(IErr,"SimulateAndFit()","WriteIterationOutput()")) RETURN
        IPreviousPrintedIteration = Iter ! reset iteration counter
      END IF
    END IF

    !===================================== ! Send the fit index to all cores
    CALL MPI_BCAST(RFigureofMerit,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
    !=====================================

  END SUBROUTINE SimulateAndFit

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  !>
  !! Procedure-description: Felixfunction simulates and produces images for each 
  !! thickness
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!  
  SUBROUTINE FelixFunction(IErr)

    ! CUgMat and smaller inputs (using bloch pixel,thickness parallels)
    ! ------> RImageSimi(x, y, LACBED_pattern_ID , thickness_ID )

    USE MyNumbers
    USE IConst, ONLY : ITHREE
    USE MyMPI
    USE message_mod; USE alert_mod

    USE bloch_mod

    ! globals - bloch inputs  
    USE CPara, ONLY : CUgMat ! from Absorption
    USE RPara, ONLY : &
      RDeltaK,& ! resolution in k space, used in calculations
      RDeltaThickness,& ! input from felix.inp
      RInitialThickness,& ! input from felix.inp
      RNormDirM,RgDotNorm,RgPool,RgPoolMag,&
      Rhkl ! from HKLMake ..fill reflection pool
    USE IPara, ONLY : &
      IHKLSelectFLAG,&
      IHolzFLAG,& ! higher order Laue
      IImageFLAG,&
      IMinStrongBeams,IMinWeakBeams,&
      INoOfLacbedPatterns,&
      IPixelCount,&
      IThicknessCount,&
      nReflections,&
      IOutputReflections
    USE BlochPara, ONLY : RBigK ! from StructureFactorInitialisation            
    
    ! globals - bloch outputs
    USE RPara, ONLY : RIndividualReflections 
    USE CPara, ONLY : CAmplitudeandPhase
    ! /\ matrix of wavefunction vectors. A vector for each thickness & pixel
    ! one pixel and some parallelises in fexlixfunction join later
    USE IPara, ONLY : IPixelComputed

    !global outputs
    USE RPara, ONLY : RImageSimi,RSimulatedPatterns !?? what is the difference
    USE IPara, ONLY : IInitialSimulationFLAG

    !global inputs
    USE RPARA, ONLY : RBlurRadius
    USE IPARA, ONLY : ICount,IDisplacements,ILocalPixelCountMax,&
          ILocalPixelCountMin,IPixelLocations

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr,ind,jnd,knd,pnd,IIterationFLAG
    INTEGER(IKIND) :: IAbsorbTag = 0
    REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RFinalMontageImageRoot
    REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: RTempImage 

    ! Reset simuation   
    RIndividualReflections = ZERO
    IPixelComputed= 0!RB what is this?

    ! Simulation (different local pixels for each core)
    CALL message(LM,"Bloch wave calculation...")
    DO knd = ILocalPixelCountMin,ILocalPixelCountMax,1
      jnd = IPixelLocations(knd,1)
      ind = IPixelLocations(knd,2)
      !?? fills array for each pixel number not x & y coordinates
      CALL BlochCoefficientCalculation(ind,jnd,knd,ILocalPixelCountMin,IErr)
      IF(l_alert(IErr,"FelixFunction","BlochCoefficientCalculation")) RETURN
    END DO

    !===================================== ! MPI gatherv into RSimulatedPatterns
    CALL MPI_GATHERV(RIndividualReflections,SIZE(RIndividualReflections),MPI_DOUBLE_PRECISION,&
         RSimulatedPatterns,ICount,IDisplacements,MPI_DOUBLE_PRECISION,&
         root,MPI_COMM_WORLD,IErr)
    !=====================================
    IF(l_alert(IErr,"SimulateAndFit()","MPI_GATHERV()")) RETURN

    ! put 1D array RSimulatedPatterns into 2D image RImageSimi
    ! remember dimensions of RSimulatedPatterns(INoOfLacbedPatterns,IThicknessCount,IPixelTotal)
    ! and RImageSimi(width, height,INoOfLacbedPatterns,IThicknessCount )
    RImageSimi = ZERO
    ind = 0
    DO jnd = 1,2*IPixelCount
       DO knd = 1,2*IPixelCount
          ind = ind+1
          RImageSimi(jnd,knd,:,:) = RSimulatedPatterns(:,:,ind)
       END DO
    END DO

    ! Gaussian blur to match experiment using global variable RBlurRadius
    IF (RBlurRadius.GT.TINY) THEN
       ALLOCATE(RTempImage(2*IPixelCount,2*IPixelCount),STAT=IErr)
       DO ind=1,INoOfLacbedPatterns
          DO jnd=1,IThicknessCount
             CALL BlurG(RImageSimi(:,:,ind,jnd),IPixelCount,RBlurRadius,IErr)
          END DO
       END DO
    END IF

    ! We have done at least one simulation now
    IInitialSimulationFLAG=0

  END SUBROUTINE FelixFunction

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  !>
  !! Procedure-description: Calculates figure of merit and determines which thickness
  !! matches best. Iteratively 
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!  
  SUBROUTINE CalculateFigureofMeritandDetermineThickness(Iter,IBestThicknessIndex,IErr)

    ! ---> IBestThicknessIndex, RFigureofMerit

    !?? JR called every SimulateAndFit()
    !?? JR called once felixrefine baseline simulation

    !?? NB core 0 only

    USE MyNumbers
    USE message_mod; USE alert_mod 

    USE utilities_mod, ONLY : PhaseCorrelate, ResidualSumofSquares, &
          Normalised2DCrossCorrelation, MaskedCorrelation
    
    ! global outputs
    USE RPARA, ONLY : RFigureofMerit

    ! global inputs
    USE IPARA, ONLY : INoOfLacbedPatterns, ICorrelationFLAG, IPixelCount, IThicknessCount, &
          IImageProcessingFLAG
    USE RPARA, ONLY : RInitialThickness, RDeltaThickness, RImageMask, &
          RImageSimi, &   ! a main input - simulated images
          RImageExpi      ! a main input - experimental images to compare

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(OUT) :: IBestThicknessIndex
    INTEGER(IKIND) :: ind,jnd,knd,IErr,IThickness,hnd,Iter
    INTEGER(IKIND),DIMENSION(INoOfLacbedPatterns) :: IBestImageThicknessIndex
    REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: RSimulatedImage,RExperimentalImage
    REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: RMaskImage
    REAL(RKIND) :: RTotalCorrelation,RBestTotalCorrelation,RImageCorrelation,RBestThickness,&
         RThicknessRange,Rradius
    REAL(RKIND),DIMENSION(INoOfLacbedPatterns) :: RBestCorrelation
    CHARACTER*200 :: SPrintString
    CHARACTER*20 :: Snum       

    IF (ICorrelationFLAG.EQ.3) THEN ! allocate mask
      ALLOCATE(RMaskImage(2*IPixelCount,2*IPixelCount),STAT=IErr)
    END IF
    ! The best correlation for each image will go in here, initialise at the maximum value
    RBestCorrelation = TEN
    RBestTotalCorrelation = TEN ! The best mean of all correlations
    ! The thickness with the lowest figure of merit for each image
    IBestImageThicknessIndex = 1 

    !\/----------------------------------------------------------------------
    DO jnd = 1,IThicknessCount
      RTotalCorrelation = ZERO ! The sum of all individual correlations, initialise at 0
      DO ind = 1,INoOfLacbedPatterns
        RSimulatedImage = RImageSimi(:,:,ind,jnd)
        RExperimentalImage = RImageExpi(:,:,ind)
        IF (ICorrelationFLAG.EQ.3) THEN ! masked correltion, update mask
          RMaskImage=RImageMask(:,:,ind)
        END IF
        
        ! image processing
        SELECT CASE (IImageProcessingFLAG)
        !CASE(0) !no processing
        CASE(1)! square root before perfoming correlation
          RSimulatedImage=SQRT(RSimulatedImage)
          RExperimentalImage=SQRT(RImageExpi(:,:,ind))
        CASE(2)! log before performing correlation
          WHERE (RSimulatedImage.GT.TINY**2)
            RSimulatedImage=LOG(RSimulatedImage)
          ELSEWHERE
            RSimulatedImage = TINY**2
          END WHERE
          WHERE (RExperimentalImage.GT.TINY**2)
            RExperimentalImage = LOG(RImageExpi(:,:,ind))
          ELSEWHERE
            RExperimentalImage =  TINY**2
          END WHERE
        
        END SELECT
        
        ! Correlation type
        SELECT CASE (ICorrelationFLAG)
          CASE(0) ! Phase Correlation
            RImageCorrelation=ONE-& ! NB Perfect Correlation = 0 not 1
                  PhaseCorrelate(RSimulatedImage,RExperimentalImage,&
                  IErr,2*IPixelCount,2*IPixelCount)
          CASE(1) ! Residual Sum of Squares (Non functional)
            RImageCorrelation = ResidualSumofSquares(&
                  RSimulatedImage,RImageExpi(:,:,ind),IErr)
          CASE(2) ! Normalised Cross Correlation
            RImageCorrelation = ONE-& ! NB Perfect Correlation = 0 not 1
                  Normalised2DCrossCorrelation(RSimulatedImage,RExperimentalImage,IErr)
          CASE(3) ! Masked Cross Correlation
            IF (Iter.LE.0) THEN
            ! we are in baseline sim or simplex initialisation: do a normalised2D CC
              RImageCorrelation = ONE-&
                    Normalised2DCrossCorrelation(RSimulatedImage,RExperimentalImage,IErr)   
            ELSE ! we are refining: do a masked CC
              RImageCorrelation = ONE-& ! NB Perfect Correlation = 0 not 1
                    MaskedCorrelation(RSimulatedImage,RExperimentalImage,RMaskImage,IErr)
            END IF
        END SELECT

        CALL message(LL,dbg6,"For Pattern ",ind,", thickness ",jnd)
        CALL message(LL,dbg6,"  the FoM = ",RImageCorrelation)
        
        ! Determines which thickness matches best for each LACBED pattern
        ! which is later used to find the range of viable thicknesses 
        IF(RImageCorrelation.LT.RBestCorrelation(ind)) THEN
          RBestCorrelation(ind) = RImageCorrelation
          IBestImageThicknessIndex(ind) = jnd
        END IF
        RTotalCorrelation = RTotalCorrelation + RImageCorrelation
      END DO
      
      ! The best total correlation
      RTotalCorrelation=RTotalCorrelation/REAL(INoOfLacbedPatterns,RKIND)

      CALL message(LL,dbg6,"Specimen thickness (Angstroms) ",jnd)
      CALL message(LL,dbg6,"Figure of merit ",RTotalCorrelation)

      ! Determines which thickness matches best
      IF(RTotalCorrelation.LT.RBestTotalCorrelation) THEN
        RBestTotalCorrelation = RTotalCorrelation
        IBestThicknessIndex = jnd
      END IF

    END DO
    !/\----------------------------------------------------------------------

    ! The figure of merit, global variable
    RFigureofMerit = RBestTotalCorrelation

    !?? Alternative method below, RWeightingCoefficients not used and has been removed 
    !?? assume that the best thickness is given by the mean of individual thicknesses  
    !IBestThicknessIndex = SUM(IBestImageThicknessIndex)/INoOfLacbedPatterns
    !RBestThickness = RInitialThickness + (IBestThicknessIndex-1)*RDeltaThickness
    !RFigureofMerit = SUM(RBestCorrelation*RWeightingCoefficients)/&
    !      REAL(INoOfLacbedPatterns,RKIND)
    
    ! Output to screen
    RBestThickness = RInitialThickness +(IBestThicknessIndex-1)*RDeltaThickness
    RThicknessRange=( MAXVAL(IBestImageThicknessIndex)-&
          MINVAL(IBestImageThicknessIndex) )*RDeltaThickness

    CALL message(LM,"Figure of merit ",RBestTotalCorrelation)
    CALL message(LM,"Specimen thickness (Angstroms) ",NINT(RBestThickness))
    CALL message(LM,"Thickness range (Angstroms) ",NINT(RThicknessRange))

    RETURN

  END SUBROUTINE CalculateFigureofMeritandDetermineThickness

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Fill the independent parameter array for a new simulation
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!  
  SUBROUTINE UpdateVariables(RIndependentVariable,IErr)

    !?? JR top level to this module, so keep, refine needs access

    !?? JR called every SimulateAndFit() (for non-Ug refinement)

    USE MyNumbers
    USE message_mod; USE alert_mod 

    ! global inputs
    USE IPARA, ONLY : INoOfVariables, IRefineMode, IIterativeVariableUniqueIDs,IAtomMoveList 
    USE RPARA, ONLY : RVector

    ! global outputs
    USE RPARA, ONLY :  RBasisOccupancy, RBasisIsoDW, RAnisotropicDebyeWallerFactorTensor, &
          RLengthX, RLengthY, RLengthZ, RAlpha, RBeta, RGamma, RConvergenceAngle, &
          RAbsorptionPercentage, RAcceleratingVoltage, RRSoSScalingFactor, RBasisAtomPosition

    IMPLICIT NONE

    REAL(RKIND),DIMENSION(INoOfVariables),INTENT(IN) :: RIndependentVariable
    INTEGER(IKIND) :: IVariableType,IVectorID,IAtomID,IErr,ind

    DO ind = 1,INoOfVariables
      IVariableType = IIterativeVariableUniqueIDs(ind,2)
      SELECT CASE (IVariableType)
      CASE(1) ! A: structure factor refinement, do in UpdateStructureFactors
        
      CASE(2)
        ! The index of the atom and vector being used
        IVectorID = IIterativeVariableUniqueIDs(ind,3)
        ! The atom being moved
        IAtomID = IAtomMoveList(IVectorID)
        ! Change in position r' = r - v*(r.v) +v*RIndependentVariable(ind)
        RBasisAtomPosition(IAtomID,:) = RBasisAtomPosition(IAtomID,:) - &
            RVector(IVectorID,:)*DOT_PRODUCT(RBasisAtomPosition(IAtomID,:),RVector(IVectorID,:)) + &
            RVector(IVectorID,:)*RIndependentVariable(ind)
      CASE(3)
        RBasisOccupancy(IIterativeVariableUniqueIDs(ind,3))=RIndependentVariable(ind) 
      CASE(4)
        RBasisIsoDW(IIterativeVariableUniqueIDs(ind,3))=RIndependentVariable(ind)
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
      END SELECT
    END DO

  END SUBROUTINE UpdateVariables

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  !>
  !! Procedure-description: Print variables
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !! 
  SUBROUTINE PrintVariables(IErr)

    !?? JR utility, however is top level to this module so maybe keep 
    !?? JR uninteresting, low priority, check if relevant performance overhead

    !?? JR called every SimulateAndFit()

    USE MyNumbers
    USE message_mod; USE alert_mod 
    
    USE IConst; USE RConst; USE SConst
    USE IPara; USE RPara; USE CPara; USE SPara;
    USE BlochPara 

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr,ind,IVariableType,jnd,knd
    REAL(RKIND),DIMENSION(3) :: RCrystalVector
    CHARACTER*200 :: SPrintString

    RCrystalVector = [RLengthX,RLengthY,RLengthZ]

    DO ind = 1,IRefinementVariableTypes
      IF (IRefineMode(ind).EQ.1) THEN
        SELECT CASE(ind)
        CASE(1)
          IF (IAbsorbFLAG.EQ.1) THEN ! proportional absorption
            CALL message(LS,"Current Absorption ",RAbsorptionPercentage)
          END IF
          CALL message(LS,"Current Structure Factors nm^-2: amplitude, phase (deg)")
          !?? RB should also put in hkl here
          !DO jnd = 1+IUgOffset,INoofUgs+IUgOffset
          !  WRITE(SPrintString,FMT='(2(1X,F7.3),2X,A1,1X,F6.3,1X,F6.2)') 100*CUniqueUg(jnd),":",&
          !        ABS(CUniqueUg(jnd)),180*ATAN2(AIMAG(CUniqueUg(jnd)),REAL(CUniqueUg(jnd)))/PI
          !END DO
          DO jnd = 1+IUgOffset,INoofUgs+IUgOffset
            CALL message(LM,"",100*CUniqueUg(jnd))
            CALL message(LM,"",(/ ABS(CUniqueUg(jnd)),180*&
                  ATAN2(AIMAG(CUniqueUg(jnd)),REAL(CUniqueUg(jnd)))/PI /))
          END DO

        CASE(2)
          CALL message(LM,"Current Atomic Coordinates")
          DO jnd = 1,SIZE(RBasisAtomPosition,DIM=1)
            CALL message(LM,SBasisAtomLabel(jnd),RBasisAtomPosition(jnd,:))
          END DO

        CASE(3)
          CALL message(LM, "Current Atomic Occupancy")
          DO jnd = 1,SIZE(RBasisOccupancy,DIM=1)
            CALL message(LM, SBasisAtomLabel(jnd),RBasisOccupancy(jnd))
          END DO

        CASE(4)
          CALL message(LM, "Current Isotropic Debye Waller Factors")
          DO jnd = 1,SIZE(RBasisIsoDW,DIM=1)
            CALL message(LM, SBasisAtomLabel(jnd),RBasisIsoDW(jnd))
          END DO

        CASE(5)
          CALL message(LM,"Current Anisotropic Debye Waller Factors")
          DO jnd = 1,SIZE(RAnisotropicDebyeWallerFactorTensor,DIM=1)
            CALL message(LM, "Tensor index = ", jnd) 
            CALL message(LM, SBasisAtomLabel(jnd),RAnisotropicDebyeWallerFactorTensor(jnd,1:3,:) )
          END DO

        CASE(6)
          CALL message(LM, "Current Unit Cell Parameters", (/ RLengthX,RLengthY,RLengthZ /) )

        CASE(7)
          CALL message(LM, "Current Unit Cell Angles", (/ RAlpha,RBeta,RGamma /) )

        CASE(8)
          CALL message(LM, "Current Convergence Angle: ", RConvergenceAngle )

        CASE(9)
          CALL message(LM, "Current Absorption Percentage", RAbsorptionPercentage )

        CASE(10)
          CALL message(LM, "Current Accelerating Voltage", RAcceleratingVoltage )

        CASE(11)
          CALL message(LM, "Current Residual Sum of Squares Scaling Factor", RRSoSScalingFactor )

        END SELECT
      END IF
    END DO

  END SUBROUTINE PrintVariables

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  !>
  !! Procedure-description: Performs a 2D Gaussian blur on the input image using 
  !! global variable RBlurRadius and renormalises the output image to have the
  !! same min and max as the input image
  !!
  !! Closed procedure, no access to global variables
  !!
  !! Major-Authors: Richard Beanland (2016)
  !! 
  PURE SUBROUTINE BlurG(RImageToBlur,IPixelsCount,RBlurringRadius,IErr)

    !?? JR niche utility, seemingly unimportant, could be moved

    !?? JR called every SimulateAndFit() indriectly via FelixFunction

    USE MyNumbers
    USE MPI
    USE message_mod; USE alert_mod !?? can't use this or l_alert if pure

    IMPLICIT NONE

    REAL(RKIND),DIMENSION(2*IPixelsCount,2*IPixelsCount),INTENT(INOUT) :: RImageToBlur
    INTEGER(IKIND),INTENT(IN) :: IPixelsCount
    REAL(RKIND),INTENT(IN) :: RBlurringRadius
    INTEGER(IKIND),INTENT(OUT) :: IErr

    REAL(RKIND),DIMENSION(2*IPixelsCount,2*IPixelsCount) :: RTempImage,RShiftImage
    INTEGER(IKIND) :: ind,jnd,IKernelRadius,IKernelSize
    REAL(RKIND),DIMENSION(:), ALLOCATABLE :: RGauss1D
    REAL(RKIND) :: Rind,Rsum,Rmin,Rmax

    

    ! get min and max of input image
    Rmin=MINVAL(RImageToBlur)
    Rmax=MAXVAL(RImageToBlur)

    ! set up a 1D kernel of appropriate size  
    IKernelRadius=NINT(3*RBlurringRadius)
    ALLOCATE(RGauss1D(2*IKernelRadius+1),STAT=IErr)!ffs
    Rsum=0
    DO ind=-IKernelRadius,IKernelRadius
      Rind=REAL(ind)
      RGauss1D(ind+IKernelRadius+1)=EXP(-(Rind**2)/((2*RBlurringRadius)**2))
      Rsum=Rsum+RGauss1D(ind+IKernelRadius+1)
      IF(ind==0) IErr=78 
    END DO
    RGauss1D=RGauss1D/Rsum!normalise
    RTempImage=RImageToBlur*0_RKIND !reset the temp image 

    ! apply the kernel in direction 1
    DO ind = -IKernelRadius,IKernelRadius
       IF (ind.LT.0) THEN
          RShiftImage(1:2*IPixelsCount+ind,:)=RImageToBlur(1-ind:2*IPixelsCount,:)
          DO jnd = 1,1-ind!edge fill on right
             RShiftImage(2*IPixelsCount-jnd+1,:)=RImageToBlur(2*IPixelsCount,:)
          END DO
       ELSE
          RShiftImage(1+ind:2*IPixelsCount,:)=RImageToBlur(1:2*IPixelsCount-ind,:)
          DO jnd = 1,1+ind!edge fill on left
             RShiftImage(jnd,:)=RImageToBlur(1,:)
          END DO
       END IF
       RTempImage=RTempImage+RShiftImage*RGauss1D(ind+IKernelRadius+1)
    END DO

    ! make the 1D blurred image the input for the next direction
    RImageToBlur=RTempImage
    RTempImage=RImageToBlur*0_RKIND ! reset the temp image

    ! apply the kernel in direction 2  
    DO ind = -IKernelRadius,IKernelRadius
       IF (ind.LT.0) THEN
          RShiftImage(:,1:2*IPixelsCount+ind)=RImageToBlur(:,1-ind:2*IPixelsCount)
          DO jnd = 1,1-ind!edge fill on bottom
             RShiftImage(:,2*IPixelsCount-jnd+1)=RImageToBlur(:,2*IPixelsCount)
          END DO
       ELSE
          RShiftImage(:,1+ind:2*IPixelsCount)=RImageToBlur(:,1:2*IPixelsCount-ind)
          DO jnd = 1,1+ind!edge fill on top
             RShiftImage(:,jnd)=RImageToBlur(:,1)
          END DO
       END IF
       RTempImage=RTempImage+RShiftImage*RGauss1D(ind+IKernelRadius+1)
    END DO
    DEALLOCATE(RGauss1D,STAT=IErr)

    ! set intensity range of outpt image to match that of the input image
    RTempImage=RTempImage-MINVAL(RTempImage)
    RTempImage=RTempImage*(Rmax-Rmin)/MAXVAL(RTempImage)+Rmin
    ! return the blurred image
    RImageToBlur=RTempImage;

  END SUBROUTINE BlurG        

END MODULE felixfunction_mod    
