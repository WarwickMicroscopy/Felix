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

!!$REAL(RKIND) FUNCTION FelixFunction(IIterationFLAG,IErr)
SUBROUTINE FelixFunction(LInitialSimulationFLAG,IErr)

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
  
  INTEGER(IKIND) :: &
       IErr,ind,jnd,knd,pnd,&
       IThicknessIndex,ILocalPixelCountMin, ILocalPixelCountMax,&
       IIterationFLAG
  INTEGER(IKIND) :: &
       IAbsorbTag = 0
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       IDisplacements,ICount
  LOGICAL,INTENT(IN) :: &
       LInitialSimulationFLAG !If function is being called during initialisation
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RIndividualReflectionsRoot,&
       RFinalMontageImageRoot
  COMPLEX(CKIND),DIMENSION(:,:,:), ALLOCATABLE :: &
       CAmplitudeandPhaseRoot 

  IF(IWriteFLAG.GE.10.AND.my_rank.EQ.0) THEN
     PRINT*,"Felix function"
  END IF

  CALL CountTotalAtoms(IErr)
  
  IDiffractionFLAG = 0

  !-------------------------------------------------------------------- 
  !Setup Experimental Variables
  !--------------------------------------------------------------------

  CALL ExperimentalSetup (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error in ExperimentalSetup()"
     RETURN
  END IF
  
  
  !--------------------------------------------------------------------
  ! Setup Image
  !--------------------------------------------------------------------


  CALL ImageSetup( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error in ImageSetup()"
     RETURN
  END IF

 
  !--------------------------------------------------------------------
  ! MAIN section
  !--------------------------------------------------------------------
 
!!$  Structure Factors must be calculated without absorption for refinement to work

  !RB IF(IAbsorbFLAG.NE.0) THEN 
  !RB   IAbsorbFLAG = 0 ! Non-absorpative structure factor calculation
  !RB   IAbsorbTAG = 1 ! Remember that IAbsorbFLAG was 1
  !RB END IF
  
  CALL StructureFactorSetup(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ",IErr,"in StructureFactorSetup()"
     RETURN
  END IF

  !RB IF(IAbsorbTAG.NE.0) IAbsorbFLAG = 1 !Reset IAbsorbFLAG to 1

  IF((IRefineModeSelectionArray(1).EQ.1).AND.(LInitialSimulationFLAG.NEQV..TRUE.)) THEN
     
     CALL ApplyNewStructureFactors(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in ApplyNewStructureFactors()"
        RETURN
     END IF
     
  END IF

  IF(IAbsorbFLAG.NE.0) THEN
     
     CALL StructureFactorsWithAbsorptionDetermination(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ",IErr,&
             "in StructureFactorsWithAbsorptionDetermination()"
        RETURN
     END IF

  END IF     
  !--------------------------------------------------------------------
  ! reserve memory for effective eigenvalue problem
  !--------------------------------------------------------------------

  !Kprime Vectors and Deviation Parameter
  
  ALLOCATE(RDevPara(nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RDevPara"
     RETURN
  END IF

  ALLOCATE(IStrongBeamList(nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables IStrongBeamList"
     RETURN
  END IF

  ALLOCATE(IWeakBeamList(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables IWeakBeamList"
     RETURN
  END IF

  !--------------------------------------------------------------------
  ! MAIN LOOP: solve for each (ind,jnd) pixel
  !--------------------------------------------------------------------

  ILocalPixelCountMin= (IPixelTotal*(my_rank)/p)+1
  ILocalPixelCountMax= (IPixelTotal*(my_rank+1)/p) 

  IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Felixfunction(", my_rank, "): starting the eigenvalue problem"
     PRINT*,"Felixfunction(", my_rank, "): for lines ", ILocalPixelCountMin, &
          " to ", ILocalPixelCountMax
  END IF
  
  IThicknessCount= (RFinalThickness- RInitialThickness)/RDeltaThickness + 1

  IF(IImageFLAG.LE.2) THEN
     ALLOCATE(RIndividualReflections(IReflectOut,IThicknessCount,&
          (ILocalPixelCountMax-ILocalPixelCountMin)+1),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Individual Images"
        RETURN
     END IF
     
     RIndividualReflections = ZERO
  ELSE
     ALLOCATE(CAmplitudeandPhase(IReflectOut,IThicknessCount,&
          (ILocalPixelCountMax-ILocalPixelCountMin)+1),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Amplitude and Phase"
        RETURN
     END IF
     CAmplitudeandPhase = CZERO
  END IF

  ALLOCATE(CFullWaveFunctions(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CFullWaveFunctions"
     RETURN
  END IF
  
  ALLOCATE(RFullWaveIntensity(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RFullWaveIntensity"
     RETURN
  END IF  

  IMAXCBuffer = 200000
  IPixelComputed= 0
  
  IF(IWriteFLAG.GE.10.AND.my_rank.EQ.0) THEN
     PRINT*,"Felixfunction(",my_rank,") Entering BlochLoop()"
  END IF

  DO knd = ILocalPixelCountMin,ILocalPixelCountMax,1
     jnd = IPixelLocations(knd,1)
     ind = IPixelLocations(knd,2)
     CALL BlochCoefficientCalculation(ind,jnd,knd,ILocalPixelCountMin,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in BlochCofficientCalculation"
        RETURN
     END IF
  END DO
  
  IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Felixfunction : ",my_rank," is exiting calculation loop"
  END IF

  !--------------------------------------------------------------------
  ! close outfiles
  !--------------------------------------------------------------------
  
  ALLOCATE(RIndividualReflectionsRoot(IReflectOut,IThicknessCount,IPixelTotal),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Root Reflections"
     RETURN
  END IF
  
  IF(IImageFLAG.GE.3) THEN
     ALLOCATE(CAmplitudeandPhaseRoot(IReflectOut,IThicknessCount,IPixelTotal),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Root Amplitude and Phase"
        RETURN
     END IF
     CAmplitudeandPhaseRoot = CZERO
  END IF

  RIndividualReflectionsRoot = ZERO

  ALLOCATE(IDisplacements(p),ICount(p),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " In ALLOCATE"
     RETURN
  END IF

  DO pnd = 1,p
     IDisplacements(pnd) = (IPixelTotal*(pnd-1)/p)
     ICount(pnd) = (((IPixelTotal*(pnd)/p) - (IPixelTotal*(pnd-1)/p)))*IReflectOut*IThicknessCount
          
  END DO
  
  DO ind = 1,p
        IDisplacements(ind) = (IDisplacements(ind))*IReflectOut*IThicknessCount
  END DO
  
  IF(IImageFLAG.LE.2) THEN
     CALL MPI_GATHERV(RIndividualReflections,SIZE(RIndividualReflections),&
          MPI_DOUBLE_PRECISION,RIndividualReflectionsRoot,&
          ICount,IDisplacements,MPI_DOUBLE_PRECISION,0,&
          MPI_COMM_WORLD,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " In MPI_GATHERV"
        RETURN
     END IF     
  ELSE     
     CALL MPI_GATHERV(CAmplitudeandPhase,SIZE(CAmplitudeandPhase),&
          MPI_DOUBLE_COMPLEX,CAmplitudeandPhaseRoot,&
          ICount,IDisplacements,MPI_DOUBLE_COMPLEX,0, &
          MPI_COMM_WORLD,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " In MPI_GATHERV"
        RETURN
     END IF   
  END IF
  
  IF(IImageFLAG.GE.3) THEN
     DEALLOCATE(CAmplitudeandPhase,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " Deallocating CAmplitudePhase"
        RETURN
     END IF   
  END IF
   
  IF(IImageFLAG.LE.2) THEN
     DEALLOCATE(RIndividualReflections,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " Deallocating RIndividualReflections"
        RETURN
     END IF   
  END IF
  
  IF(my_rank.EQ.0.AND.IImageFLAG.GE.3) THEN
     RIndividualReflectionsRoot = &
          CAmplitudeandPhaseRoot * CONJG(CAmplitudeandPhaseRoot)
  END IF
  
  IF(my_rank.EQ.0) THEN
     ALLOCATE(RIndividualReflections(IReflectOut,IThicknessCount,IPixelTotal),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Root Reflections"
        RETURN
     END IF
     
     RIndividualReflections = RIndividualReflectionsRoot
  END IF
  
  !--------------------------------------------------------------------
  ! free memory
  !--------------------------------------------------------------------
  
  !Dellocate Global Variables
  
  DEALLOCATE(RMask,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error  ", IErr, &
          " in Deallocation of RMask etc"
     RETURN
  END IF

  DEALLOCATE(RDevPara,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating RDevPara"
     RETURN
  END IF

  DEALLOCATE(IPixelLocations,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating IPixelLocations"
     RETURN
  END IF

  DEALLOCATE(IStrongBeamList,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating IStrongBeamList"
     RETURN
  END IF

  DEALLOCATE(IWeakBeamList,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating IWeakBeamList"
     RETURN
  END IF

  DEALLOCATE(IDisplacements,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating IDisplacements"
     RETURN
  END IF

  DEALLOCATE(ICount,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating ICount"
     RETURN
  END IF
  
  DEALLOCATE(CFullWaveFunctions, & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CFullWaveFunctions"
     RETURN
  END IF
  
  DEALLOCATE(RFullWaveIntensity,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RFullWaveIntensity"
     RETURN
  END IF  

  IF(IImageFLAG.GE.3) THEN
     DEALLOCATE(CAmplitudeandPhaseRoot,STAT=IErr)     
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error in Deallocation of CAmplitudeandPhase"
        RETURN  
     END IF
  END IF

  DEALLOCATE(RIndividualReflectionsRoot,STAT=IErr) 
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error in Deallocation of RIndividualReflectionsRoot "
     RETURN  
  END IF
  
  IF((my_rank.NE.0).AND.(LInitialSimulationFLAG.NEQV..TRUE.)) THEN     
     DEALLOCATE(Rhkl,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in Deallocation Rhkl"
        RETURN
     END IF
  END IF
  
!!$  IF (my_rank.NE.0) THEN
  
  DEALLOCATE(Rhklpositions,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating Rhklpositions"
     RETURN
  END IF
!!$   END IF
  
  DEALLOCATE(MNP,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation MNP"
     RETURN
  END IF
      
  DEALLOCATE(SMNP,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation SMNP"
     RETURN
  END IF
       
  DEALLOCATE(RDWF,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RDWF"
     RETURN
  END IF
       
  DEALLOCATE(ROcc,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation ROcc"
     RETURN
  END IF
       
  DEALLOCATE(IAtoms,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation IAtoms"
     RETURN
  END IF
       
  DEALLOCATE(IAnisoDWFT,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation IAnisoDWFT"
     RETURN
  END IF
       
  DEALLOCATE(RFullAtomicFracCoordVec,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RFullAtomicFracCoordVec"
     RETURN
  END IF
       
  DEALLOCATE(SFullAtomicNameVec,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation SFullAtomicNameVec"
     RETURN
  END IF
       
  DEALLOCATE(RFullPartialOccupancy,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RFullPartialOccupancy"
     RETURN
  END IF
       
  DEALLOCATE(RFullIsotropicDebyeWallerFactor,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RFullIsotropicDebyeWallerFactor"
     RETURN
  END IF
       
  DEALLOCATE(IFullAtomNumber,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation IFullAtomNumber"
     RETURN
  END IF
       
  DEALLOCATE(IFullAnisotropicDWFTensor,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation IFullAnisotropicDWFTensor"
     RETURN
  END IF
       
  DEALLOCATE(RgVecMatT,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMatT"
     RETURN
  END IF
       
  DEALLOCATE(RgVecMag,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMag"
     RETURN
  END IF
       
  DEALLOCATE(RgVecVec,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMag"
     RETURN
  END IF
!XX
!  DEALLOCATE(RgSumMat,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, " in DEALLOCATE of RrVecMat"
     RETURN
  END IF

END SUBROUTINE FelixFunction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CalculateFigureofMeritandDetermineThickness(IThicknessCountFinal,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: ind,jnd,knd,IErr,ICountedPixels,IThickness,hnd
  INTEGER(IKIND),DIMENSION(IReflectOut) :: IThicknessByReflection
  INTEGER(IKIND),INTENT(OUT) :: IThicknessCountFinal
  REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: RSimulatedImageForPhaseCorrelation,RExperimentalImage
  REAL(RKIND) :: RCrossCorrelationOld,RIndependentCrossCorrelation,RThickness,&
	   PhaseCorrelate,Normalised2DCrossCorrelation,ResidualSumofSquares
  REAL(RKIND),DIMENSION(IReflectOut) :: RReflectionCrossCorrelations,RReflectionThickness
  CHARACTER*200 :: SPrintString
       
  
  IF(IWriteFLAG.GE.10.AND.my_rank.EQ.0) THEN
     PRINT*,"CalculateFigureofMeritandDetermineThickness(",my_rank,")"
  END IF

  RReflectionCrossCorrelations = ZERO

  DO hnd = 1,IReflectOut
     RCrossCorrelationOld = 1.0E15 !A large Number
     RThickness = ZERO
     DO ind = 1,IThicknessCount
        
        ICountedPixels = 0

        RSimulatedImageForPhaseCorrelation = ZERO

        RIndependentCrossCorrelation = ZERO

        DO jnd = 1,2*IPixelCount
           DO knd = 1,2*IPixelCount
              IF(ABS(RMask(jnd,knd)).GT.TINY) THEN
                 ICountedPixels = ICountedPixels+1
                                  
                 RSimulatedImageForPhaseCorrelation(jnd,knd) = &
                      RIndividualReflections(hnd,ind,ICountedPixels)
              END IF
           END DO
        END DO
                
        SELECT CASE (IImageProcessingFLAG)
        CASE(0)
           RExperimentalImage = RImageExpi(:,:,hnd)
        CASE(1)
           RSimulatedImageForPhaseCorrelation = &
                SQRT(RSimulatedImageForPhaseCorrelation)
           RExperimentalImage = &
                SQRT(RImageExpi(:,:,hnd))
        CASE(2)
           WHERE (RSimulatedImageForPhaseCorrelation.GT.TINY**2)
              RSimulatedImageForPhaseCorrelation = &
                   LOG(RSimulatedImageForPhaseCorrelation)
           ELSEWHERE
              RSimulatedImageForPhaseCorrelation = &
                   TINY**2
           END WHERE
              
           WHERE (RExperimentalImage.GT.TINY**2)
              RExperimentalImage = &
                   LOG(RImageExpi(:,:,hnd))
           ELSEWHERE
              RExperimentalImage = &
                   TINY**2
           END WHERE
              
        END SELECT

        SELECT CASE (ICorrelationFLAG)
           
        CASE(0) ! Phase Correlation
           
           RIndependentCrossCorrelation = &
                ONE-& ! So Perfect Correlation = 0 not 1
                PhaseCorrelate(&
                RSimulatedImageForPhaseCorrelation,RExperimentalImage,&
                IErr,2*IPixelCount,2*IPixelCount)
           
        CASE(1) ! Residual Sum of Squares (Non functional)
           RIndependentCrossCorrelation = &
                ResidualSumofSquares(&
                RSimulatedImageForPhaseCorrelation,RImageExpi(:,:,hnd),IErr)
           
        CASE(2) ! Normalised Cross Correlation

           RIndependentCrossCorrelation = &
                ONE-& ! So Perfect Correlation = 0 not 1
                Normalised2DCrossCorrelation(&
                RSimulatedImageForPhaseCorrelation,RExperimentalImage,&
                (/2*IPixelCount, 2*IPixelCount/),IPixelTotal,IErr)
           
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

  RCrossCorrelation = &
       SUM(RReflectionCrossCorrelations*RWeightingCoefficients)/&
       REAL(IReflectOut,RKIND)
!RB assume that the thickness is given by the mean of individual thicknesses  
  IThicknessCountFinal = SUM(IThicknessByReflection)/IReflectOut

  RThickness = RInitialThickness + (IThicknessCountFinal-1)*RDeltaThickness 
  

  IF(my_rank.eq.0) THEN
!RB     PRINT*,"Thicknesses",RReflectionThickness
!RB     PRINT*,"Correlation",RCrossCorrelation
!RB     PRINT*,"Thickness Final",IThicknessCountFinal
     WRITE(SPrintString,FMT='(A18,I4,A10)') "Specimen thickness ",NINT(RThickness)," Angstroms"
     PRINT*,TRIM(ADJUSTL(SPrintString))
     !XXPRINT*,"---------------------------------------------------------"
!     PRINT*,"Specimen thickness",NINT(RThickness),"Angstroms"
  END IF

  IF (RCrossCorrelation.NE.RCrossCorrelation) THEN!RB what?
     IErr = 1
  END IF

END SUBROUTINE CalculateFigureofMeritandDetermineThickness

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REAL(RKIND) FUNCTION SimplexFunction(RIndependentVariableValues,IIterationCount,IExitFLAG,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,IExitFLAG,IThickness
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(INOUT) :: RIndependentVariableValues
  INTEGER(IKIND),INTENT(IN) :: IIterationCount
  LOGICAL :: LInitialSimulationFLAG = .FALSE.
  
  IF(IWriteFLAG.GE.10.AND.my_rank.EQ.0) THEN
     PRINT*,"SimplexFunction(",my_rank,")"
  END IF
  

  IF(IRefineModeSelectionArray(1).EQ.1) THEN     
     CALL UpdateStructureFactors(RIndependentVariableValues,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in UpdateStructureFactors"
        RETURN
     END IF     
  ELSE
     CALL UpdateVariables(RIndependentVariableValues,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in UpdateVariables"
        RETURN
     END IF
!RB     WHERE(RAtomSiteFracCoordVec.LT.0) RAtomSiteFracCoordVec=RAtomSiteFracCoordVec+ONE
!RB     WHERE(RAtomSiteFracCoordVec.GT.1) RAtomSiteFracCoordVec=RAtomSiteFracCoordVec-ONE
  END IF
  WHERE(RAtomSiteFracCoordVec.LT.0) RAtomSiteFracCoordVec=RAtomSiteFracCoordVec+ONE
  WHERE(RAtomSiteFracCoordVec.GT.1) RAtomSiteFracCoordVec=RAtomSiteFracCoordVec-ONE

  IF (my_rank.EQ.0) THEN
     CALL PrintVariables(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in PrintVariables"
        RETURN
     END IF
  END IF

  CALL FelixFunction(LInitialSimulationFLAG,IErr) ! Simulate !!  
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
          " in FelixFunction"
     RETURN
  END IF

  IF(my_rank.EQ.0) THEN   
     CALL CreateImagesAndWriteOutput(IIterationCount,IExitFLAG,IErr) 
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in CreateImagesAndWriteOutput"
        RETURN
     ENDIF
     
     SimplexFunction = RCrossCorrelation     
  END IF

!RB   NB Also deallocated in felixrefine!!!
  DEALLOCATE(RgSumMat,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation of RgSumMat"
     RETURN
  ENDIF
!RB   PRINT*,"Deallocating CUgMatNoAbs,CUgMatPrime,CUgMat in SimplexFunction" 
  DEALLOCATE(CUgMatNoAbs,STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation (", my_rank, ") error in Deallocation()"
     RETURN
  ENDIF
 
 DEALLOCATE(CUgMatPrime,STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation (", my_rank, ") error in Deallocation()"
     RETURN
  ENDIF
 
  DEALLOCATE(CUgMat,STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation (", my_rank, ") error in Deallocation()"
     RETURN
  ENDIF
  
END FUNCTION SimplexFunction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CreateImagesAndWriteOutput(IIterationCount,IExitFLAG,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,IThicknessIndex,IIterationCount,IExitFLAG
  
  ALLOCATE(RMask(2*IPixelCount,2*IPixelCount),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateImagesAndWriteOutput(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variable RMask"
     RETURN
  ENDIF
  
  CALL ImageMaskInitialisation(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateImagesAndWriteOutput(", my_rank, ") error ", IErr, &
          " in ImageMaskInitialisation"
     RETURN
  ENDIF
  
  CALL CalculateFigureofMeritandDetermineThickness(IThicknessIndex,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateImagesAndWriteOutput(", my_rank, ") error ", IErr, &
          "Calling function CalculateFigureofMeritandDetermineThickness"
     RETURN
  ENDIF
  
  DEALLOCATE(RMask,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateImagesAndWriteOutput(", my_rank, ") error ", IErr, &
          " in DEALLOCATE() of DYNAMIC variable RMask"
     RETURN
  ENDIF
  
!!$     OUTPUT -------------------------------------
  
  CALL WriteIterationOutput(IIterationCount,IThicknessIndex,IExitFLAG,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateImagesAndWriteOutput(", my_rank, ") error ", IErr, &
          " in WriteIterationOutput"
     RETURN
  ENDIF

!!$     FINISH OUTPUT  --------------------------------
  
  DEALLOCATE(RIndividualReflections,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateImagesAndWriteOutput(", my_rank, ") error ", IErr, &
          " Deallocating RIndividualReflections"
     RETURN
  ENDIF
  
  DEALLOCATE(IPixelLocations,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateImagesAndWriteOutput(", my_rank, ") error ", IErr, &
          " Deallocating IPixelLocations"
     RETURN
  ENDIF
       
  DEALLOCATE(Rhkl,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateImagesAndWriteOutput(", my_rank, ") error ", IErr, &
          " in Deallocation Rhkl"
     RETURN
  ENDIF
END SUBROUTINE CreateImagesAndWriteOutput

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE UpdateVariables(RIndependentVariableValues,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IVariableType,IErr,ind
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(IN) :: RIndependentVariableValues

  !!$  Fill the Independent Value array with values
  
  
  IF(IRefineModeSelectionArray(2).EQ.1) THEN     
     RAtomSiteFracCoordVec = RInitialAtomSiteFracCoordVec
  END IF

  DO ind = 1,IIndependentVariables
     IVariableType = IIterativeVariableUniqueIDs(ind,2)
     SELECT CASE (IVariableType)
     CASE(1)
	    !RB structure factor refinement, do in UpdateStructureFactors
     CASE(2)
        CALL ConvertVectorMovementsIntoAtomicCoordinates(ind,RIndependentVariableValues,IErr)
     CASE(3)
        RAtomicSitePartialOccupancy(IIterativeVariableUniqueIDs(ind,3)) = &
             RIndependentVariableValues(ind)
     CASE(4)
        RIsotropicDebyeWallerFactors(IIterativeVariableUniqueIDs(ind,3)) = &
             RIndependentVariableValues(ind)
     CASE(5)
        RAnisotropicDebyeWallerFactorTensor(&
             IIterativeVariableUniqueIDs(ind,3),&
             IIterativeVariableUniqueIDs(ind,4),&
             IIterativeVariableUniqueIDs(ind,5)) = & 
             RIndependentVariableValues(ind)
     CASE(6)
        SELECT CASE(IIterativeVariableUniqueIDs(ind,3))
        CASE(1)
           RLengthX = RIndependentVariableValues(ind)
        CASE(2)
           RLengthY = RIndependentVariableValues(ind)
        CASE(3)
           RLengthZ = RIndependentVariableValues(ind)
        END SELECT
     CASE(7)
        SELECT CASE(IIterativeVariableUniqueIDs(ind,3))
        CASE(1)
           RAlpha = RIndependentVariableValues(ind)
        CASE(2)
           RBeta = RIndependentVariableValues(ind)
        CASE(3)
           RGamma = RIndependentVariableValues(ind)
        END SELECT
     CASE(8)
        RConvergenceAngle = RIndependentVariableValues(ind)
     CASE(9)
        RAbsorptionPercentage = RIndependentVariableValues(ind)
     CASE(10)
        RAcceleratingVoltage = RIndependentVariableValues(ind)
     CASE(11)
        RRSoSScalingFactor = RIndependentVariableValues(ind)
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
  !REAL(RKIND) :: &!RB
  !     RUgAmplitude,RUgPhase!RB
  CHARACTER*200 :: SPrintString

  RCrystalVector = [RLengthX,RLengthY,RLengthZ]

  DO ind = 1,IRefinementVariableTypes
     IF (IRefineModeSelectionArray(ind).EQ.1) THEN
        SELECT CASE(ind)
        CASE(1)
           WRITE(SPrintString,FMT='(A18,1X,F9.4)') "Current Absorption",RAbsorptionPercentage
           PRINT*,TRIM(ADJUSTL(SPrintString))
!           PRINT*,"Current Absorption",RAbsorptionPercentage
           PRINT*,"Current Structure Factors"!RB should also put in hkl here
           DO jnd = 2,INoofUgs+1!yy since no.1 is 000
   !           RUgAmplitude=( REAL(CSymmetryStrengthKey(jnd))**2 + AIMAG(CSymmetryStrengthKey(jnd))**2 )**0.5!RB
   !           RUgPhase=ATAN2(AIMAG(CSymmetryStrengthKey(jnd)),REAL(CSymmetryStrengthKey(jnd)))*180/PI!RB
              WRITE(SPrintString,FMT='(2(1X,F9.4))') REAL(CSymmetryStrengthKey(jnd)),AIMAG(CSymmetryStrengthKey(jnd))
              PRINT*,TRIM(ADJUSTL(SPrintString))
!XX              PRINT*,CSymmetryStrengthKey(jnd)!,": Amplitude ",RUgAmplitude,", phase ",RUgPhase
           END DO           
        CASE(2)
           PRINT*,"Current Atomic Coordinates"
           DO jnd = 1,SIZE(RAtomSiteFracCoordVec,DIM=1)
              WRITE(SPrintString,FMT='(A2,3(1X,F9.4))') SAtomName(jnd),RAtomSiteFracCoordVec(jnd,:)
              PRINT*,TRIM(ADJUSTL(SPrintString))              
           END DO
        CASE(3)
           PRINT*,"Current Atomic Occupancy"
           DO jnd = 1,SIZE(RAtomicSitePartialOccupancy,DIM=1)
              WRITE(SPrintString,FMT='(A2,1X,F9.6)') SAtomName(jnd),RAtomicSitePartialOccupancy(jnd)
              PRINT*,TRIM(ADJUSTL(SPrintString))
           END DO
        CASE(4)
           PRINT*,"Current Isotropic Debye Waller Factors"
           DO jnd = 1,SIZE(RIsotropicDebyeWallerFactors,DIM=1)
              WRITE(SPrintString,FMT='(A2,1X,F9.6)') SAtomName(jnd),RIsotropicDebyeWallerFactors(jnd)
              PRINT*,TRIM(ADJUSTL(SPrintString))
           END DO
        CASE(5)
           PRINT*,"Current Anisotropic Debye Waller Factors"
           DO jnd = 1,SIZE(RAnisotropicDebyeWallerFactorTensor,DIM=1)
              DO knd = 1,3
                 WRITE(SPrintString,FMT='(A2,3(1X,F9.4))') SAtomName(jnd),RAnisotropicDebyeWallerFactorTensor(jnd,knd,:)
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

SUBROUTINE UpdateStructureFactors(RIndependentVariableValues,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,ind
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(IN) :: RIndependentVariableValues

  IF(IRefineModeSelectionArray(1).EQ.1) THEN
     DO ind = 1,INoofUgs
        CSymmetryStrengthKey(ind+1) = &!yy ind+1 instead of ind
             CMPLX(RIndependentVariableValues((ind-1)*2+1),RIndependentVariableValues((ind-1)*2+2),CKIND)
     END DO
	 RAbsorptionPercentage = RIndependentVariableValues(2*INoofUgs+1)!RB
  END IF
  
END SUBROUTINE UpdateStructureFactors

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ConvertVectorMovementsIntoAtomicCoordinates(IVariableID,RIndependentVariableValues,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,jnd,IVariableID,IVectorID,IAtomID
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(IN) :: &
       RIndependentVariableValues

!!$  Use IVariableID to determine which vector is being applied (IVectorID)

  IVectorID = IIterativeVariableUniqueIDs(IVariableID,3)

!!$  Use IVectorID to determine which atomic coordinate the vector is to be applied to (IAtomID)

  IAtomID = IAllowedVectorIDs(IVectorID)

!!$  Use IAtomID to applied the IVectodID Vector to the IAtomID atomic coordinate
    
  RAtomSiteFracCoordVec(IAtomID,:) = RAtomSiteFracCoordVec(IAtomID,:) + &
       RIndependentVariableValues(IVariableID)*RAllowedVectors(IVectorID,:)
  
END SUBROUTINE ConvertVectorMovementsIntoAtomicCoordinates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE InitialiseWeightingCoefficients(IErr)
  
  USE MyNumbers
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RWeightingCoefficientsDummy

  ALLOCATE(RWeightingCoefficients(IReflectOut),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"InitialiseWeightingCoefficients(", my_rank, ") error ", IErr, &
          " in allocation RWeightingCoefficients"
     RETURN
  ENDIF

  ALLOCATE(RWeightingCoefficientsDummy(IReflectOut),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"InitialiseWeightingCoefficients(", my_rank, ") error ", IErr, &
          " in allocation RWeightingCoefficients"
     RETURN
  ENDIF
  
  SELECT CASE (IWeightingFLAG)
  CASE(0)
     RWeightingCoefficients = ONE
  CASE(1)
     RWeightingCoefficientsDummy = RgVecMag(IOutputReflections)/MAXVAL(RgVecMag(IOutputReflections))
     IF(SIZE(RWeightingCoefficients).GT.1) THEN
        RWeightingCoefficientsDummy(1) = RWeightingCoefficients(2)/TWO 
     END IF
     DO ind = 1,IReflectOut
        RWeightingCoefficients(ind) = RWeightingCoefficientsDummy(IReflectOut-(ind-1))
     END DO
!!$     RWeightingCoefficients = 1/RgVecMag(IOutputReflections)
  CASE(2)
     RWeightingCoefficients = RgVecMag(IOutputReflections)/MAXVAL(RgVecMag(IOutputReflections))
     IF(SIZE(RWeightingCoefficients).GT.1) THEN
        RWeightingCoefficients(1) = RWeightingCoefficients(2)/TWO 
     END IF
  END SELECT

END SUBROUTINE InitialiseWeightingCoefficients

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

  INTEGER(IKIND) :: &
      IErr
  REAL(RKIND),INTENT(INOUT) :: &
       RStandardDeviation,RMean
  REAL(RKIND),INTENT(IN) :: &
       RFigureofMerit
  
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
