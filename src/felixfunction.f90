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
!!$  IIterationCount
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       IDisplacements,ICount
  LOGICAL,INTENT(IN) :: &
       LInitialSimulationFLAG !If function is being called during initialisation
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RIndividualReflectionsRoot,&
       RFinalMontageImageRoot
  COMPLEX(CKIND),DIMENSION(:,:,:), ALLOCATABLE :: &
       CAmplitudeandPhaseRoot 
!!$  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(INOUT) :: &
!!$       RIndependentVariableValues 

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
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
  ENDIF
  
  
  !--------------------------------------------------------------------
  ! Setup Image
  !--------------------------------------------------------------------


  CALL ImageSetup( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error in ImageSetup()"
     RETURN
  ENDIF

 
  !--------------------------------------------------------------------
  ! MAIN section
  !--------------------------------------------------------------------
 

  IF(IAbsorbFLAG.NE.0) THEN
     
     IAbsorbFLAG = 0
     
     CALL StructureFactorSetup(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ",IErr,"in StructureFactorSetup()"
        RETURN
     ENDIF
     
     IAbsorbFLAG = 1
  ELSE
     
     CALL StructureFactorSetup(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ",IErr,"in StructureFactorSetup()"
        RETURN
     ENDIF

  END IF
  
  IF(IRefineModeSelectionArray(1).EQ.1.AND.LInitialSimulationFLAG.NEQV..TRUE.) THEN
     
     CALL ApplyNewStructureFactors(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in ApplyNewStructureFactors()"
        RETURN
     ENDIF
     
  END IF

  IF(IAbsorbFLAG.NE.0) THEN
     
     ALLOCATE(&
          CUgMatPrime(nReflections,nReflections),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables CUgMatPrime"
        RETURN
     ENDIF
     
     CALL StructureFactorsWithAbsorptionDetermination(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ",IErr,&
             "in StructureFactorsWithAbsorptionDetermination()"
        RETURN
     ENDIF

     CUgMat = CUgMat + CUgMatPrime

     DEALLOCATE( &
          CUgMatPrime,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " Deallocating CUgMatPrime"
        RETURN
     ENDIF
     
  END IF     
  !--------------------------------------------------------------------
  ! reserve memory for effective eigenvalue problem
  !--------------------------------------------------------------------

  !Kprime Vectors and Deviation Parameter
  
  ALLOCATE( &
       RDevPara(nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RDevPara"
     RETURN
  ENDIF

  ALLOCATE( &
       IStrongBeamList(nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables IStrongBeamList"
     RETURN
  ENDIF

  ALLOCATE( &
       IWeakBeamList(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables IWeakBeamList"
     RETURN
  ENDIF

  !--------------------------------------------------------------------
  ! MAIN LOOP: solve for each (ind,jnd) pixel
  !--------------------------------------------------------------------

  ILocalPixelCountMin= (IPixelTotal*(my_rank)/p)+1
  ILocalPixelCountMax= (IPixelTotal*(my_rank+1)/p) 

  IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Felixfunction(", my_rank, "): starting the eigenvalue problem"
     PRINT*,"Felixfunction(", my_rank, "): for lines ", ILocalPixelCountMin, &
          " to ", ILocalPixelCountMax
  ENDIF
  
  IThicknessCount= (RFinalThickness- RInitialThickness)/RDeltaThickness + 1

  IF(IImageFLAG.LE.2) THEN
     ALLOCATE( &
          RIndividualReflections(IReflectOut,IThicknessCount,&
          (ILocalPixelCountMax-ILocalPixelCountMin)+1),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Individual Images"
        RETURN
     ENDIF
     
     RIndividualReflections = ZERO
  ELSE
     ALLOCATE( &
          CAmplitudeandPhase(IReflectOut,IThicknessCount,&
          (ILocalPixelCountMax-ILocalPixelCountMin)+1),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Amplitude and Phase"
        RETURN
     ENDIF
     CAmplitudeandPhase = CZERO
  END IF

  ALLOCATE( &
       CFullWaveFunctions(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CFullWaveFunctions"
     RETURN
  ENDIF
  
  ALLOCATE( &
       RFullWaveIntensity(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RFullWaveIntensity"
     RETURN
  ENDIF  

  IMAXCBuffer = 200000
  IPixelComputed= 0
  
  IF(IWriteFLAG.GE.0.AND.my_rank.EQ.0) THEN
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
     ENDIF
  END DO
  
  IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Felixfunction : ",my_rank," is exiting calculation loop"
  END IF

  !--------------------------------------------------------------------
  ! close outfiles
  !--------------------------------------------------------------------
  
  ALLOCATE( &
       RIndividualReflectionsRoot(IReflectOut,IThicknessCount,IPixelTotal),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Root Reflections"
     RETURN
  ENDIF
  
  IF(IImageFLAG.GE.3) THEN
     ALLOCATE(&
          CAmplitudeandPhaseRoot(IReflectOut,IThicknessCount,IPixelTotal),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Root Amplitude and Phase"
        RETURN
     ENDIF
     CAmplitudeandPhaseRoot = CZERO
  END IF

  RIndividualReflectionsRoot = ZERO

  ALLOCATE(&
       IDisplacements(p),ICount(p),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " In ALLOCATE"
     RETURN
  ENDIF

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
     ENDIF     
  ELSE     
     CALL MPI_GATHERV(CAmplitudeandPhase,SIZE(CAmplitudeandPhase),&
          MPI_DOUBLE_COMPLEX,CAmplitudeandPhaseRoot,&
          ICount,IDisplacements,MPI_DOUBLE_COMPLEX,0, &
          MPI_COMM_WORLD,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " In MPI_GATHERV"
        RETURN
     ENDIF   
  END IF
  
  IF(IImageFLAG.GE.3) THEN
     DEALLOCATE(&
          CAmplitudeandPhase,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " Deallocating CAmplitudePhase"
        RETURN
     ENDIF   
  END IF
   
  IF(IImageFLAG.LE.2) THEN
     DEALLOCATE( &
          RIndividualReflections,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " Deallocating RIndividualReflections"
        RETURN
     ENDIF   
  END IF
  
  IF(my_rank.EQ.0.AND.IImageFLAG.GE.3) THEN
     RIndividualReflectionsRoot = &
          CAmplitudeandPhaseRoot * CONJG(CAmplitudeandPhaseRoot)
  END IF
  
  IF(my_rank.EQ.0) THEN
     ALLOCATE( &
          RIndividualReflections(IReflectOut,IThicknessCount,IPixelTotal),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Root Reflections"
        RETURN
     ENDIF
     
     RIndividualReflections = RIndividualReflectionsRoot
  END IF
  
  !--------------------------------------------------------------------
  ! free memory
  !--------------------------------------------------------------------
  
  !Dellocate Global Variables
  
  DEALLOCATE( &
       RgVecMatT,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating RgVecMatT"
     RETURN
  ENDIF

  DEALLOCATE( &
       RMask,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error  ", IErr, &
          " in Deallocation of RMask etc"
     RETURN
  ENDIF

  DEALLOCATE( &
       RDevPara,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating RDevPara"
     RETURN
  ENDIF

  DEALLOCATE( &
       IPixelLocations,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating IPixelLocations"
     RETURN
  ENDIF

  DEALLOCATE( &
       IStrongBeamList,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating IStrongBeamList"
     RETURN
  ENDIF

  DEALLOCATE( &
       IWeakBeamList,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating IWeakBeamList"
     RETURN
  ENDIF

  DEALLOCATE( &
       IDisplacements,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating IDisplacements"
     RETURN
  ENDIF

  DEALLOCATE( &
       ICount,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating ICount"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       CFullWaveFunctions, & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CFullWaveFunctions"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       RFullWaveIntensity, & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RFullWaveIntensity"
     RETURN
  ENDIF  

  IF(IImageFLAG.GE.3) THEN
     DEALLOCATE(&
          CAmplitudeandPhaseRoot,STAT=IErr) 
     
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error in Deallocation of CAmplitudeandPhase"
        RETURN  
     ENDIF
  END IF

  DEALLOCATE(&
       RIndividualReflectionsRoot,STAT=IErr) 
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error in Deallocation of RIndividualReflectionsRoot "
     RETURN  
  ENDIF
  
  DEALLOCATE( &
       MNP,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation MNP"
     RETURN
  ENDIF

  DEALLOCATE( &
       SMNP, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation SMNP"
     RETURN
  ENDIF

  DEALLOCATE( &
       RFullAtomicFracCoordVec, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RFullAtomicFracCoordVec"
     RETURN
  ENDIF

  DEALLOCATE( &
       SFullAtomicNameVec,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation SFullAtomicNameVec"
     RETURN
  ENDIF

  DEALLOCATE( &
       IFullAnisotropicDWFTensor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation IFullAnisotropicDWFTensor"
     RETURN
  ENDIF

  DEALLOCATE( &
       IFullAtomNumber,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation IFullAtomNumber"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       RFullIsotropicDebyeWallerFactor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RFullIsotropicDebyeWallerFactor"
     RETURN
  ENDIF

  DEALLOCATE( &
       RFullPartialOccupancy,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RFullPartialOccupancy"
     RETURN
  ENDIF

  DEALLOCATE( &
       RDWF,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RDWF"
     RETURN
  ENDIF

  DEALLOCATE( &
       ROcc,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation ROcc"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       IAtoms,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation IAtoms"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       IAnisoDWFT,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation IAnisoDWFT"
     RETURN
  ENDIF

  IF(LInitialSimulationFLAG.NEQV..TRUE.) THEN !Need these for Simplex Initialisation
     
     DEALLOCATE( &
          RgVecMag,&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in Deallocation RgVecMag"
        RETURN
     ENDIF
     
     DEALLOCATE( &
          CUgMat,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " Deallocating CUgMat"
        RETURN
     ENDIF
     
  END IF

  IF((my_rank.NE.0).AND.(LInitialSimulationFLAG.NEQV..TRUE.)) THEN     
     DEALLOCATE( &
          Rhkl,&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in Deallocation Rhkl"
        RETURN
     ENDIF     
  END IF

!!$  IF (my_rank.NE.0) THEN

      DEALLOCATE( &
           Rhklpositions,STAT=IErr)
      IF( IErr.NE.0 ) THEN
         PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
              " Deallocating Rhklpositions"
         RETURN
      ENDIF
!!$   END IF

  DEALLOCATE( &
       RgVecVec,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecVec"
     RETURN
  ENDIF

END SUBROUTINE FelixFunction

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

  INTEGER(IKIND) :: &
       ind,jnd,knd,IErr,ICountedPixels,IThickness,hnd
  INTEGER(IKIND),DIMENSION(IReflectOut) :: &
       IThicknessByReflection
  INTEGER(IKIND),INTENT(OUT) :: &
       IThicknessCountFinal
  REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: &
       RSimulatedImageForPhaseCorrelation,RExperimentalImage
  REAL(RKIND) :: &
       RCrossCorrelationOld,RIndependentCrossCorrelation,RThickness,PhaseCorrelate
!!$  REAL(RKIND),DIMENSION(IReflectOut,IThicknessCount,IPixelTotal) :: &
!!$       RSimulatedImages
  REAL(RKIND),DIMENSION(IReflectOut) :: &
       RReflectionCrossCorrelations
  REAL(RKIND) :: &
       ResidualSumofSquares
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
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
        
        
        IF(ICorrelationFLAG.EQ.0) THEN
           
            RIndependentCrossCorrelation = &
                 1.0_RKIND-&
                 PhaseCorrelate(&
                 RSimulatedImageForPhaseCorrelation,RExperimentalImage,&
                 IErr,2*IPixelCount,2*IPixelCount)

        ELSE
           RIndependentCrossCorrelation = &
                ResidualSumofSquares(&
                RSimulatedImageForPhaseCorrelation,RImageExpi(:,:,hnd),IErr)
        END IF
                
        IF(ABS(RIndependentCrossCorrelation).LT.RCrossCorrelationOld) THEN

           RCrossCorrelationOld = RIndependentCrossCorrelation

           IThicknessByReflection(hnd) = ind

        END IF
     END DO
     
     RReflectionCrossCorrelations(hnd) = RCrossCorrelationOld
     
  END DO

  RCrossCorrelation = &
       SUM(RReflectionCrossCorrelations*RWeightingCoefficients)/&
       REAL(IReflectOut,RKIND)
  
  IThicknessCountFinal = SUM(IThicknessByReflection)/IReflectOut

  RThickness = RInitialThickness + (IThicknessCountFinal-1)*RDeltaThickness 

  IF(my_rank.eq.0) THEN
     PRINT*,"Thicknesses",IThicknessByReflection
     PRINT*,"Correlation",RCrossCorrelation
     PRINT*,"Thickness Final",IThicknessCountFinal
     PRINT*,"Thickness",RThickness
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
  
  INTEGER(IKIND) :: &
       IErr,IExitFLAG,IThickness
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(INOUT) :: &
       RIndependentVariableValues
  INTEGER(IKIND),INTENT(IN) :: &
       IIterationCount
  LOGICAL :: &
       LInitialSimulationFLAG = .FALSE.
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"SimplexFunction(",my_rank,")"
  END IF
  

  IF(IRefineModeSelectionArray(1).EQ.1) THEN     
     CALL UpdateStructureFactors(RIndependentVariableValues,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in UpdateStructureFactors"
        RETURN
     ENDIF     
  ELSE
     CALL UpdateVariables(RIndependentVariableValues,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in UpdateVariables"
        RETURN
     ENDIF
     WHERE(RAtomSiteFracCoordVec.LT.0) RAtomSiteFracCoordVec=RAtomSiteFracCoordVec+ONE
     WHERE(RAtomSiteFracCoordVec.GT.1) RAtomSiteFracCoordVec=RAtomSiteFracCoordVec-ONE
  END IF

  IF (my_rank.EQ.0) THEN
     CALL PrintVariables(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in PrintVariables"
        RETURN
     ENDIF
  END IF

  CALL FelixFunction(LInitialSimulationFLAG,IErr) ! Simulate !!  
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
          " in FelixFunction"
     RETURN
  ENDIF

  IF(my_rank.EQ.0) THEN   
     CALL CreateImagesAndWriteOutput(IIterationCount,IExitFLAG,IErr) 
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in CreateImagesAndWriteOutput"
        RETURN
     ENDIF
     SimplexFunction = RCrossCorrelation     
  END IF
    
END FUNCTION SimplexFunction

SUBROUTINE CreateImagesAndWriteOutput(IIterationCount,IExitFLAG,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       IErr,IThicknessIndex,IIterationCount,IExitFLAG
  
  ALLOCATE( &
       RMask(2*IPixelCount,2*IPixelCount),&
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
  
  DEALLOCATE( &
       RMask,&
       STAT=IErr)
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
  
  DEALLOCATE( &
       RIndividualReflections,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateImagesAndWriteOutput(", my_rank, ") error ", IErr, &
          " Deallocating RIndividualReflections"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       IPixelLocations,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateImagesAndWriteOutput(", my_rank, ") error ", IErr, &
          " Deallocating IPixelLocations"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       Rhkl,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateImagesAndWriteOutput(", my_rank, ") error ", IErr, &
          " in Deallocation Rhkl"
     RETURN
  ENDIF
END SUBROUTINE CreateImagesAndWriteOutput

SUBROUTINE WriteIterationOutput(IIterationCount,IThicknessIndex,IExitFlag,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,IIterationCount,IThickness
  INTEGER(IKIND),INTENT(IN) :: &
       IThicknessIndex,IExitFLAG
  CHARACTER*200 :: &
       path
  
  IF(IExitFLAG.EQ.1.OR.(IIterationCount.GE.(IPreviousPrintedIteration+IPrint))) THEN
     

     IThickness = RInitialThickness + (IThicknessIndex-1)*RDeltaThickness 
     
     
     WRITE(path,"(A2,I1.1,I1.1,I1.1,I1.1,A2,I5.5,A2,I5.5,A2,I5.5,A10,I5.5)") &
          "f-",&
          IScatterFactorMethodFLAG, &
          IZolzFLAG, &
          IAbsorbFLAG, &
          IAnisoDebyeWallerFactorFlag,&
          "-T",IThickness,&
          "-P",2*IPixelcount,&
          "-P",2*IPixelcount,&
          "_Iteration",IIterationCount
     
     call system('mkdir ' // path)
     
     PRINT*,"I am Printing Because IExitFLAG = ",IExitFLAG,"and im",&
          IIterationCount-IPreviousPrintedIteration,"Iterations from my last print"
     
     IPreviousPrintedIteration = IIterationCount
     
     CALL WriteIterationImages(path,IThicknessIndex,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationOutput(", my_rank, ") error ", IErr, &
             " in WriteIterationImages"
        RETURN
     ENDIF
     
     DEALLOCATE( &
          Rhkl,&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationOutput(", my_rank, ") error ", IErr, &
             " in Deallocation Rhkl"
        RETURN
     ENDIF

     CALL WriteIterationStructure(path,IErr) 
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationOutput(", my_rank, ") error ", IErr, &
             " in WriteIterationStructure"
        RETURN
     ENDIF
      
  ELSE

  END IF
     
END SUBROUTINE WriteIterationOutput

SUBROUTINE WriteIterationStructure(path,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,jnd
  CHARACTER*200,INTENT(IN) :: &
       path
  CHARACTER*200 :: &
       SPrintString,filename,fullpath

!!$  Write out non symmetrically related atomic positions

  WRITE(filename,*) "StructureCif.txt"
  WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))

  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
        FILE=TRIM(ADJUSTL(fullpath)))
  
  DO jnd = 1,SIZE(RAtomSiteFracCoordVec,DIM=1)
     WRITE(IChOutSimplex,FMT='(A2,1X,3(F9.6,1X))') SAtomName(jnd),RAtomSiteFracCoordVec(jnd,:)
  END DO
  
  CLOSE(IChOutSimplex)

!!$  Write out full atomic positions

  CALL ExperimentalSetup(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in ExperimentalSetup"
     RETURN
  ENDIF

  WRITE(filename,*) "StructureFull.txt"
  WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))

  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
        FILE=TRIM(ADJUSTL(fullpath)))
  
  DO jnd = 1,SIZE(MNP,DIM=1)
     WRITE(IChOutSimplex,FMT='(A2,1X,3(F9.6,1X))') SMNP(jnd),MNP(jnd,1:3)
  END DO
  
  CLOSE(IChOutSimplex)
  
  DEALLOCATE( &
       MNP,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation MNP"
     RETURN
  ENDIF
      
  DEALLOCATE( &
        SMNP, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation SMNP"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RFullAtomicFracCoordVec, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RFullAtomicFracCoordVec"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       SFullAtomicNameVec,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation SFullAtomicNameVec"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IFullAnisotropicDWFTensor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation IFullAnisotropicDWFTensor"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IFullAtomNumber,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation IFullAtomNumber"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RFullIsotropicDebyeWallerFactor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RFullIsotropicDebyeWallerFactor"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RFullPartialOccupancy,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RFullPartialOccupancy"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RDWF,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RDWF"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       ROcc,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation ROcc"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IAtoms,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation IAtoms"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IAnisoDWFT,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation IAnisoDWFT"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RgVecMatT,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMatT"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RgVecMag,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMag"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RgVecVec,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMag"
     RETURN
  ENDIF

  DEALLOCATE( &
       RrVecMat, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, " in DEALLOCATE of RrVecMat"
     RETURN
  ENDIF

END SUBROUTINE WriteIterationStructure

SUBROUTINE WriteIterationImages(path,IThicknessIndex,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind,jnd,hnd,gnd,IThicknessIndex
  REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: &
       RImage
  CHARACTER*200,INTENT(IN) :: &
       path

  DO ind = 1,IReflectOut
     CALL OpenReflectionImage(IChOutWIImage,path,IErr,ind,2*IPixelCount,2_IKIND)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationImages(", my_rank, ") error in OpenReflectionImage()"
        RETURN
     ENDIF
     
     RImage = ZERO
     DO jnd = 1,IPixelTotal
        gnd = IPixelLocations(jnd,1)
        hnd = IPixelLocations(jnd,2)
        RImage(gnd,hnd) = RIndividualReflections(ind,IThicknessIndex,jnd)
     END DO
     
     CALL WriteReflectionImage(IChOutWIImage,&
          RImage,IErr,2*IPixelCount,2*IPixelCount,2_IKIND)       
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationImages(", my_rank, ") error in WriteReflectionImage()"
        RETURN
     ENDIF
     
     CLOSE(IChOutWIImage,IOSTAT=IErr)       
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationImages(", my_rank, ") error Closing Reflection Image()"
        RETURN
     ENDIF
     
  END DO

END SUBROUTINE WriteIterationImages

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

  INTEGER(IKIND) :: &
       IVariableType,IErr,ind
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(IN) :: &
       RIndependentVariableValues

  !!$  Fill the Independent Value array with values
  
  
  IF(IRefineModeSelectionArray(2).EQ.1) THEN     
     RAtomSiteFracCoordVec = RInitialAtomSiteFracCoordVec
  END IF

  DO ind = 1,IIndependentVariables
     IVariableType = IIterativeVariableUniqueIDs(ind,2)
     SELECT CASE (IVariableType)
     CASE(1)
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

  INTEGER(IKIND) :: &
       IErr,ind,IVariableType,jnd,knd
  REAL(RKIND),DIMENSION(3) :: &
       RCrystalVector
  CHARACTER*200 :: &
       SPrintString

  RCrystalVector = [RLengthX,RLengthY,RLengthZ]

  DO ind = 1,IRefinementVariableTypes
     IF (IRefineModeSelectionArray(ind).EQ.1) THEN
        SELECT CASE(ind)
        CASE(1)
           PRINT*,"Current Structure Factors"
           DO jnd = 1,INoofUgs
              PRINT*,CSymmetryStrengthKey(jnd)
           END DO           
        CASE(2)
           PRINT*,"Current Atomic Coordinates"
           DO jnd = 1,SIZE(RAtomSiteFracCoordVec,DIM=1)
              WRITE(SPrintString,FMT='(3(F9.3,1X))') RAtomSiteFracCoordVec(jnd,:)*RCrystalVector
              PRINT*,TRIM(ADJUSTL(SPrintString))              
           END DO
        CASE(3)
           PRINT*,"Current Atomic Occupancy"
           DO jnd = 1,SIZE(RAtomicSitePartialOccupancy,DIM=1)
              WRITE(SPrintString,FMT='((F9.6,1X))') RAtomicSitePartialOccupancy(jnd)
              PRINT*,TRIM(ADJUSTL(SPrintString))
           END DO
        CASE(4)
           PRINT*,"Current Isotropic Debye Waller Factors"
           DO jnd = 1,SIZE(RIsotropicDebyeWallerFactors,DIM=1)
              WRITE(SPrintString,FMT='((F9.6,1X))') RIsotropicDebyeWallerFactors(jnd)
              PRINT*,TRIM(ADJUSTL(SPrintString))
           END DO
        CASE(5)
           PRINT*,"Current Anisotropic Debye Waller Factors"
           DO jnd = 1,SIZE(RAnisotropicDebyeWallerFactorTensor,DIM=1)
              DO knd = 1,3
                 WRITE(SPrintString,FMT='((F9.6,1X))') RAnisotropicDebyeWallerFactorTensor(jnd,knd,:)
              PRINT*,TRIM(ADJUSTL(SPrintString))
              END DO
           END DO
        CASE(6)
           PRINT*,"Current Unit Cell Parameters"
           WRITE(SPrintString,FMT='(3(F9.6,1X))') RLengthX,RLengthY,RLengthZ
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
  
  INTEGER(IKIND) :: &
       IErr,ind
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(IN) :: &
       RIndependentVariableValues

  IF(IRefineModeSelectionArray(1).EQ.1) THEN
     DO ind = 1,INoofUgs
        CSymmetryStrengthKey(ind) = &
             CMPLX(RIndependentVariableValues((ind-1)*2+1),RIndependentVariableValues((ind-1)*2+2),CKIND)
     END DO
  END IF
  
END SUBROUTINE UpdateStructureFactors

SUBROUTINE ConvertVectorMovementsIntoAtomicCoordinates(IVariableID,RIndependentVariableValues,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind,jnd,IVariableID,IVectorID,IAtomID
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

SUBROUTINE InitialiseWeightingCoefficients(IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       RWeightingCoefficientsDummy

  ALLOCATE( &
       RWeightingCoefficients(IReflectOut),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"InitialiseWeightingCoefficients(", my_rank, ") error ", IErr, &
          " in allocation RWeightingCoefficients"
     RETURN
  ENDIF
  ALLOCATE( &
       RWeightingCoefficientsDummy(IReflectOut),&
       STAT=IErr)
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

  DEALLOCATE(&
       RgVecMag,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"InitialiseWeightingCoefficients(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMag"
     RETURN
  ENDIF


END SUBROUTINE InitialiseWeightingCoefficients

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
