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
SUBROUTINE FelixFunction(RIndependentVariableValues,IIterationCount,IErr)

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
       IIterationFLAG,IIterationCount
  INTEGER, DIMENSION(:), ALLOCATABLE :: &
       IDisplacements,ICount
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RIndividualReflectionsRoot,&
       RFinalMontageImageRoot
  COMPLEX(CKIND),DIMENSION(:,:,:), ALLOCATABLE :: &
       CAmplitudeandPhaseRoot 
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(INOUT) :: &
       RIndependentVariableValues 

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

  END IF

  CALL ApplyNewStructureFactors(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in ApplyNewStructureFactors()"
     RETURN
  ENDIF


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

  Deallocate( &
       RgMatMat,RgMatMag,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation of RgMat"
     RETURN
  ENDIF

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
  
  DEALLOCATE( &
       RrVecMat, Rsg, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in DEALLOCATE() "
     RETURN
  ENDIF
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.6) THEN
     PRINT*,"Felixfunction(",my_rank,") Entering BlochLoop()"
  END IF

  DO knd = ILocalPixelCountMin,ILocalPixelCountMax,1
     ind = IPixelLocations(knd,1)
     jnd = IPixelLocations(knd,2)
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
       Rhklpositions,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating Rhklpositions"
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
       CUgMat,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating CUgMat"
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

  DEALLOCATE( &
       RgVecMag,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMag"
     RETURN
  ENDIF

  DEALLOCATE( &
       RGn,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation RGn"
     RETURN
  ENDIF

END SUBROUTINE FelixFunction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ReadExperimentalImages(IErr)

 USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       ind,IErr
  CHARACTER*34 :: &
       filename

  DO ind = 1,IReflectOut
     
     WRITE(filename,"(A6,I3.3,A4)") "felix.",ind,".img"
     
     CALL OpenImageForReadIn(IErr,filename)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadExperimentalImages (", my_rank, ") error in OpenImageForReadIn()"
        RETURN
     END IF
     
     ALLOCATE( &
          RImageIn(2*IPixelCount,2*IPixelCount), &
          STAT=IErr)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadExperimentalImages (", my_rank, ") error in Allocation()"
        RETURN
     ENDIF
     
     CALL ReadImageForRefinement(IErr)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadExperimentalImages (", my_rank, ") error in ReadImageForRefinement()"
        RETURN
     ELSE
        IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"Image Read In Successful"
        END IF
     ENDIF
     
     RImageExpi(:,:,ind) = RImageIn
     
     DEALLOCATE( &
          RImageIn, &
          STAT=IErr)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadExperimentalImages (", my_rank, ") error in deAllocation()"
        RETURN
     ENDIF

     
     CLOSE(IChInImage,IOSTAT=IErr) 
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadExperimentalImages (", my_rank, ") error in CLOSE()"
        RETURN
     END IF

  END DO


END SUBROUTINE ReadExperimentalImages

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CalculateFigureofMeritandDetermineThickness(RSimulatedImages,IThicknessCountFinal,IErr)
  
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
  INTEGER(IKIND),INTENT(OUT) :: &
       IThicknessCountFinal
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RSimulatedImageForPhaseCorrelation
  REAL(RKIND) :: &
       RCrossCorrelationOld,RIndependentCrossCorrelation,RThickness
  REAL(RKIND),DIMENSION(IReflectOut,IThicknessCount,IPixelTotal) :: &
       RSimulatedImages
  REAL(RKIND),DIMENSION(IReflectOut) :: &
       RReflectionCrossCorrelations
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"CalculateFigureofMeritandDetermineThickness(",my_rank,")"
  END IF

  ALLOCATE(&
       RSimulatedImageForPhaseCorrelation(2*IPixelCount,2*IPixelCount),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CalculateFigureofMeritandDetermineThickness(", my_rank, ") error ", IErr, &
          " In ALLOCATE"
     RETURN
  ENDIF

  RReflectionCrossCorrelations = ZERO

  DO hnd = 1,IReflectOut
     RCrossCorrelationOld = ZERO
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
                      RSimulatedImages(hnd,ind,ICountedPixels)
              END IF
           END DO
        END DO

        CALL PhaseCorrelate(RSimulatedImageForPhaseCorrelation,RImageExpi(:,:,hnd),&
             IErr,2*IPixelCount,2*IPixelCount)
        RIndependentCrossCorrelation = RCrossCorrelation       
        
        IF(RIndependentCrossCorrelation.GT.RCrossCorrelationOld) THEN

           RCrossCorrelationOld = RIndependentCrossCorrelation

           RThickness = RInitialThickness + (ind-1)*RDeltaThickness 

           IThicknessCountFinal = ind

        END IF
     END DO
     
     RReflectionCrossCorrelations(hnd) = RCrossCorrelationOld

  END DO

  RCrossCorrelation = SUM(RReflectionCrossCorrelations*RWeightingCoefficients)/REAL(IReflectOut,RKIND)

  DEALLOCATE(&
       RSimulatedImageForPhaseCorrelation,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CalculateFigureofMeritandDetermineThickness(", my_rank, ") error ", IErr, &
          " In DEALLOCATE of RSimulatedImageForPhaseCorrelation"
     RETURN
  ENDIF  
  
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
       IErr,IPixels,CountPixels,ind,INewVariableID,IExitFLAG,gnd,knd,IThickness,hnd,jnd,&
       IThicknessIndex,IVectorID,IAtomID
  REAL(RKIND) :: &
       RNewValue
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(IN) :: &
       RIndependentVariableValues
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RImage
  INTEGER(IKIND),INTENT(IN) :: &
       IIterationCount
  CHARACTER*40 :: &
       path
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"SimplexFunction(",my_rank,")"
  END IF
  
  IFelixCount = IFelixCount + 1
  IPixels = 0
  IPixels = CountPixels(IErr) ! Count The Number of Pixels
  IThicknessCount= (RFinalThickness- RInitialThickness)/RDeltaThickness + 1
  
  CALL UpdateVariables(RIndependentVariableValues,IErr)
  IF( IErr.NE.0 ) THEN
     
     PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
          " in UpdateVariables"
     RETURN
  ENDIF
  
  CALL UpdateStructureFactors(RIndependentVariableValues,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
          " in UpdateStructureFactors"
     RETURN
  ENDIF

  CALL FelixFunction(RIndependentVariableValues,IIterationCount,IErr) ! Simulate !!  
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
          " in FelixFunction"
     RETURN
  ENDIF
  
  IF(my_rank.EQ.0) THEN    
     
     ALLOCATE( &
          RMask(2*IPixelCount,2*IPixelCount),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variable RMask"
        RETURN
     ENDIF
     
     CALL ImageMaskInitialisation(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in ImageMaskInitialisation"
        RETURN
     ENDIF
     
     CALL CalculateFigureofMeritandDetermineThickness(RIndividualReflections,IThicknessIndex,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             "Calling function CalculateFigureofMeritandDetermineThickness"
        RETURN
     ENDIF

     DEALLOCATE( &
          RMask,&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in DEALLOCATE() of DYNAMIC variable RMask"
        RETURN
     ENDIF
     
!!$     OUTPUT AN IMAGE -------------------------------------
     
     IF(IExitFLAG.EQ.1.OR.(IIterationCount.GE.(IPreviousPrintedIteration+IPrint))) THEN
        PRINT*,"I am Printing Because IExitFLAG = ",IExitFLAG,"and im",&
             IPreviousPrintedIteration+IPrint,"Iterations from my last print"
        IPreviousPrintedIteration = IIterationCount
        IThickness = RInitialThickness + (IThicknessIndex-1)*RDeltaThickness 

        ALLOCATE( &
             RImage(2*IPixelCount,2*IPixelCount), &
             STAT=IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"WriteOutput(", my_rank, ") error ", IErr, &
                " in ALLOCATE() of DYNAMIC variables RImage"
           RETURN
        ENDIF
                
        WRITE(path,"(A2,A1,I1.1,A2,I1.1,A2,I1.1,A2,I4.4,A2,I5.5,A10,I5.5)") &
             "F-",&
             "S", IScatterFactorMethodFLAG, &
             "_B", ICentralBeamFLAG, &
             "_M", IMaskFLAG, &
             "_P", IPixelCount, &
             "_T", IThickness, &
             "_Iteration",IIterationCount
        
        call system('mkdir ' // path)
        
        DO ind = 1,IReflectOut
           CALL OpenReflectionImage(IChOutWIImage,path,IErr,ind,2*IPixelCount,2_IKIND)
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error in OpenReflectionImage()"
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
              PRINT*,"WriteOutput(", my_rank, ") error in WriteReflectionImage()"
              RETURN
           ENDIF

           CLOSE(IChOutWIImage,IOSTAT=IErr)       
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error Closing Reflection Image()"
              RETURN
           ENDIF

        END DO
        
        DEALLOCATE( &
             RImage, &
             STAT=IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"WriteOutput(", my_rank, ") error ", IErr, &
                " in DEALLOCATE() of DYNAMIC variables RImage"
           RETURN
        ENDIF
     END IF
     
!!$     FINISH OUT PUTTING IMAGE --------------------------------
     
     DEALLOCATE( &
          RIndividualReflections,&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " Deallocating RIndividualReflections"
        RETURN
     ENDIF
     
     DEALLOCATE( &
          IPixelLocations,&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " Deallocating IPixelLocations"
        RETURN
     ENDIF
     
     SimplexFunction = 1.0_RKIND-RCrossCorrelation
     
  END IF
    
  DEALLOCATE( &
     Rhkl,&
     STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation Rhkl"
     RETURN
  ENDIF

END FUNCTION SimplexFunction

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
  
  RAtomSiteFracCoordVec = RInitialAtomSiteFracCoordVec

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
     END SELECT
  END DO
  
END SUBROUTINE UpdateVariables

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
