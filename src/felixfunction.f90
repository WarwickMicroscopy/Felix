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
SUBROUTINE FelixFunction(RIndependentVariableValues,IErr)

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
 
  CALL StructureFactorSetup(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error in StructureFactorSetup()"
     RETURN
  ENDIF
  
  Deallocate( &
       RgMatMat,RgMatMag,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation of RgMat"
     RETURN
  ENDIF

  CALL RankSymmetryRelatedStructureFactor(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error in StructureFactorRefinementSetup()"
     RETURN
  ENDIF
         
  CALL StructureFactorRefinementSetup(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error in StructureFactorRefinementSetup()"
     RETURN
  ENDIF
  
  CALL UpdateStructureFactors(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"FelixFunction(", my_rank, ") error ", IErr, &
          " in UpdateStructureFactors"
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
  
  IF(IWriteFLAG.GE.10) THEN
     
     PRINT*,"REDUCING Reflections",my_rank
     
  END IF

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
     CALL MPI_GATHERV(RIndividualReflections,ICount,&
          MPI_DOUBLE_PRECISION,RIndividualReflectionsRoot,&
          ICount,IDisplacements,MPI_DOUBLE_PRECISION,0,&
          MPI_COMM_WORLD,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " In MPI_GATHERV"
        RETURN
     ENDIF     
  ELSE     
     CALL MPI_GATHERV(CAmplitudeandPhase,ICount,&
          MPI_DOUBLE_COMPLEX,CAmplitudeandPhaseRoot,&
          ICount,IDisplacements,MPI_DOUBLE_COMPLEX,0, &
          MPI_COMM_WORLD,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " In MPI_GATHERV"
        RETURN
     ENDIF   
  END IF

  IF(IWriteFLAG.GE.10) THEN
     PRINT*,"REDUCED Reflections",my_rank
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
       CUgMatPrime,STAT=IErr)
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
          " Deallocating IPixelLocations"
     RETURN
  ENDIF

  DEALLOCATE( &
       IWeakBeamList,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating IPixelLocations"
     RETURN
  ENDIF

  DEALLOCATE( &
       IDisplacements,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating IPixelLocations"
     RETURN
  ENDIF

  DEALLOCATE( &
       ICount,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " Deallocating IPixelLocations"
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
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       SMNP, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       RFullAtomicFracCoordVec, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       SFullAtomicNameVec,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       IFullAnisotropicDWFTensor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       IFullAtomNumber,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       RFullIsotropicDebyeWallerFactor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       RFullPartialOccupancy,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       RDWF,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       ROcc,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       IAtoms,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       IAnisoDWFT,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       Rhkl,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       RgVecMag,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       RGn,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       ISymmetryRelations,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       ISymmetryStrengthKey,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  DEALLOCATE( &
       CSymmetryStrengthKey,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
          " in Deallocation"
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

SUBROUTINE CalculateFigureofMeritandDetermineThickness(RSimulatedImages,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       ind,jnd,knd,IErr,ICountedPixels,IThickness,&
       IThicknessCountFinal
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RSimulatedImageForPhaseCorrelation
  REAL(RKIND) :: &
       RCrossCorrelationOld,RThickness
  REAL(RKIND),DIMENSION(IReflectOut,IThicknessCount,IPixelTotal) :: &
       RSimulatedImages
  
  
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
  
  RCrossCorrelationOld = ZERO
  RThickness = ZERO

  DO ind = 1,IThicknessCount
     ICountedPixels = 0
     RSimulatedImageForPhaseCorrelation = ZERO
        RCrossCorrelation = ZERO
     DO jnd = 1,2*IPixelCount
        DO knd = 1,2*IPixelCount
           IF(ABS(RMask(jnd,knd)).GT.TINY) THEN
              ICountedPixels = ICountedPixels+1
              
              
              RSimulatedImageForPhaseCorrelation(jnd,knd) = &
                   RSimulatedImages(1,ind,ICountedPixels)
           END IF
        END DO
     END DO
     
     WHERE(RSimulatedImageForPhaseCorrelation.LT.TINY) RImageExpi(:,:,1) = ZERO
     WHERE(RImageExpi(:,:,1).LT.TINY) RSimulatedImageForPhaseCorrelation = ZERO
     

!!$        PRINT*,"---------------------------------------------------------"
!!$     PRINT*,RSimulatedImageForPhaseCorrelation(64:65,64)
!!$     PRINT*,RSimulatedImageForPhaseCorrelation(64:65,65)
!!$     PRINT*,RImageExpi(64:65,64,1)
!!$     PRINT*,RImageExpi(64:65,65,1)
!!$
!!$        PRINT*,"---------------------------------------------------------"


     CALL PhaseCorrelate(RSimulatedImageForPhaseCorrelation,RImageExpi(:,:,1),&
          IErr,2*IPixelCount,2*IPixelCount)

     PRINT*,RInitialThickness + (ind-1)*RDeltaThickness,&
          SUM(RSimulatedImageForPhaseCorrelation),&
          SUM(RImageExpi(:,:,1)),&
          MAXVAL(RSimulatedImageForPhaseCorrelation),&
          RCrossCorrelation
     
!!$     PRINT*,"RCrossCorrelation = ",RCrossCorrelation

     IF(RCrossCorrelation.GT.RCrossCorrelationOld) THEN
        RCrossCorrelationOld = RCrossCorrelation
        RThickness = RInitialThickness + (ind-1)*RDeltaThickness 
        IThicknessCountFinal = ind
     END IF
  END DO
  
  RCrossCorrelation = RCrossCorrelationOld

  DEALLOCATE(&
       RSimulatedImageForPhaseCorrelation,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CalculateFigureofMeritandDetermineThickness(", my_rank, ") error ", IErr, &
          " In DEALLOCATE"
     RETURN
  ENDIF
  
  
END SUBROUTINE CalculateFigureofMeritandDetermineThickness

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REAL(RKIND) FUNCTION SimplexFunction(RIndependentVariableValues,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       IErr,IPixels,CountPixels,ind,INewVariableID
  REAL(RKIND) :: &
       RNewValue
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(IN) :: &
       RIndependentVariableValues


!!$  RValue = RNewValue
!!$  IVariableID = INewVariableID
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"SimplexFunction(",my_rank,")"
  END IF

  IPixels = 0
  IPixels = CountPixels(IErr) ! Count The Number of Pixels
  IThicknessCount= (RFinalThickness- RInitialThickness)/RDeltaThickness + 1

  CALL UpdateVariables(RIndependentVariableValues,IErr) 
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
          " in UpdateVariables"
     RETURN
  ENDIF   

  CALL FelixFunction(RIndependentVariableValues,IErr) ! Simulate !!  
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

     CALL CalculateFigureofMeritandDetermineThickness(RIndividualReflections,IErr)
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
     
     DEALLOCATE( &
          IPixelLocations,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " Deallocating IPixelLocations"
        RETURN
     ENDIF

     DEALLOCATE( &
          RIndividualReflections,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " Deallocating IPixelLocations"
        RETURN
     ENDIF
     
     SimplexFunction = RCrossCorrelation
  END IF
  

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

!!$  RIndependentVariableValues(IVariableID) = RValue

  !!$  Fill the Independent Value array with values
  
  DO ind = 1,IIndependentVariables
     IVariableType = IIterativeVariableUniqueIDs(ind,2)
     SELECT CASE (IVariableType)
     CASE(1)
     CASE(2)
        RAtomSiteFracCoordVec(&
             IIterativeVariableUniqueIDs(ind,3),&
             IIterativeVariableUniqueIDs(ind,4)) = &
             RIndependentVariableValues(ind)
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

  IIterationCount = 2
  
  IF(IRefineModeSelectionArray(1).EQ.1.AND.IIterationCount.NE.1) THEN
     DO ind = 1,INoofUgs
        CSymmetryStrengthKey(ind) = &
             CMPLX(RIndependentVariableValues(ind),CKIND)
     ENd DO
  END IF
  
END SUBROUTINE UpdateStructureFactors
