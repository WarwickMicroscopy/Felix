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

!!$REAL(RKIND) FUNCTION FelixFunction(VariableId,VariableValue)
SUBROUTINE FelixFunction(IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind,jnd,knd,pnd,&
       IThicknessIndex,ILocalPixelCountMin, ILocalPixelCountMax
  INTEGER, DIMENSION(:), ALLOCATABLE :: &
       IDisplacements,ICount
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RIndividualReflectionsRoot
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RFinalMontageImageRoot
  COMPLEX(CKIND),DIMENSION(:,:,:), ALLOCATABLE :: &
       CAmplitudeandPhaseRoot

  !--------------------------------------------------------------------
  ! local variable definitions
  !--------------------------------------------------------------------
  

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
       RScattFactors, &
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

  IF(my_rank.EQ.0) THEN
     ALLOCATE( &
          RFinalMontageImageRoot(MAXVAL(IImageSizeXY),&
          MAXVAL(IImageSizeXY),IThicknessCount),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Root Montage"
        RETURN
     ENDIF

    RFinalMontageImageRoot = ZERO		
  END IF

  IF(my_rank.EQ.0.AND.IImageFLAG.GE.3) THEN
     RIndividualReflectionsRoot = &
          CAmplitudeandPhaseRoot * CONJG(CAmplitudeandPhaseRoot)
  END IF

  IF(my_rank.EQ.ZERO) THEN
     CALL MontageSetup(IThicknessIndex,knd,ind,jnd,RFinalMontageImageRoot, &
          RIndividualReflectionsRoot,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in MontageSetup"
        RETURN
     ENDIF
  END IF
  !--------------------------------------------------------------------
  ! Write out Images
  !--------------------------------------------------------------------
  
  IF (my_rank.EQ.0) THEN

     CALL WriteOutput(CAmplitudeandPhaseRoot,RIndividualReflectionsRoot,RFinalMontageImageRoot,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ", IErr, &
             " in WriteOutput"
        RETURN
     ENDIF
          
  END IF
  
  !--------------------------------------------------------------------
  ! free memory
  !--------------------------------------------------------------------
  
  !Dellocate Global Variables
  
  DEALLOCATE( &
       RgVecMatT, &
       Rhklpositions, RMask,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error in Deallocation of RgVecMatT etc"
     RETURN
  ENDIF
  DEALLOCATE( &
       CUgMat,IPixelLocations, &
       RDevPara,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error in Deallocation of CUgmat etc"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       RIndividualReflectionsRoot,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixfunction(", my_rank, ") error in Deallocation of RIndividualReflectionsRoot"
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
END SUBROUTINE FelixFunction
