!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: FelixSim.f90,v 1.89 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PROGRAM felixsim
 
  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  !--------------------------------------------------------------------
  ! local variable definitions
  !--------------------------------------------------------------------
  
  IMPLICIT NONE

  REAL(RKIND) :: &
       RThickness,&
       Duration, &
       time, norm
  INTEGER(IKIND) :: &
       ind,jnd,hnd,knd,pnd,gnd,IErr, &
       IHours,IMinutes,ISeconds,IMilliSeconds,&
       IThicknessIndex, &
       ILocalPixelCountMin, ILocalPixelCountMax
  INTEGER :: IStartTime, ICurrentTime ,IRate
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       IDisplacements,ICount
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RIndividualReflectionsRoot,&
       RFinalMontageImageRoot
  COMPLEX(CKIND),DIMENSION(:,:,:), ALLOCATABLE :: &
       CAmplitudeandPhaseRoot

  CHARACTER*40 surname, my_rank_string 
  CHARACTER*1000  SLocalPixelCountMin, SLocalPixelCountMax


  !-------------------------------------------------------------------
  ! constants
  !-------------------------------------------------------------------

  CALL Init_Numbers
  
  !-------------------------------------------------------------------
  ! set the error value to zero, will change upon error
  !-------------------------------------------------------------------

  IErr=0

  !--------------------------------------------------------------------
  ! MPI initialization
  !--------------------------------------------------------------------

  ! Initialise MPI  
  CALL MPI_Init(IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error in MPI_Init()"
     GOTO 9999
  ENDIF

  ! Get the rank of the current process
  CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error in MPI_Comm_rank()"
     GOTO 9999
  ENDIF

  ! Get the size of the current communicator
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error in MPI_Comm_size()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! protocal feature startup
  !--------------------------------------------------------------------
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"--------------------------------------------------------------"
     PRINT*,"felixsim: ", RStr
     PRINT*,"          ", DStr
     PRINT*,"          ", AStr
     PRINT*,"          on rank= ", my_rank, " of ", p, " in total."
     PRINT*,"--------------------------------------------------------------"
  END IF  
    
  !--------------------------------------------------------------------
  ! timing startup
  !--------------------------------------------------------------------

  CALL SYSTEM_CLOCK(count_rate=IRate)
  CALL SYSTEM_CLOCK(IStarttime)

  !--------------------------------------------------------------------
  ! INPUT section 
  !--------------------------------------------------------------------
  
  ISoftwareMode = 0 ! felixsimmode
  
  !Read from input files
  CALL ReadInput (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error in ReadInput()"
     GOTO 9999
  ENDIF
  
  CALL Message("felixsim",IMust,IErr)   
  CALL Message("felixsim",IInfo,IErr, MessageVariable = "ITotalAtoms", &
       IVariable = ITotalAtoms)

  !-------------------------------------------------------------------- 
  !Setup Experimental Variables
  !--------------------------------------------------------------------

  CALL ExperimentalSetup (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error in ExperimentalSetup()"
     GOTO 9999
  ENDIF


  !--------------------------------------------------------------------
  ! Setup Image
  !--------------------------------------------------------------------

  CALL ImageSetup( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error in ImageSetup()"
     GOTO 9999
  ENDIF

 
  !--------------------------------------------------------------------
  ! MAIN section
  !--------------------------------------------------------------------
 
  CALL StructureFactorSetup(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error in StructureFactorSetup()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! reserve memory for effective eigenvalue problem
  !--------------------------------------------------------------------

  !Kprime Vectors and Deviation Parameter
  
  ALLOCATE( &
       RDevPara(nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RDevPara"
     GOTO 9999
  ENDIF

  ALLOCATE( &
       IStrongBeamList(nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables IStrongBeamList"
     GOTO 9999
  ENDIF

  ALLOCATE( &
       IWeakBeamList(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables IWeakBeamList"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! MAIN LOOP: solve for each (ind,jnd) pixel
  !--------------------------------------------------------------------

  ILocalPixelCountMin= (IPixelTotal*(my_rank)/p)+1
  ILocalPixelCountMax= (IPixelTotal*(my_rank+1)/p)
  
  WRITE(SLocalPixelCountMin,"(I6.1)")ILocalPixelCountMin
  WRITE(SLocalPixelCountMax,"(I6.1)")ILocalPixelCountMax


  CALL Message("felixsim",IAllInfo,IErr,MessageString=": starting the eigenvalue problem")
  CALL Message("felixsim",IAllInfo,IErr,MessageString="for lines " // &
       TRIM(ADJUSTL(SLocalPixelCountMin)) // " to "// TRIM(ADJUSTL(SLocalPixelCountMax)))
       

!!$  IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
!!$     PRINT*,"felixsim(", my_rank, "): starting the eigenvalue problem"
!!$     PRINT*,"felixsim(", my_rank, "): for lines ", ILocalPixelCountMin, &
!!$          " to ", ILocalPixelCountMax
!!$  ENDIF
  

  IThicknessCount= (RFinalThickness- RInitialThickness)/RDeltaThickness + 1

  IF(IImageFLAG.LE.2) THEN
     ALLOCATE( &
          RIndividualReflections(IReflectOut,IThicknessCount,&
          (ILocalPixelCountMax-ILocalPixelCountMin)+1),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Individual Images"
        GOTO 9999
     ENDIF
     
     RIndividualReflections = ZERO
  ELSE
     ALLOCATE( &
          CAmplitudeandPhase(IReflectOut,IThicknessCount,&
          (ILocalPixelCountMax-ILocalPixelCountMin)+1),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixsim(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Amplitude and Phase"
        GOTO 9999
     ENDIF
     CAmplitudeandPhase = CZERO
  END IF

  ALLOCATE( &
       CFullWaveFunctions(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CFullWaveFunctions"
     GOTO 9999
  ENDIF
  
  ALLOCATE( &
       RFullWaveIntensity(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RFullWaveIntensity"
     GOTO 9999
  ENDIF  

  IMAXCBuffer = 200000
  IPixelComputed= 0
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0.AND.ISoftwareMode.LT.2) &
       .OR.IWriteFLAG.GE.10.AND.ISoftwareMode .LT. 2) THEN

     PRINT*,"*********************************"
     CALL Message("felixsim",ISilent,IErr,MessageString = " Entering BlochLoop")   
     PRINT*,"*********************************"
     
  END IF
  
  DO knd = ILocalPixelCountMin,ILocalPixelCountMax,1
     ind = IPixelLocations(knd,2)
     jnd = IPixelLocations(knd,1)
     CALL BlochCoefficientCalculation(ind,jnd,knd,ILocalPixelCountMin,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixsim(", my_rank, ") error ", IErr, &
             " in BlochCofficientCalculation"
        GOTO 9999
     ENDIF
  END DO
  
!!$     reset message counter
  IMessageCounter = 0

  CALL Message("felixsim",IAllInfo,IErr,&
       MessageString="is exiting calculation loop")

  !IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
  !   PRINT*,"felixsim : ",my_rank," is exiting calculation loop"
  !END IF

  ALLOCATE( &
       RIndividualReflectionsRoot(IReflectOut,IThicknessCount,IPixelTotal),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Root Reflections"
     GOTO 9999
  ENDIF
  
  IF(IImageFLAG.GE.3) THEN
     ALLOCATE(&
          CAmplitudeandPhaseRoot(IReflectOut,IThicknessCount,IPixelTotal),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixsim(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Root Amplitude and Phase"
        GOTO 9999
     ENDIF
     CAmplitudeandPhaseRoot = CZERO
  END IF

  RIndividualReflectionsRoot = ZERO

  ALLOCATE(&
       IDisplacements(p),ICount(p),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " In ALLOCATE"
     GOTO 9999
  ENDIF

  DO pnd = 1,p
     IDisplacements(pnd) = (IPixelTotal*(pnd-1)/p)*IReflectOut*IThicknessCount
     ICount(pnd) = (((IPixelTotal*(pnd)/p) - (IPixelTotal*(pnd-1)/p)))*IReflectOut*IThicknessCount
  END DO 

  IF(IImageFLAG.LE.2) THEN
     CALL MPI_GATHERV(RIndividualReflections,SIZE(RIndividualReflections),&
          MPI_DOUBLE_PRECISION,RIndividualReflectionsRoot,&
          ICount,IDisplacements,MPI_DOUBLE_PRECISION,0,&
          MPI_COMM_WORLD,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixsim(", my_rank, ") error ", IErr, &
             " In MPI_GATHERV"
        GOTO 9999
     ENDIF     
  ELSE     
     CALL MPI_GATHERV(CAmplitudeandPhase,SIZE(CAmplitudeandPhase),&
          MPI_DOUBLE_COMPLEX,CAmplitudeandPhaseRoot,&
          ICount,IDisplacements,MPI_DOUBLE_COMPLEX,0, &
          MPI_COMM_WORLD,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixsim(", my_rank, ") error ", IErr, &
             " In MPI_GATHERV"
        GOTO 9999
     ENDIF   
  END IF

  IF(IImageFLAG.GE.3) THEN
     DEALLOCATE(&
          CAmplitudeandPhase,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixsim(", my_rank, ") error ", IErr, &
             " Deallocating CAmplitudePhase"
        GOTO 9999
     ENDIF   
  END IF
   
  IF(IImageFLAG.LE.2) THEN
     DEALLOCATE( &
          RIndividualReflections,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixsim(", my_rank, ") error ", IErr, &
             " Deallocating RIndividualReflections"
        GOTO 9999
     ENDIF   
  END IF

  IF(my_rank.EQ.0) THEN
     ALLOCATE( &
          RFinalMontageImageRoot(MAXVAL(IImageSizeXY),&
          MAXVAL(IImageSizeXY),IThicknessCount),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixsim(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Root Montage"
        GOTO 9999
     ENDIF

     RFinalMontageImageRoot = ZERO		
  END IF

  IF(my_rank.EQ.0.AND.IImageFLAG.GE.3) THEN
     RIndividualReflectionsRoot = &
          CAmplitudeandPhaseRoot * CONJG(CAmplitudeandPhaseRoot)
  END IF

  IF(my_rank.EQ.0) THEN
     CALL MontageSetup(RFinalMontageImageRoot, &
          RIndividualReflectionsRoot,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixsim(", my_rank, ") error ", IErr, &
             " in MontageSetup"
        GOTO 9999
     ENDIF
  END IF
  !--------------------------------------------------------------------
  ! Write out Images
  !--------------------------------------------------------------------

  IF (my_rank.EQ.0) THEN

     CALL WriteOutput(CAmplitudeandPhaseRoot,RIndividualReflectionsRoot,RFinalMontageImageRoot,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixsim(", my_rank, ") error ", IErr, &
             " in WriteOutput"
        GOTO 9999
     ENDIF
          
  END IF
  
  !--------------------------------------------------------------------
  ! free memory
  !--------------------------------------------------------------------
  
  !Dellocate Global Variables  

  DEALLOCATE( &
       RgVecMatT,STAT=IErr)
  IF( IErr.NE.0 ) THEN

     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " Deallocating RgVecMatT"
     GOTO 9999

  ENDIF
  
  DEALLOCATE( &
       Rhklpositions,STAT=IErr)
  IF( IErr.NE.0 ) THEN

     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " Deallocating Rhklpositions"
     GOTO 9999
  ENDIF 

  DEALLOCATE( &
       RMask,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error  ", IErr, &
          " in Deallocation of RMask etc"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       CUgMat,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " Deallocating CUgMat"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       CUgMatNoAbs,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " Deallocating CUgMatNoAbs"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       CUgMatPrime,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " Deallocating CUgMat"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       RDevPara,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " Deallocating RDevPara"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       IPixelLocations,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " Deallocating IPixelLocations"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       IStrongBeamList,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " Deallocating IStrongBeamList"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       IWeakBeamList,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " Deallocating IWeakBeamList"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       IDisplacements,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " Deallocating IDisplacements"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       ICount,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " Deallocating ICount"
     GOTO 9999

  ENDIF
  
  DEALLOCATE( &
       CFullWaveFunctions, & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CFullWaveFunctions"
     GOTO 9999
  ENDIF
  
  DEALLOCATE( &
       RFullWaveIntensity, & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RFullWaveIntensity"
     GOTO 9999
  ENDIF  

  IF(IImageFLAG.GE.3) THEN
     DEALLOCATE(&
          CAmplitudeandPhaseRoot,STAT=IErr) 
     
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixsim(", my_rank, ") error in Deallocation of CAmplitudeandPhaseRoot"
        GOTO 9999  
     ENDIF
  END IF

  DEALLOCATE(&
       RIndividualReflectionsRoot,STAT=IErr) 
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error in Deallocation of RIndividualReflectionsRoot "
     GOTO 9999  
  ENDIF
  
  DEALLOCATE( &
       MNP,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation MNP"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       SMNP, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation SMNP"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       RFullAtomicFracCoordVec, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation RFullAtomicFracCoordVec"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       SFullAtomicNameVec,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation SFullAtomicNameVec"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       IFullAnisotropicDWFTensor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation IFullAnisotropicDWFTensor"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       IFullAtomNumber,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation IFullAtomNumber"
     GOTO 9999
  ENDIF
  
  DEALLOCATE( &
       RFullIsotropicDebyeWallerFactor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation RFullIsotropicDebyeWallerFactor"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       RFullPartialOccupancy,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation RFullPartialOccupancy"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       RDWF,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation RDWF"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       ROcc,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation ROcc"
     GOTO 9999
  ENDIF
  
  DEALLOCATE( &
       IAtoms,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation IAtoms"
     GOTO 9999
  ENDIF
  
  DEALLOCATE( &
       IAnisoDWFT,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation IAnisoDWFT"
     GOTO 9999
  ENDIF
  
  DEALLOCATE( &
       Rhkl,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation Rhkl"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       RgVecMag,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMag"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       RgVecVec,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecVec"
     GOTO 9999
  ENDIF
     
  !--------------------------------------------------------------------
  ! finish off
  !--------------------------------------------------------------------

  WRITE(my_rank_string,*) my_rank
    
  CALL SYSTEM_CLOCK(ICurrentTime)
  Duration=REAL(ICurrentTime-IStartTime)/REAL(IRate)
  IHours = FLOOR(Duration/3600.0D0)
  IMinutes = FLOOR(MOD(Duration,3600.0D0)/60.0D0)
  ISeconds = MOD(Duration,3600.0D0)-IMinutes*60
  IMilliSeconds = INT((Duration-(IHours*3600+IMinutes*60+ISeconds))*1000,IKIND)

  PRINT*, "felixsim( ", TRIM(ADJUSTL(my_rank_string)), " ) ", &
       RStr, ", used time=", IHours, "hrs ", &
       IMinutes,"mins ",ISeconds,"secs ", IMilliSeconds,"millisecs"

  CALL MPI_Barrier(MPI_COMM_WORLD,IErr)

  !--------------------------------------------------------------------
  ! Shut down MPI
  !--------------------------------------------------------------------
9999 &
  CALL MPI_Finalize(IErr)
  IF( IErr.NE.0 ) THEN
     
     PRINT*,"felixsim(", my_rank, ") error ", IErr, " in MPI_Finalize()"
     STOP
  ENDIF

  ! clean shutdown
  STOP

END PROGRAM felixsim
