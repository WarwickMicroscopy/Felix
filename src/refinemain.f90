!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! FelixSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
!
! (C) 2013/14, all right reserved
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This file is part of FelixSim.
!
!  FelixSim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  FelixSim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with FelixSim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PROGRAM FelixRefine
 
  USE MyNumbers
  
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
       time, norm,StartTime, CurrentTime, Duration, TotalDurationEstimate
  INTEGER(IKIND) :: &
       IErr,ind,jnd,hnd,knd,gnd, &
       IHours,IMinutes,ISeconds,IX,IY
  CHARACTER*34 :: &
       filename
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RImageSim,RImageExpi
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RImageExpiCentred
  COMPLEX(CKIND) sumC, sumD
  REAL(RKIND) :: RCrossCorrelationOld
  REAL(RKIND) :: RTestCase
  !REAL(RKIND) StartTime, CurrentTime, Duration, TotalDurationEstimate
  
  !--------------------------------------------------------------------
  ! image related variables	
  REAL(RKIND) Rx0,Ry0, RImageRadius,Rradius, Rthickness
  
  INTEGER(IKIND) ILocalPixelCountMin, ILocalPixelCountMax
  COMPLEX(RKIND) CVgij
  
  !INTEGER ind,jnd,hnd,knd,pnd
  INTEGER, DIMENSION(:), ALLOCATABLE :: &
       IWeakBeamVec,IRankArraySize,IRankArraySizeRoot
  REAL(RKIND),DIMENSION(:,:,:,:),ALLOCATABLE :: &
       RIndividualReflectionsRoot
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RFinalMontageImageRoot
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RImage
  COMPLEX(CKIND),DIMENSION(:,:,:,:), ALLOCATABLE :: &
       CAmplitudeandPhaseRoot
  INTEGER IRootArraySize, IPixelPerRank
  CHARACTER*40 surname, path
  CHARACTER*25 CThickness 
  CHARACTER*25 CThicknessLength
 
  INTEGER(IKIND),DIMENSION(2,2) :: ITest
  
  INTEGER(IKIND):: IThickness, IThicknessIndex, ILowerLimit, &
       IUpperLimit
  !REAL(RKIND) StartTime, CurrentTime, Duration, TotalDurationEstimate

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
     PRINT*,"main(", my_rank, ") error in MPI_Init()"
     GOTO 9999
  ENDIF

  ! Get the rank of the current process
  CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in MPI_Comm_rank()"
     GOTO 9999
  ENDIF

  ! Get the size of the current communicator
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in MPI_Comm_size()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! protocal feature startup
  !--------------------------------------------------------------------
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"FelixSim: ", RStr,DStr,AStr, ", process ", my_rank, " of ", p
     PRINT*,"--------------------------------------------------------------"
  END IF

  !--------------------------------------------------------------------
  ! timing startup
  !--------------------------------------------------------------------

  CALL cpu_time(StartTime)

  !--------------------------------------------------------------------
  ! INPUT section
  !--------------------------------------------------------------------

  IImageSizeXY(1) = 128
  IImageSizeXY(2) = 128
  
  
  ALLOCATE( &
       RImageSim(IImageSizeXY(1),IImageSizeXY(2)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Refinemain (", my_rank, ") error in Allocation()"
     GOTO 9999
  ENDIF
  
  ALLOCATE( &
       RImageExpi(IImageSizeXY(1),IImageSizeXY(2)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Refinemain (", my_rank, ") error in Allocation()"
     GOTO 9999
  ENDIF
  
  DO ind = 1,1
     
     WRITE(filename,"(A6,I3.3,A4)") "Felix.",ind,".img"

     CALL OpenImageForReadIn(IErr,filename)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"Refinemain (", my_rank, ") error in OpenImageForReadIn()"
        GOTO 9999
     END IF
     
     ALLOCATE( &
          RImageIn(IImageSizeXY(1),IImageSizeXY(2)), &
          STAT=IErr)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"Refinemain (", my_rank, ") error in Allocation()"
        GOTO 9999
     ENDIF
     
     CALL ReadImageForRefinement(IErr)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"Refinemain (", my_rank, ") error in ReadImageForRefinement()"
        GOTO 9999
     ELSE
        IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"Image Read In Successful"
        END IF
     ENDIF

     RImageExpi = RImageIn

     DEALLOCATE( &
          RImageIn, &
          STAT=IErr)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"Refinemain (", my_rank, ") error in deAllocation()"
        GOTO 9999
     ENDIF
     
     CLOSE(IChInImage,IOSTAT=IErr) 
     IF( IErr.NE.0 ) THEN
        PRINT*,"Refinemain (", my_rank, ") error in CLOSE()"
        GOTO 9999
     END IF

  END DO

  CALL Input( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in Input()"
     GOTO 9999
  ENDIF

  CALL InputScatteringFactors( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in InputScatteringFactors()"
     GOTO 9999
  ENDIF

  CALL InpCIF(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in InpCIF()"
     GOTO 9999
  ENDIF

  IF (ITotalAtoms.EQ.0) THEN
     CALL CountTotalAtoms(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error in CountTotalAtoms()"
        GOTO 9999
     ENDIF
  END IF
     
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ITotalAtoms = ",ITotalAtoms
  END IF
  
  !--------------------------------------------------------------------
  ! open outfiles
  !--------------------------------------------------------------------

  WRITE(surname,'(A1,I1.1,A1,I1.1,A1,I1.1,A2,I4.4)') &
       "S", IScatterFactorMethodFLAG, &
       "B", ICentralBeamFLAG, &
       "M", IMaskFLAG, &
       "_P", IPixelCount
  
  ! eigensystem
  IF(IOutputFLAG.GE.1) THEN
     CALL OpenData_MPI(IChOutES_MPI, "ES", surname, IErr)
  ENDIF
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in OpenData_MPI(EigenSystem)"
     GOTO 9999
  ENDIF
  
  ! UgMatEffective
  IF(IOutputFLAG.GE.2) THEN
     CALL OpenData_MPI(IChOutUM_MPI, "UM", surname, IErr)
  ENDIF
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in OpenDataMPI()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! Allocate Crystallography Variables
  !--------------------------------------------------------------------
       
  ALLOCATE( &
       RrVecMat(ITotalAtoms,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, " in ALLOCATE()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! microscopy settings
  !--------------------------------------------------------------------

  CALL MicroscopySettings( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in MicroscopySettings()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! crystallography settings
  !--------------------------------------------------------------------

  CALL Crystallography( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in Crystallography()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! diffraction initialization
  !--------------------------------------------------------------------

  CALL DiffractionPatternDefinitions( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in DiffractionPatternDefinitions()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! allocate memory for DYNAMIC variables according to nReflections
  !--------------------------------------------------------------------

  ! Image initialisation 
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: nReflections=", nReflections
  END IF
  
  ALLOCATE( &
       Rhklpositions(nReflections,2), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variable Rhklpositions"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! image initialization
  !--------------------------------------------------------------------

  CALL ImageInitialization( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in ImageInitializtion()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! define image masks
  !--------------------------------------------------------------------
      
  !Allocate Memory for Masking Image

  ALLOCATE( &
       RMask(2*IPixelCount,2*IPixelCount),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variable RMask"
     GOTO 9999
  ENDIF

  CALL ImageMask(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ImageMask"
     GOTO 9999
  END IF

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: main(", my_rank, ") IPixelTotal=", IPixelTotal
  END IF
 
  !--------------------------------------------------------------------
  ! MAIN section
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! Calculate Reflection Matrix
  !--------------------------------------------------------------------

  ALLOCATE( &  
       RgMatMat(nReflections,nReflections,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     GOTO 9999
  ENDIF
       
  ALLOCATE( &  
       RgMatMag(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     GOTO 9999
  ENDIF

  CALL GMatrixInitialisation (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in GMatrixInitialisation"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! calculating Ug matrix
  !--------------------------------------------------------------------

  !Allocate memory for Ug Matrix

  ALLOCATE( & 
       CUgMat(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     GOTO 9999
  ENDIF       

  ALLOCATE( & 
       CUgMatPrime(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     GOTO 9999
  ENDIF       

  CALL UgCalculation (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in UgCalculation"
     GOTO 9999
  ENDIF
  
  ALLOCATE( &  
       ISymmetryRelations(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     GOTO 9999
  ENDIF

  CALL DetermineSymmetryRelatedUgs (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in DetermineSymmetryRelatedUgs"
     GOTO 9999
  ENDIF

  IF(IDevFLAG.EQ.1) THEN
     WHERE (ISymmetryRelations.EQ.INoofUgs) 
        CUgMat = CUgMat*(ONE+RPercentageUgChange/100.0_RKIND)
     END WHERE
  END IF

  CALL UgAddAbsorption(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in UgCalculation"
     GOTO 9999
  ENDIF


  IF(IAbsorbFLAG.EQ.1) THEN
     CUgMat =  CUgMat+CUgMatPrime
  end IF

  Deallocate( &
       RgMatMat,RgMatMag,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in Deallocation of RgMat"
     GOTO 9999
  ENDIF
  
  !--------------------------------------------------------------------
  ! high-energy approximation (not HOLZ compatible)
  !--------------------------------------------------------------------
  
  RBigK= SQRT(RElectronWaveVectorMagnitude**2 + RMeanInnerCrystalPotential)

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: main(", my_rank, ") BigK=", RBigK
  END IF
  
         
  !--------------------------------------------------------------------
  ! reserve memory for effective eigenvalue problem
  !--------------------------------------------------------------------

  !Kprime Vectors and Deviation Parameter
  
  ALLOCATE( &
       RDevPara(nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RDevPara"
     GOTO 9999
  ENDIF

  ALLOCATE( &
       IStrongBeamList(nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables IStrongBeamList"
     GOTO 9999
  ENDIF

  ALLOCATE( &
       IWeakBeamList(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables IWeakBeamList"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! MAIN LOOP: solve for each (ind,jnd) pixel
  !--------------------------------------------------------------------

  ILocalPixelCountMin= (IPixelTotal*(my_rank)/p)+1
  ILocalPixelCountMax= (IPixelTotal*(my_rank+1)/p) 
  
  IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"main(", my_rank, "): starting the eigenvalue problem"
     PRINT*,"main(", my_rank, "): for lines ", ILocalPixelCountMin, &
          " to ", ILocalPixelCountMax
  ENDIF
  
  IThicknessCount= (RFinalThickness- RInitialThickness)/RDeltaThickness + 1
  
  IF(IImageFLAG.LE.1) THEN
     ALLOCATE( &
          RIndividualReflections(2*IPixelCount,&
          2*IPixelCount,IReflectOut,IThicknessCount),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Individual Images"
        GOTO 9999
     ENDIF
     
     RIndividualReflections = ZERO
  ELSE
     ALLOCATE( &
          CAmplitudeandPhase(2*IPixelCount,&
          2*IPixelCount,IReflectOut,IThicknessCount),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Amplitude and Phase"
        GOTO 9999
     ENDIF
     CAmplitudeandPhase = CZERO
  END IF
  
  ALLOCATE( &
       CFullWaveFunctions(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CFullWaveFunctions"
     GOTO 9999
  ENDIF
  
  ALLOCATE( &
       RFullWaveIntensity(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RFullWaveIntensity"
     GOTO 9999
  ENDIF
  
  IMAXCBuffer = 200000
  IPixelComputed= 0
  
  DEALLOCATE( &
       RScattFactors, &
       RrVecMat, Rsg, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in DEALLOCATE() "
     GOTO 9999
  ENDIF
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.6) THEN
     PRINT*,"main(",my_rank,") Entering BlochLoop()"
  END IF
  
  DO knd = ILocalPixelCountMin,ILocalPixelCountMax,1
     ind = IPixelLocations(knd,1)
     jnd = IPixelLocations(knd,2)
     CALL BlochCoefficientCalculation(ind,jnd,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " in BlochCofficientCalculation"
        GOTO 9999
     ENDIF
  END DO
  
  IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"MAIN : ",my_rank," is exiting calculation loop"
  END IF
  
  !--------------------------------------------------------------------
  ! close outfiles
  !--------------------------------------------------------------------
  
  ! eigensystem
  IF(IOutputFLAG.GE.1) THEN
     CALL MPI_FILE_CLOSE(IChOutES_MPI, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " Closing IChOutES"
        GOTO 9999
     ENDIF
  ENDIF
  
  ! UgMatEffective
  IF(IOutputFLAG.GE.2) THEN
     CALL MPI_FILE_CLOSE(IChOutUM_MPI, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " Closing IChOutUM"
        GOTO 9999
     ENDIF
  ENDIF
  
  ALLOCATE( &
       RIndividualReflectionsRoot(2*IPixelCount,&
       2*IPixelCount,IReflectOut,IThicknessCount),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Root Reflections"
     GOTO 9999
  ENDIF
  
  IF(IImageFLAG.GE.2) THEN
     ALLOCATE(&
          CAmplitudeandPhaseRoot(2*IPixelCount,&
          2*IPixelCount,IReflectOut,IThicknessCount),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Root Amplitude and Phase"
        GOTO 9999
     ENDIF
     CAmplitudeandPhaseRoot = CZERO
  END IF
  
  RIndividualReflectionsRoot = ZERO
  
  IRootArraySize = SIZE(RIndividualReflectionsRoot)
  
  IF(IWriteFLAG.GE.10) THEN
     
     PRINT*,"REDUCING Reflections",my_rank
     
  END IF
  
  IF(IImageFLAG.LE.1) THEN
     CALL MPI_REDUCE(RIndividualReflections,RIndividualReflectionsRoot,&
          IRootArraySize,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
          MPI_COMM_WORLD,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " In MPI_REDUCE"
        GOTO 9999
     ENDIF
     
     
  ELSE     
     CALL MPI_REDUCE(CAmplitudeandPhase,CAmplitudeandPhaseRoot,&
          IRootArraySize,MPI_DOUBLE_COMPLEX,MPI_SUM,0, &
          MPI_COMM_WORLD,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " In MPI_REDUCE"
        GOTO 9999
     ENDIF
  END IF
  

  IF(my_rank.eq.0) THEN
     PRINT*,RIndividualReflectionsRoot(32:33,32:33,1,1)
  END IF

  IF(IWriteFLAG.GE.10) THEN
     PRINT*,"REDUCED Reflections",my_rank
  END IF

  IF (my_rank.EQ.0) THEN
     IImageSizeXY(1) = 2*IPixelCount
     IImageSizeXY(2) = 2*IPixelCount
     RCrossCorrelationOld = ZERO

     DO ind = 1,IThicknessCount
        CALL PhaseCorrelate(RIndividualReflectionsRoot(:,:,1,ind),RImageExpi,IErr)

        IF(RCrossCorrelation.GT.RCrossCorrelationOld) THEN
           RCrossCorrelationOld = RCrossCorrelation
           RThickness = RInitialThickness + (ind-1)*RDeltaThickness 
        END IF
        
        
     END DO
     
     PRINT*,"Thickness = ",RThickness," Angstoms"
  END IF
  
  IF(IImageFLAG.LE.1) THEN
     DEALLOCATE( &
          RIndividualReflections,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " Deallocating RIndividualReflections"
        GOTO 9999
     ENDIF
  END IF 

  !Dellocate Global Variables
  
  DEALLOCATE( &
       RgVecMatT, &
       Rhklpositions, RMask,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in Deallocation of RgVecMatT etc"
     GOTO 9999
  ENDIF
  DEALLOCATE( &
       CUgMat,IPixelLocations, &
       RDevPara,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in Deallocation of CUgmat etc"
     GOTO 9999
  ENDIF
  
  DEALLOCATE( &
       RIndividualReflectionsRoot,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in Deallocation of RIndividualReflectionsRoot"
     GOTO 9999
  ENDIF
  
  IF(IImageFLAG.GE.2) THEN
     DEALLOCATE(&
          CAmplitudeandPhaseRoot,STAT=IErr) 
     
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error in Deallocation of CAmplitudeandPhase"
        GOTO 9999
     ENDIF
  END IF
  
  !--------------------------------------------------------------------
  ! finish off
  !--------------------------------------------------------------------
    
  CALL cpu_time(CurrentTime)
  Duration=(CurrentTime-StartTime)
  IHours = FLOOR(Duration/3600.0D0)
  IMinutes = FLOOR(MOD(Duration,3600.0D0)/60.0D0)
  ISeconds = MOD(Duration,3600.0D0)-IMinutes*60.0D0

  PRINT*, "FelixSim(", my_rank, ") ", RStr, ", used time=", IHours, "hrs ",IMinutes,"mins ",ISeconds,"Seconds "

  !--------------------------------------------------------------------
  ! Shut down MPI
  !--------------------------------------------------------------------
9999 &
  CALL MPI_Finalize(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, " in MPI_Finalize()"
     STOP
  ENDIF

  ! clean shutdown
  STOP
  

END PROGRAM FelixRefine
