!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! FelixSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
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
       IErr,ind,jnd,hnd,knd,gnd,pnd,fnd, &
       IHours,IMinutes,ISeconds,IX,IY
  CHARACTER*34 :: &
       filename
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RImageSim,RImageExpi
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RImageExpiCentred
  !REAL(RKIND) StartTime, CurrentTime, Duration, TotalDurationEstimate
  
  !--------------------------------------------------------------------
  ! image related variables	
  REAL(RKIND) :: &
       Rthickness,RDeltaDebyeWallerFactorPerElement,RCrossCorrelationOld, &
       RIterationTolerance
  
  INTEGER(IKIND) ILocalPixelCountMin, ILocalPixelCountMax,IAtomNo,&
       ILocalFluxCountMin,ILocalFluxCountMax,ISubgroupNo,IFluxStepsPerSubgroup,&
       ISubgroups,ISubgroupRootProcessContainingAnswer,ISubgroupContainingAnswer, &
       ISubgroupSize

  INTEGER my_status(MPI_STATUS_SIZE)
  INTEGER(IKIND) ICorrelationMaximum,IRefinementIteration

  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       IWeakBeamVec,IDisplacements,ICount,IFluxIterationIndices, &
       IRanks
  REAL(IKIND), DIMENSION(:), ALLOCATABLE :: &
       RFinalDWFConfig,RFinalDebyeWallerFactorPerElement
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RIndividualReflectionsRoot,RBestCorrelationImage
  REAL(RKIND),DIMENSION(:,:,:,:),ALLOCATABLE :: &
       RReflectionImagesForPhaseCorrelation
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RFinalMontageImageRoot
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RImage,RFluxIterationCorrelations,RFluxIterationCorrelationsRoot
  COMPLEX(CKIND),DIMENSION(:,:,:), ALLOCATABLE :: &
       CAmplitudeandPhaseRoot
  COMPLEX(CKIND),DIMENSION(:,:), ALLOCATABLE :: &
       CZeroMat
  INTEGER(IKIND) IThicknessCountFinal,INoofDWFs,IRemainder,IMaxiter, &
       InDWFs,a,b,c
  CHARACTER*40 surname, path
  CHARACTER*25 CThickness 
  CHARACTER*25 CThicknessLength

  INTEGER(IKIND),DIMENSION(2) :: ILoc
  
  INTEGER(IKIND):: IThickness, IThicknessIndex, ILowerLimit, &
       IUpperLimit,IWorldGrp,INewGroup,Inewcomm,my_newrank,q
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
     PRINT*,"RefineMain(", my_rank, ") error in MPI_Init()"
     GOTO 9999
  ENDIF
  ! Get the rank of the current process
  CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in MPI_Comm_rank()"
     GOTO 9999
  ENDIF

  ! Get the size of the current communicator
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in MPI_Comm_size()"
     GOTO 9999
  ENDIF

  PRINT*,"MPI comm world is ",p," Processors in size"

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
  
  ISoftwareMode = 2 ! FelixRefineMode

  CALL Input( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in Input()"
     GOTO 9999
  ENDIF

  CALL ReadInHKLs(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in ReadInHKLs()"
     GOTO 9999
  ENDIF  
  
  SELECT CASE (IRefineModeFLAG)
  CASE(0)
     
     !DebyeWallerLoop
     PRINT*,IElementList
     PRINT*,
     
     INoofDWFs = NINT((RFinalDebyeWallerFactor-&
          RInitialDebyeWallerFactor)/&
          RDeltaDebyeWallerFactor +1)**SIZE(IElementList,DIM=1)
  CASE(1)
  CASE(2)
     INoofDWFs = 1
     RIterationTolerance = RDeltaThickness
  END SELECT
  
  !---------------------------------------------------------------------
  !CLEAN THIS CODE UP

  a = p
  b = INoofDWFs
  PRINT*,"No of DWFs = ",b
  
!!$  DO                    ! now we have a <= b
!!$     c = MOD(a, b)      !    compute c, the reminder
!!$     IF (c == 0) EXIT   !    if c is zero, we are done.  GCD = b
!!$     a = b              !    otherwise, b becomes a
!!$     b = c              !    and c becomes b
!!$  END DO                !    go back

  CALL GreatestCommonDivisor(p,INoofDWFs,ISubgroups)
  
  PRINT*,"No of DWFs = ",ISubgroups

  IF(my_rank.EQ.0) Then
     
     PRINT*,"Running ",INoofDWFs, "Debye Waller Factor Configurations On ",p," Cores,Using ",p/ISubgroups," Cores per Configuration"
  END IF
  ALLOCATE(&
       IRanks(p/ISubgroups),&
       STAT=IErr)
  
  CALL MPI_COMM_GROUP(MPI_COMM_WORLD,IWorldGrp,IErr)
  
  jnd = floor(REAL(my_rank/(p/ISubgroups)))
  DO ind = 1,SIZE(IRanks)
     IRanks(ind) = (ind-1)+(jnd)*SIZE(IRanks)
  end DO
  
  CALL MPI_GROUP_INCL(IWorldGrp,p/ISubgroups,IRanks,INewGroup,IErr)
  
  CALL MPI_COMM_CREATE(MPI_COMM_WORLD,INewGroup,Inewcomm,IErr)
  
  CALL MPI_Comm_rank(Inewcomm,my_newrank,IErr)
  
  CALL MPI_Comm_size(Inewcomm,ISubgroupSize,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in MPI_Comm_size()"
     GOTO 9999
  ENDIF
     
  PRINT*,"Hello i am ",my_rank," of ",p," in comm_world and ",my_newrank,"of ",ISubgroupSize," in new comm"
  !---------------------------------------------------------------------


  CALL InputScatteringFactors( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in InputScatteringFactors()"
     GOTO 9999
  ENDIF

  CALL InpCIF(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in InpCIF()"
     GOTO 9999
  ENDIF

  IF (ITotalAtoms.EQ.0) THEN
     CALL CountTotalAtoms(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error in CountTotalAtoms()"
        GOTO 9999
     ENDIF
  END IF
     
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ITotalAtoms = ",ITotalAtoms
  END IF
  
  
  IImageSizeXY(1) = 2*IPixelCount
  IImageSizeXY(2) = 2*IPixelCount
  
  
  ALLOCATE( &
       RImageSim(IImageSizeXY(1),IImageSizeXY(2), &
       IReflectOut),&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain (", my_rank, ") error in Allocation()"
     GOTO 9999
  ENDIF
  
  ALLOCATE( &
       RImageExpi(IImageSizeXY(1),IImageSizeXY(2), &
       IReflectOut),&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain (", my_rank, ") error in Allocation()"
     GOTO 9999
  ENDIF
  
  DO ind = 1,IReflectOut
     
     WRITE(filename,"(A6,I3.3,A4)") "Felix.",ind,".img"

     CALL OpenImageForReadIn(IErr,filename)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain (", my_rank, ") error in OpenImageForReadIn()"
        GOTO 9999
     END IF
     
     ALLOCATE( &
          RImageIn(IImageSizeXY(1),IImageSizeXY(2)), &
          STAT=IErr)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain (", my_rank, ") error in Allocation()"
        GOTO 9999
     ENDIF
     
     CALL ReadImageForRefinement(IErr)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain (", my_rank, ") error in ReadImageForRefinement()"
        GOTO 9999
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
        PRINT*,"RefineMain (", my_rank, ") error in deAllocation()"
        GOTO 9999
     ENDIF
     
     CLOSE(IChInImage,IOSTAT=IErr) 
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain (", my_rank, ") error in CLOSE()"
        GOTO 9999
     END IF

  END DO


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
     PRINT*,"RefineMain(", my_rank, ") error in OpenData_MPI(EigenSystem)"
     GOTO 9999
  ENDIF
  
  ! UgMatEffective
  IF(IOutputFLAG.GE.2) THEN
     CALL OpenData_MPI(IChOutUM_MPI, "UM", surname, IErr)
  ENDIF
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in OpenDataMPI()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! Allocate Crystallography Variables
  !--------------------------------------------------------------------
       
  ALLOCATE( &
       RrVecMat(ITotalAtoms,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error ", IErr, " in ALLOCATE()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! microscopy settings
  !--------------------------------------------------------------------

  CALL MicroscopySettings( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in MicroscopySettings()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! crystallography settings
  !--------------------------------------------------------------------

  CALL Crystallography( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in Crystallography()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! diffraction initialization
  !--------------------------------------------------------------------

  CALL DiffractionPatternDefinitions( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in DiffractionPatternDefinitions()"
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
     PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variable Rhklpositions"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! image initialization
  !--------------------------------------------------------------------

  CALL ImageInitialization( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in ImageInitializtion()"
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
     PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variable RMask"
     GOTO 9999
  ENDIF

  CALL ImageMask(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
          " in ImageMask"
     GOTO 9999
  END IF

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RefineMain(", my_rank, ") IPixelTotal=", IPixelTotal
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
     PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     GOTO 9999
  ENDIF
       
  ALLOCATE( &  
       RgMatMag(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     GOTO 9999
  ENDIF

  CALL GMatrixInitialisation (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
          " in GMatrixInitialisation"
     GOTO 9999
  ENDIF
  
  
  
  !Allocate memory for Ug Matrix
  
  ALLOCATE( & 
       CUgMat(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CUgMat"
     GOTO 9999
  ENDIF
  
  ALLOCATE( & 
       CZeroMat(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CZeroMat"
     GOTO 9999
  ENDIF
  
  ALLOCATE( & 
       CUgMatPrime(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CUgMatPrime"
     GOTO 9999
  ENDIF
  
  
  ALLOCATE( &
       RFinalDWFConfig(IElements),STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in allocation of RFinalDWFConfig "
     GOTO 9999
  END IF
  
  ALLOCATE( &
       RFinalDebyeWallerFactorPerElement(IElements),STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in allocation of RFinalDWFConfig "
     GOTO 9999
  END IF


  

  !--------------------------------------------------------------------
  ! Map iteration steps for flux loop
  !--------------------------------------------------------------------
  
  CALL DetermineFluxSteps(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
          " in DetermineFluxSteps"
     GOTO 9999
  ENDIF
  
  IRefinementIteration = 0
  
  DO 
     IRefinementIteration = IRefinementIteration + 1
     IF(my_rank.EQ.0) THEN
        PRINT*,"Entering Iteration",IRefinementIteration
     END IF
     
     SELECT CASE (IRefineModeFLAG)
     CASE(0)
        
        !DebyeWallerLoop
        
        IF(IRefinementIteration.GT.1) THEN
           RFinalDebyeWallerFactorPerElement = RFinalDWFConfig+RDeltaDebyeWallerFactorPerElement
           RFinalDebyeWallerFactorPerElement = RFinalDWFConfig-RDeltaDebyeWallerFactorPerElement
           RDeltaDebyeWallerFactor = RDeltaDebyeWallerFactor/10
        ELSE
           
        END IF
        
        INoofDWFs = ((RFinalDebyeWallerFactor-&
             RInitialDebyeWallerFactor)/&
             RDeltaDebyeWallerFactor +1)
     CASE(1)
        
     CASE(2)
        !Thickness Determination Doesnt Require the Flux Loop
     CASE DEFAULT     
        IErr = 1
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                " Refinement Mode Not Defined"
           GOTO 9999
        ENDIF
     END SELECT
     
     
     
     ALLOCATE(&
          IFluxIterationIndices(IElements),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
             " in ALLOCATE IFluxIterationIndices"
        GOTO 9999
     ENDIF
     
     IFluxIterationIndices = 0        
     
     ALLOCATE(&
          RFluxIterationCorrelations(IFluxIterationSteps,2),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
             " in ALLOCATE RFluxIterationCorrelations"
        GOTO 9999
     ENDIF
     
     RFluxIterationCorrelations = Zero
     
     ALLOCATE(&
          RFluxIterationCorrelationsRoot(IFluxIterationSteps,2),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
             " in ALLOCATE RFluxIterationCorrelationsRoot"
        GOTO 9999
     ENDIF
     
     RFluxIterationCorrelationsRoot = Zero
     
     ALLOCATE(&
          RBestCorrelationImage(2*IPixelCount,2*IPixelCount,IReflectOut),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
             " in ALLOCATE RBestCorrelationImage"
        GOTO 9999
     ENDIF
     
     RBestCorrelationImage = Zero
     
     !--------------------------------------------------------------------
     ! ENTER REFINEMENT LOOP
     !--------------------------------------------------------------------  
     
     ISubgroups = ISubgroups
     ISubgroupNo = FLOOR(REAL(my_rank,RKIND)/(REAL(p,RKIND)/REAL(ISubgroups,RKIND)))
     IFluxStepsPerSubgroup = IFluxIterationSteps/ISubgroups
!!$  PRINT*,"I am rank",my_rank," In, Subgroup",ISubgroupNo,"Of",ISubgroups
     
     IMaxiter = IElements
!!$  PRINT*,"Total Flux Steps",IFluxIterationSteps
     ILocalFluxCountMin= (ISubgroupNo*IFluxStepsPerSubgroup)+1
     ILocalFluxCountMax= (ISubgroupNo+1)*IFluxStepsPerSubgroup
     
     IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"RefineMain(", my_rank, "): starting the eigenvalue problem"
        PRINT*,"RefineMain(", my_rank, "): for lines ", ILocalFluxCountMin, &
             " to ", ILocalFluxCountMax
     ENDIF
     DO fnd = ILocalFluxCountMin,ILocalFluxCountMax

        CUgMat = CZERO ! Initialise Ugs at the beginning of each iteration
        
        SELECT CASE (IRefineModeFLAG)
           
        CASE(0)
           
           DO ind = 1,IMaxiter
              
              IF(ind.EQ.1) THEN
                 
                 IFluxIterationIndices(ind) = FLOOR(REAL((fnd-1))/REAL(INoofDWFs**(SIZE(IElementlist)-ind)))+1
                 IRemainder = MOD(fnd-1,INoofDWFs**(SIZE(IElementlist)-ind))
              ELSE
                 IF(ind.EQ.IMaxiter) THEN
                    IFluxIterationIndices(ind) = IRemainder+1             
                 ELSE
                    IFluxIterationIndices(ind) = FLOOR(REAL(IRemainder)/REAL(INoofDWFs**(SIZE(IElementlist)-ind)))+1
                    IRemainder = MOD(IRemainder,INoofDWFs**(SIZE(IElementlist)-ind))
                 END IF
              END IF
           END DO
           
           DO ind = 1,IMaxIter
              WHERE(IAtoms.EQ.IElementList(ind))
                 RDWF = ((IFluxIterationIndices(ind)-1)*RDeltaDebyeWallerFactor)+RInitialDebyeWallerFactor
              END WHERE
           END DO
           
           !--------------------------------------------------------------------
           ! calculating Ug matrix with variable DebyeWallerFactor
           !--------------------------------------------------------------------
           

           CALL UgCalculation (IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " in UgCalculation"
              GOTO 9999
           ENDIF
           
        CASE(1)
           
           !--------------------------------------------------------------------
           ! calculating Ug matrix for Ug Alteration
           !--------------------------------------------------------------------
           
           CALL UgCalculation (IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " in UgCalculation"
              GOTO 9999
           ENDIF
           
           CZeroMAT = CZERO
           
           DO ind = 1,nReflections
              DO jnd = 1,ind
                 CZeroMAT(ind,jnd) = CONE
              END DO
           END DO
           
           ALLOCATE( &  
                ISymmetryRelations(nReflections,nReflections), &
                STAT=IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
              GOTO 9999
           ENDIF
           
           ISymmetryRelations = ISymmetryRelations*CZeroMat  
           
           CALL DetermineSymmetryRelatedUgs (IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " in DetermineSymmetryRelatedUgs"
              GOTO 9999
           ENDIF
           
           DO ind = 1,(SIZE(ISymmetryStrengthKey,DIM=1))
              ILoc = MINLOC(ABS(ISymmetryRelations-ind))
              ISymmetryStrengthKey(ind,1) = ind
              ISymmetryStrengthKey(ind,2) = ind
              CSymmetryStrengthKey(ind) = CUgMat(ILoc(1),ILoc(2))
           END DO
           
           CALL ReSortUgs(ISymmetryStrengthKey,CSymmetryStrengthKey,SIZE(CSymmetryStrengthKey,DIM=1))
           
           IF(IDevFLAG.EQ.1) THEN
              CUgMat = CUgMat*CZeroMat
              WHERE (ISymmetryRelations.EQ.ISymmetryStrengthKey(INoofUgs,2)) 
                 CUgMat = CUgMat*(ONE+RPercentageUgChange/100.0_RKIND)
              END WHERE
              CUgMat = CUgmat + CONJG(TRANSPOSE(CUgMat))
           END IF
           
        CASE(2)
           
           !--------------------------------------------------------------------
           ! calculating Ug matrix for Thickness Determination
           !--------------------------------------------------------------------
           
           CALL UgCalculation (IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " in UgCalculation"
              GOTO 9999
           ENDIF
        CASE DEFAULT
           
           IErr = 1
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " Refinement Mode Not Defined"
              GOTO 9999
           ENDIF
           
        END SELECT
        
        CALL UgAddAbsorption(IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                " in UgCalculation"
           GOTO 9999
        ENDIF
        
        IF(IAbsorbFLAG.EQ.1) THEN
           CUgMat =  CUgMat+CUgMatPrime
        end IF

        !--------------------------------------------------------------------
        ! high-energy approximation (not HOLZ compatible)
        !--------------------------------------------------------------------
        
        RBigK= SQRT(RElectronWaveVectorMagnitude**2 + RMeanInnerCrystalPotential)
        
        IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"RefineMain(", my_rank, ") BigK=", RBigK
        END IF
        
        
        !--------------------------------------------------------------------
        ! reserve memory for effective eigenvalue problem
        !--------------------------------------------------------------------
        
        !Kprime Vectors and Deviation Parameter
        
        ALLOCATE( &
             RDevPara(nReflections), &
             STAT=IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                " in ALLOCATE() of DYNAMIC variables RDevPara"
           GOTO 9999
        ENDIF
        
        ALLOCATE( &
             IStrongBeamList(nReflections), &
             STAT=IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                " in ALLOCATE() of DYNAMIC variables IStrongBeamList"
           GOTO 9999
        ENDIF
        
        ALLOCATE( &
             IWeakBeamList(nReflections), & 
             STAT=IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                " in ALLOCATE() of DYNAMIC variables IWeakBeamList"
           GOTO 9999
        ENDIF
        
        !--------------------------------------------------------------------
        ! MAIN LOOP: solve for each (ind,jnd) pixel
        !--------------------------------------------------------------------
        
        ILocalPixelCountMin= (IPixelTotal*(my_newrank)/ISubgroupSize)+1
        ILocalPixelCountMax= (IPixelTotal*(my_newrank+1)/ISubgroupSize) 
        
        IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"RefineMain(", my_rank, "): starting the eigenvalue problem"
           PRINT*,"RefineMain(", my_rank, "): for lines ", ILocalPixelCountMin, &
                " to ", ILocalPixelCountMax
        ENDIF
        
        IThicknessCount= (RFinalThickness- RInitialThickness)/RDeltaThickness +1 
        
        IF(IImageFLAG.LE.2) THEN
           ALLOCATE( &
                RIndividualReflections(IReflectOut,IThicknessCount,&
                (ILocalPixelCountMax-ILocalPixelCountMin)+1),&
                STAT=IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
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
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " in ALLOCATE() of DYNAMIC variables Amplitude and Phase"
              GOTO 9999
           ENDIF
           CAmplitudeandPhase = CZERO
        END IF
        
        ALLOCATE( &
             CFullWaveFunctions(nReflections), & 
             STAT=IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                " in ALLOCATE() of DYNAMIC variables CFullWaveFunctions"
           GOTO 9999
        ENDIF
        
        ALLOCATE( &
             RFullWaveIntensity(nReflections), & 
             STAT=IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                " in ALLOCATE() of DYNAMIC variables RFullWaveIntensity"
           GOTO 9999
        ENDIF
        
        IMAXCBuffer = 200000
        IPixelComputed= 0
        
        IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.6) THEN
           PRINT*,"RefineMain(",my_rank,") Entering BlochLoop()"
        END IF
        
        DO knd = ILocalPixelCountMin,ILocalPixelCountMax,1
           ind = IPixelLocations(knd,1)
           jnd = IPixelLocations(knd,2)
           CALL BlochCoefficientCalculation(ind,jnd,knd,ILocalPixelCountMin,IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " in BlochCofficientCalculation"
              GOTO 9999
           ENDIF
        END DO
        
        PRINT*,"Exiting Bloch Loop"
        
        IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"REFINEMAIN : ",my_rank," is exiting calculation loop"
        END IF
        
        !--------------------------------------------------------------------
        ! close outfiles
        !--------------------------------------------------------------------
        
        ! eigensystem
        IF(IOutputFLAG.GE.1) THEN
           CALL MPI_FILE_CLOSE(IChOutES_MPI, IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " Closing IChOutES"
              GOTO 9999
           ENDIF
        ENDIF
        
        ! UgMatEffective
        IF(IOutputFLAG.GE.2) THEN
           CALL MPI_FILE_CLOSE(IChOutUM_MPI, IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " Closing IChOutUM"
              GOTO 9999
           ENDIF
        ENDIF
        
        ALLOCATE( &
             RIndividualReflectionsRoot(IReflectOut,IThicknessCount,IPixelTotal),&
             STAT=IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                " in ALLOCATE() of DYNAMIC variables Root Reflections"
           GOTO 9999
        ENDIF
        
        IF(IImageFLAG.GE.3) THEN
           ALLOCATE(&
                CAmplitudeandPhaseRoot(IReflectOut,IThicknessCount,IPixelTotal),&
                STAT=IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " in ALLOCATE() of DYNAMIC variables Root Amplitude and Phase"
              GOTO 9999
           ENDIF
           CAmplitudeandPhaseRoot = CZERO
        END IF
        
        RIndividualReflectionsRoot = ZERO
        
        IF(IWriteFLAG.GE.10) THEN
           
           PRINT*,"REDUCING Reflections",my_rank
           
        END IF
        
        IF(IWriteFLAG.GE.10) THEN
           PRINT*,"REDUCED Reflections",my_rank
        END IF
        
        ALLOCATE(&
             IDisplacements(ISubgroupSize),ICount(ISubgroupSize),&
             STAT=IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                " In ALLOCATE of IDisplacements"
           GOTO 9999
        ENDIF
        
        DO pnd = 1,ISubgroupSize
           IDisplacements(pnd) = (IPixelTotal*(pnd-1)/ISubgroupSize)
           ICount(pnd) = (((IPixelTotal*(pnd)/ISubgroupSize) - (IPixelTotal*(pnd-1)/ISubgroupSize)))*IReflectOut*IThicknessCount
           
        END DO
        
        DO ind = 1,ISubgroupSize
           IDisplacements(ind) = (IDisplacements(ind))*IReflectOut*IThicknessCount
        END DO
        
        IF(IImageFLAG.LE.2) THEN
           CALL MPI_GATHERV(RIndividualReflections,ICount,&
                MPI_DOUBLE_PRECISION,RIndividualReflectionsRoot,&
                ICount,IDisplacements,MPI_DOUBLE_PRECISION,0,&
                Inewcomm,IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " In MPI_GATHERV"
              GOTO 9999
           ENDIF
        ELSE     
           CALL MPI_GATHERV(CAmplitudeandPhase,ICount,&
                MPI_DOUBLE_COMPLEX,CAmplitudeandPhaseRoot,&
                ICount,IDisplacements,MPI_DOUBLE_COMPLEX,0, &
                Inewcomm,IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " In MPI_GATHERV"
              GOTO 9999
           ENDIF
        END IF
        
        
        IF(IImageFLAG.LE.2) THEN
           DEALLOCATE( &
                RIndividualReflections,STAT=IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " Deallocating RIndividualReflections"
              GOTO 9999
           ENDIF
        END IF
        IF (my_newrank.EQ.0) THEN !This HAS to be WRONG
           RCrossCorrelationOld = ZERO
           
           ALLOCATE( &
                RReflectionImagesForPhaseCorrelation(2*IPixelCount,2*IPixelCount,&
                IReflectOut,IThicknessCount),&
                STAT=IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " in ALLOCATE() of DYNAMIC variables Root Reflections for phase correlation"
              GOTO 9999
           ENDIF
           
           RReflectionImagesForPhaseCorrelation = ZERO
           DO ind = 1,IReflectOut
              DO knd = 1,IThicknessCount
                 DO jnd = 1,IPixelTotal
                    gnd = IPixelLocations(jnd,1)
                    hnd = IPixelLocations(jnd,2)
                    
                    RReflectionImagesForPhaseCorrelation(gnd,hnd,ind,knd) = RIndividualReflectionsRoot(ind,knd,jnd)
                 END DO
              END DO
           END DO
           
           DO ind = 1,IThicknessCount
              CALL PhaseCorrelate(RReflectionImagesForPhaseCorrelation(:,:,IOutputReflections(1),ind),&
                   RImageExpi(:,:,1),IErr,2*IPixelCount,2*IPixelCount)
              
              IF(RCrossCorrelation.GT.RCrossCorrelationOld) THEN
                 RCrossCorrelationOld = RCrossCorrelation
                 RThickness = RInitialThickness + (ind-1)*RDeltaThickness 
                 IThicknessCountFinal = ind
                 RBestCorrelationImage = RReflectionImagesForPhaseCorrelation(:,:,IOutputReflections(:),ind)
              END IF
              
              
           END DO
           DEALLOCATE(&
                RReflectionImagesForPhaseCorrelation,STAT=IErr)       
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error in Deallocation of RReflectionImagesForPhaseCorrelation"
              GOTO 9999
           END IF
           
           RFluxIterationCorrelations(fnd,1) = RThickness
           RFluxIterationCorrelations(fnd,2) = RCrossCorrelationOld/(2*IPixelCount)**2
           
!!$     PRINT*,"Thickness = ",RThickness," Angstoms with a correlation of ",RCrossCorrelationOld/(2*IPixelCount)**2
        END IF
        
        
        DEALLOCATE(RIndividualReflectionsRoot,&
             STAT=IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                " in DEALLOCATE() of DYNAMIC variables Root Reflections"
           GOTO 9999
        ENDIF
        DEALLOCATE(&
             RDevPara,STAT=IErr)       
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error in Deallocation of RDevPara etc"
           GOTO 9999
        ENDIF
        DEALLOCATE(&
             IWeakBeamList,STAT=IErr)       
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error in Deallocation of IWeakBeamList etc"
           GOTO 9999
        ENDIF
        DEALLOCATE(&
             IStrongBeamList,STAT=IErr)       
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error in Deallocation of IStrongBeamList etc"
           GOTO 9999
        ENDIF
        DEALLOCATE(&
             CFullWaveFunctions,STAT=IErr)       
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error in Deallocation of CFullWaveFunctions etc"
           GOTO 9999
        ENDIF
        DEALLOCATE(&
             RFullWaveIntensity,STAT=IErr)       
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error in Deallocation of RFullWaveIntensity etc"
           GOTO 9999
        ENDIF
        DEALLOCATE(&
             IDisplacements,STAT=IErr)       
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error in Deallocation of IDisplacements etc"
           GOTO 9999
        ENDIF
        DEALLOCATE(&
             ICount,STAT=IErr)       
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error in Deallocation of ICount etc"
           GOTO 9999
        END IF

        PRINT*,"Exiting 1"
        
     END DO
     
     CALL MPI_ALLREDUCE(RFluxIterationCorrelations,RFluxIterationCorrelationsRoot,&
          2*IFluxIterationSteps,MPI_DOUBLE_PRECISION,MPI_SUM,&
          MPI_COMM_WORLD,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
             " In MPI_REDUCE"
        GOTO 9999
     ENDIF
     
     ICorrelationMaximum = MAXLOC(RFluxIterationCorrelationsRoot(:,2),1)
     
     RThickness = RFluxIterationCorrelationsRoot(ICorrelationMaximum,1)
     
     SELECT CASE (IRefineModeFLAG)
     CASE(0)
        
        DO ind = 1,IMaxiter
           
           IF(ind.EQ.1) THEN
              
              IFluxIterationIndices(ind) = FLOOR(REAL((ICorrelationMaximum-1))/REAL(INoofDWFs**(SIZE(IElementlist)-ind)))+1
              IRemainder = MOD(ICorrelationMaximum-1,INoofDWFs**(SIZE(IElementlist)-ind))
           ELSE
              IF(ind.EQ.IMaxiter) THEN
                 IFluxIterationIndices(ind) = IRemainder+1             
              ELSE
                 IFluxIterationIndices(ind) = FLOOR(REAL(IRemainder)/REAL(INoofDWFs**(SIZE(IElementlist)-ind)))+1
                 IRemainder = MOD(IRemainder,INoofDWFs**(SIZE(IElementlist)-ind))
              END IF
           END IF
        END DO
        
        DO ind = 1,IMaxIter
           RFinalDWFConfig(ind) = &
                ((IFluxIterationIndices(ind)-1)*&
                RDeltaDebyeWallerFactor)+&
                RInitialDebyeWallerFactor
        END DO
        
        IF(my_rank.EQ.0) THEN
           PRINT*,"Maximum Correlation was ",&
                RFluxIterationCorrelationsRoot(ICorrelationMaximum,2),&
                "from a thickness of ",&
                RFluxIterationCorrelationsRoot(ICorrelationMaximum,1),&
                "and a debyewaller Factor configuration of"
           DO ind = 1,IMaxIter
              PRINT*,RFinalDWFConfig(ind)
           END DO
        END IF
        
     CASE(1)
     CASE(2)
        IF(my_rank.EQ.0) THEN
           PRINT*,"Maximum Correlation was ",&
                RFluxIterationCorrelationsRoot(ICorrelationMaximum,2),&
                "from a thickness of ",&
                RFluxIterationCorrelationsRoot(ICorrelationMaximum,1)
        END IF
        
     END SELECT
     
     DEALLOCATE(&
          IFluxIterationIndices,STAT=IErr)       
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error in Deallocation of IFluxIterationIndices etc"
        GOTO 9999
     END IF
     DEALLOCATE(&
          RFluxIterationCorrelations,STAT=IErr)       
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error in Deallocation of RFluxIterationCorrelations etc"
        GOTO 9999
     END IF
     DEALLOCATE(&
          RFluxIterationCorrelationsRoot,STAT=IErr)       
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error in Deallocation of RFluxIterationCorrelationsRoot etc"
        GOTO 9999
     END IF
     
     SELECT CASE (IRefineModeFLAG)
     CASE(0)
        IF(RDeltaDebyeWallerFactor.LE.RIterationTolerance) THEN
           EXIT
        END IF
     CASE(1)
     CASE(2)
        IF(RDeltaThickness.LE.RIterationTolerance) THEN
           EXIT
        END IF
        
     CASE DEFAULT 
     END SELECT

     PRINT*,"Exiting 2"
     
  END DO

  ISubgroupContainingAnswer = &
       FLOOR(REAL(ICorrelationMaximum,RKIND)/REAL(IFluxStepsPerSubgroup,RKIND))
  ISubgroupRootProcessContainingAnswer = (ISubgroupContainingAnswer-1_IKIND)*ISubgroupSize

  PRINT*,ISubgroupRootProcessContainingAnswer,ISubgroupContainingAnswer
  
  IF(my_rank.EQ.ISubgroupRootProcessContainingAnswer) THEN
    PRINT*,"Hi, Im rank",my_rank,"im about to write out the result"
     ALLOCATE(&
          RFinalMontageImage(MAXVAL(IImageSizeXY),&
          MAXVAL(IImageSizeXY),1),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Montage"
        GOTO 9999
     ENDIF

     PRINT*,"Memory Allocated"
     RFinalMontageImage = ZERO
     DO ind = 1,2*IPixelCount
        DO jnd = 1,2*IPixelCount
           CALL MakeMontagePixel(ind,jnd,1,&
                RFinalMontageImage,&
                RBestCorrelationImage(ind,jnd,:),IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                   " in MakeMontagePixel"
              GOTO 9999
           ENDIF
        END DO
        !PRINT*,RFinalMontageImage(ind,:,1)
     END DO
     PRINT*,"Have Created Montage IMage"

     
     DEALLOCATE(&
          RBestCorrelationImage,STAT=IErr)       
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error in Deallocation of RBestCorrelationImage etc"
        GOTO 9999
     END IF
     
     WRITE(surname,"(A2,A1,I5.5,A2,I5.5)") &
          "M-","T",NINT(RThickness),"-P",MAXVAL(IImageSizeXY)  
!!$     PRINT*,surname
     
     PRINT*,"Openning File"

     CALL MPI_FILE_OPEN( MPI_COMM_SELF, TRIM(surname), &
          MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, IChOutWI_MPI, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
             " in MPI_FILE_OPEN"
        GOTO 9999
     ENDIF
     PRINT*,"File Open"
     
     DO ind = 1,MAXVAL(IImageSizeXY)
        CALL MPI_FILE_WRITE(IChOutWI_MPI,RFinalMontageImage(ind,:,1),MAXVAL(IImageSizeXY),&
             MPI_DOUBLE_PRECISION,my_status,IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
                " in MPI_FILE_WRITE"
           GOTO 9999
        ENDIF
     END DO
     PRINT*,"FILE WRITTEN"
     CALL MPI_FILE_CLOSE(IChOutWI_MPI,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"RefineMain(", my_rank, ") error ", IErr, &
             " in MPI_FILE_CLOSE"
        GOTO 9999
     ENDIF
     PRINT*,"FILE CLOSED"
  END IF
  
  IF(my_rank.EQ.0) THEN
     
     IF(IImageFLAG.GE.3) THEN
        DEALLOCATE(&
             CAmplitudeandPhaseRoot,STAT=IErr) 
        
        IF( IErr.NE.0 ) THEN
           PRINT*,"RefineMain(", my_rank, ") error in Deallocation of CAmplitudeandPhase"
           GOTO 9999
        ENDIF
     END IF
  END IF
  
  DEALLOCATE( &
       RgVecMatT,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in Deallocation of RgVecMatT etc"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       Rhklpositions,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in Deallocation of Rhklpositions etc"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       RMask,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in Deallocation of RMask etc"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       CUgMat,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in Deallocation of CUgmat etc"
     GOTO 9999
  ENDIF

  DEALLOCATE( &
       IPixelLocations,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error in Deallocation of IPixelLocations etc"
     GOTO 9999
  ENDIF
  
  !--------------------------------------------------------------------
  ! finish off
  !--------------------------------------------------------------------
    
  CALL cpu_time(CurrentTime)
  Duration=(CurrentTime-StartTime)
  IHours = FLOOR(Duration/3600.0D0)
  IMinutes = FLOOR(MOD(Duration,3600.0D0)/60.0D0)
  ISeconds = MOD(Duration,3600.0D0)-IMinutes*60.0D0

  PRINT*, "FelixRefine(", my_rank, ") ", RStr, ", used time=", IHours, "hrs ",IMinutes,"mins ",ISeconds,"Seconds "

  !--------------------------------------------------------------------
  ! Shut down MPI
  !--------------------------------------------------------------------
9999 &
  CALL MPI_Finalize(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RefineMain(", my_rank, ") error ", IErr, " in MPI_Finalize()"
     STOP
  ENDIF

  ! clean shutdown
  STOP
  

END PROGRAM FelixRefine
