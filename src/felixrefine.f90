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

PROGRAM Felixrefine
 
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

  INTEGER(IKIND) :: &
       IHours,IMinutes,ISeconds,IErr,IMilliSeconds,IIterationFLAG,&
       ind,IIterationCount,ISpaceGrp,ICount,jnd
  REAL(RKIND) :: &
       StartTime, CurrentTime, Duration, TotalDurationEstimate,&
       RFigureOfMerit,SimplexFunction  
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RSimplexVolume
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       RSimplexFoM,RIndependentVariableValues
  REAL(RKIND) :: &
       RBCASTREAL
  INTEGER(IKIND),DIMENSION(:),ALLOCATABLE :: &
       IVectors
  REAL(RKIND) :: &
       RStandardDeviation,RMean

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
     PRINT*,"Felixrefine(", my_rank, ") error in MPI_Init()"
     GOTO 9999
  ENDIF

  ! Get the rank of the current process
  CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in MPI_Comm_rank()"
     GOTO 9999
  ENDIF

  ! Get the size of the current communicator
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in MPI_Comm_size()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! protocal feature startup
  !--------------------------------------------------------------------
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"--------------------------------------------------------------"
     PRINT*,"Felixrefine: ", RStr
     PRINT*,"          ", DStr
     PRINT*,"          ", AStr
     PRINT*,"          on rank= ", my_rank, " of ", p, " in total."
     PRINT*,"--------------------------------------------------------------"
  END IF

  !--------------------------------------------------------------------
  ! timing startup
  !--------------------------------------------------------------------

  CALL cpu_time(StartTime)

  !--------------------------------------------------------------------
  ! INPUT section 
  !--------------------------------------------------------------------
  
  ISoftwareMode =2 ! felixrefinemode

  !Read from input files
  CALL ReadInput (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in ReadInput()"
     GOTO 9999
  ENDIF

  IF(IRefineModeSelectionArray(2).EQ.1) THEN 
     
     CALL ConvertSpaceGroupToNumber(ISpaceGrp,IErr)
     
     ALLOCATE(&
          IVectors(SIZE(SWyckoffSymbols)),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixrefine (", my_rank, ") error in Allocation() of IVectors"
        GOTO 9999
     ENDIF
     
     DO ind = 1,SIZE(SWyckoffSymbols)
!!$     SWyckoffSymbol = SWyckoffSymbols(ind)
        CALL CountAllowedMovements(ISpaceGrp,SWyckoffSymbols(ind),IVectors(ind),IErr)
     END DO
     
     IAllowedVectors = SUM(IVectors)
     
     ALLOCATE(&
          IAllowedVectorIDs(IAllowedVectors),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixrefine (", my_rank, ") error in Allocation() of IVectors"
        GOTO 9999
     ENDIF
     
     ICount = 0
     
     DO ind = 1,SIZE(SWyckoffSymbols)
        DO jnd = 1,IVectors(ind)
           ICount = ICount + 1
           IAllowedVectorIDs(ICount) = IAtomicSitesToRefine(ind)
        END DO
     END DO
     
     ALLOCATE(&
          RAllowedVectors(IAllowedVectors,THREEDIM),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixrefine (", my_rank, ") error in Allocation() of IVectors"
        GOTO 9999
     ENDIF
     
     ALLOCATE(&
          RAllowedVectorMagnitudes(IAllowedVectors),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixrefine (", my_rank, ") error in Allocation() of IVectors"
        GOTO 9999
     ENDIF
     
     RAllowedVectorMagnitudes = ZERO
     
     DO ind = 1,SIZE(SWyckoffSymbols)
!!$     SWyckoffSymbol = SWyckoffSymbols(ind)
        CALL DetermineAllowedMovements(ISpaceGrp,SWyckoffSymbols(ind),&
             RAllowedVectors(SUM(IVectors(:(ind-1)))+1:SUM(IVectors(:(ind))),:),&
             IVectors(ind),IErr)
     END DO

  END IF
     
!!$Calculate Number of Independent Refinement Variables
  
  IIndependentVariables = &
       IRefineModeSelectionArray(1)*INoofUgs*2+&
       IRefineModeSelectionArray(2)*IAllowedVectors+&
       IRefineModeSelectionArray(3)*SIZE(IAtomicSitesToRefine)+&
       IRefineModeSelectionArray(4)*SIZE(IAtomicSitesToRefine)+&
       IRefineModeSelectionArray(5)*SIZE(IAtomicSitesToRefine)*6+&
       IRefineModeSelectionArray(6)*3+&
       IRefineModeSelectionArray(7)*3+&
       IRefineModeSelectionArray(8)
  
  ALLOCATE( &
       RImageExpi(2*IPixelCount,2*IPixelCount, &
       IReflectOut),&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Allocation()"
     GOTO 9999
  ENDIF

  CALL ReadExperimentalImages(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in ReadExperimentalImages()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! Save Atomic Coordinates  
  !--------------------------------------------------------------------

  ALLOCATE(&
       RInitialAtomSiteFracCoordVec(&
       SIZE(RAtomSiteFracCoordVec,DIM=1),&
       SIZE(RAtomSiteFracCoordVec,DIM=2)),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in ALLOCATE()RInitialAtomSiteFracCoordVec "
     GOTO 9999
  ENDIF  

  RInitialAtomSiteFracCoordVec = RAtomSiteFracCoordVec

  !--------------------------------------------------------------------
  ! Setup Simplex Variables
  !--------------------------------------------------------------------

  CALL AssignIterativeIDs(IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in AssignIterativeIDs()"
     GOTO 9999
  ENDIF
  
  ALLOCATE(&
       RIndependentVariableValues(IIndependentVariables),&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in allocation()"
     GOTO 9999
  ENDIF

  IF(my_rank.EQ.0) THEN
     PRINT*,"IErr = ",IErr
  END IF

  CALL RefinementVariableSetup(RIndependentVariableValues,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in RefinementVariableSetup()"
     GOTO 9999
  ENDIF

  
  !--------------------------------------------------------------------
  ! Initialise Simplex
  !--------------------------------------------------------------------

  ALLOCATE( &
       RSimplexVolume(IIndependentVariables+1,IIndependentVariables),&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Allocation()"
     GOTO 9999
  ENDIF

  ALLOCATE( &
       RSimplexFoM(IIndependentVariables),&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Allocation()"
     GOTO 9999
  ENDIF
  
  IIterationCount = 0

  CALL SimplexInitialisation(RSimplexVolume,RSimplexFoM,RIndependentVariableValues,IIterationCount,RStandardDeviation,RMean,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in SimplexInitialisation()"
     GOTO 9999
  ENDIF
     
  !--------------------------------------------------------------------
  ! Apply Simplex Method
  !--------------------------------------------------------------------

  CALL NDimensionalDownhillSimplex(RSimplexVolume,RSimplexFoM,&
       IIndependentVariables+1,&
       IIndependentVariables,IIndependentVariables,&
       0.0001d0,IIterationCount,RStandardDeviation,RMean,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in NDimensionalDownhillSimplex()"
     GOTO 9999
  ENDIF

  CALL MPI_BARRIER(MPI_COMM_WORLD,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error ", IErr, " in MPI_BARRIER()"
     STOP
  ENDIF

  PRINT*,"Im rank",my_rank

  STOP

  !--------------------------------------------------------------------
  ! Deallocate Memory
  !--------------------------------------------------------------------

  DEALLOCATE( &
       RImageExpi,&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Deallocation()"
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
  IMilliSeconds = INT((Duration-(IHours*3600+IMinutes*60+ISeconds))*100,IKIND)

  PRINT*, "Felixrefine(", my_rank, ") ", RStr, ", used time=", IHours, "hrs ",IMinutes,"mins ",ISeconds,"Seconds ",&
       IMilliSeconds,"Milliseconds"
  !--------------------------------------------------------------------
  ! Shut down MPI
  !--------------------------------------------------------------------

  CALL MPI_BARRIER(MPI_COMM_WORLD,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error ", IErr, " in MPI_BARRIER()"
     STOP
  ENDIF

9999 &
  CALL MPI_Finalize(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error ", IErr, " in MPI_Finalize()"
     STOP
  ENDIF
  
  ! clean shutdown
  STOP
  
END PROGRAM Felixrefine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE AssignIterativeIDs(IErr)

USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: &
       ind,jnd,IErr,ICalls
  INTEGER(IKIND),DIMENSION(IRefinementVariableTypes) :: &
       INoofelementsforeachrefinementtype

  CALL DetermineNumberofRefinementVariablesPerType(INoofelementsforeachrefinementtype,IErr)

  ALLOCATE(&
       IIterativeVariableUniqueIDs(IIndependentVariables,5),&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Allocation() of IIterativeVariableUniqueIDs"
     RETURN
  ENDIF

  IIterativeVariableUniqueIDs = 0
  ICalls = 0

  DO ind = 1,IRefinementVariableTypes !Loop over all possible iterative variables
     IF(IRefineModeSelectionArray(ind).EQ.1) THEN
        DO jnd = 1,INoofelementsforeachrefinementtype(ind)
           ICalls = ICalls + 1
           IIterativeVariableUniqueIDs(ICalls,1) = ICalls
           CALL AssignArrayLocationsToIterationVariables(ind,jnd,IIterativeVariableUniqueIDs,IErr)
        END DO
     END IF
  END DO
  
END SUBROUTINE AssignIterativeIDs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE AssignArrayLocationsToIterationVariables(IIterativeVariableType,IVariableNo,IArrayToFill,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IIterativeVariableType,IVariableNo,IErr,IArrayIndex,&
       IAnisotropicDebyeWallerFactorElementNo
  INTEGER(IKIND),DIMENSION(IIndependentVariables,5),INTENT(OUT) :: &
       IArrayToFill  
  INTEGER(IKIND),DIMENSION(IRefinementVariableTypes) :: &
       INoofelementsforeachrefinementtype

!!$  Calculate How Many of Each Variable Type There are

  CALL DetermineNumberofRefinementVariablesPerType(INoofelementsforeachrefinementtype,IErr)
  
!!$  Where am I in the Array Right Now?

  IArrayIndex = SUM(INoofelementsforeachrefinementtype(:(IIterativeVariableType-1)))+IVariableNo

  SELECT CASE(IIterativeVariableType)

  CASE(1) ! Ugs

     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = &
          NINT(REAL(INoofUgs,RKIND)*(REAL(IVariableNo/REAL(INoofUgs,RKIND),RKIND)-&
          CEILING(REAL(IVariableNo/REAL(INoofUgs,RKIND),RKIND)))+REAL(INoofUgs,RKIND))

  CASE(2) ! Coordinates (x,y,z)

     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
!!$     IArrayToFill(IArrayIndex,3) = IAtomicSitesToRefine(INT(CEILING(REAL(IVariableNo/3.0D0,RKIND))))
!!$     IArrayToFill(IArrayIndex,4) = &
!!$          NINT(3.D0*(REAL(IVariableNo/3.0D0,RKIND)-CEILING(REAL(IVariableNo/3.0D0,RKIND)))+3.0D0)
!!$     IArrayToFill(IArrayIndex,5) = 0

     IArrayToFill(IArrayIndex,3) = IVariableNo

  CASE(3) ! Atomic Site Occupancies

     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = IAtomicSitesToRefine(IVariableNo)

  CASE(4) ! Isotropic Debye Waller Factors 

     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = IAtomicSitesToRefine(IVariableNo)

  CASE(5) ! Anisotropic Debye Waller Factors (a11-a33)

     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = IAtomicSitesToRefine(INT(CEILING(REAL(IVariableNo/6.0D0,RKIND))))
     IAnisotropicDebyeWallerFactorElementNo = &
          NINT(6.D0*(REAL(IVariableNo/6.0D0,RKIND)-CEILING(REAL(IVariableNo/6.0D0,RKIND)))+6.0D0)

     SELECT CASE(IAnisotropicDebyeWallerFactorElementNo)

        CASE(1)
           IArrayToFill(IArrayIndex,4:5) = [1,1]
        CASE(2)
           IArrayToFill(IArrayIndex,4:5) = [2,1]
        CASE(3)
           IArrayToFill(IArrayIndex,4:5) = [2,2]
        CASE(4)
           IArrayToFill(IArrayIndex,4:5) = [3,1]
        CASE(5)
           IArrayToFill(IArrayIndex,4:5) = [3,2]
        CASE(6)
           IArrayToFill(IArrayIndex,4:5) = [3,3]

        END SELECT

  CASE(6) ! Lattice Parameters

     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = &
          NINT(3.D0*(REAL(IVariableNo/3.0D0,RKIND)-CEILING(REAL(IVariableNo/3.0D0,RKIND)))+3.0D0)
     
  CASE(7) ! Lattice Angles

     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = &
          NINT(3.D0*(REAL(IVariableNo/3.0D0,RKIND)-CEILING(REAL(IVariableNo/3.0D0,RKIND)))+3.0D0)

  CASE(8)
     
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType

  CASE(9)
     
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     
  END SELECT
  
END SUBROUTINE AssignArrayLocationsToIterationVariables

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE RefinementVariableSetup(RIndependentVariableValues,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       IErr,ind,IVariableType
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(OUT) :: &
       RIndependentVariableValues
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RefinementVariableSetup(",my_rank,")"
  END IF
  
!!$  Fill the Independent Value array with values
  
  DO ind = 1,IIndependentVariables
     IVariableType = IIterativeVariableUniqueIDs(ind,2)
     SELECT CASE (IVariableType)
     CASE(1)
     CASE(2)
        RIndependentVariableValues(ind) = &
             RAllowedVectorMagnitudes(IIterativeVariableUniqueIDs(ind,3))
!!$        RIndependentVariableValues(ind) = &
!!$             RAtomSiteFracCoordVec(&
!!$             IIterativeVariableUniqueIDs(ind,3),&
!!$             IIterativeVariableUniqueIDs(ind,4))
     CASE(3)
        RIndependentVariableValues(ind) = &
             RAtomicSitePartialOccupancy(IIterativeVariableUniqueIDs(ind,3))
     CASE(4)
        RIndependentVariableValues(ind) = &
             RIsotropicDebyeWallerFactors(IIterativeVariableUniqueIDs(ind,3))
     CASE(5)
        RIndependentVariableValues(ind) = &
             RAnisotropicDebyeWallerFactorTensor(&
             IIterativeVariableUniqueIDs(ind,3),&
             IIterativeVariableUniqueIDs(ind,4),&
             IIterativeVariableUniqueIDs(ind,5))
     CASE(6)
        SELECT CASE(IIterativeVariableUniqueIDs(ind,3))
        CASE(1)
           RIndependentVariableValues(ind) = RLengthX
        CASE(2)
           RIndependentVariableValues(ind) = RLengthY
        CASE(3)
           RIndependentVariableValues(ind) = RLengthZ
        END SELECT
     CASE(7)
        SELECT CASE(IIterativeVariableUniqueIDs(ind,3))
        CASE(1)
           RIndependentVariableValues(ind) = RAlpha
        CASE(2)
           RIndependentVariableValues(ind) = RBeta
        CASE(3)
           RIndependentVariableValues(ind) = RGamma
        END SELECT
     CASE(8)
        RIndependentVariableValues(ind) = &
             RConvergenceAngle
     CASE(9)
        RIndependentVariableValues(ind) = &
             RAbsorptionPercentage
     END SELECT
  END DO

END SUBROUTINE RefinementVariableSetup
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE StructureFactorRefinementSetup(RIndependentVariableValues,IIterationCount,IErr)
  
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
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(OUT) :: &
       RIndependentVariableValues
  INTEGER(IKIND),INTENT(IN) :: &
       IIterationCount

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"StructureFactorRefinementSetup(",my_rank,")"
  END IF

  IF(IRefineModeSelectionArray(1).EQ.1.AND.IIterationCount.EQ.1) THEN
     DO ind = 1,INoofUgs
        RIndependentVariableValues((ind-1)*2+1) = &
             REAL(CSymmetryStrengthKey(ind),RKIND)
        RIndependentVariableValues((ind-1)*2+2) = &
             AIMAG(CSymmetryStrengthKey(ind))
     END DO
  END IF

END SUBROUTINE StructureFactorRefinementSetup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE RankSymmetryRelatedStructureFactor(IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 
  
  INTEGER(IKIND) :: &
       IErr,ind
  INTEGER(IKIND),DIMENSION(2) :: &
       ILoc

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RankSymmetryRelatedStructureFactor(",my_rank,")"
  END IF
  
  ALLOCATE( &  
       ISymmetryRelations(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RankSymmetryRelatedStructureFactor(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables ISymmetryRelations"
     RETURN
  ENDIF
  
  CALL SymmetryRelatedStructureFactorDetermination (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RankSymmetryRelatedStructureFactor(", my_rank, ") error ", IErr, &
          " in SymmetryRelatedStructureFactorDetermination"
     RETURN
  ENDIF
  
  DO ind = 1,(SIZE(ISymmetryStrengthKey,DIM=1))
     ILoc = MINLOC(ABS(ISymmetryRelations-ind))
     ISymmetryStrengthKey(ind,1) = ind
     ISymmetryStrengthKey(ind,2) = ind
     CSymmetryStrengthKey(ind) = CUgMat(ILoc(1),ILoc(2))
  END DO
  
  CALL ReSortUgs(ISymmetryStrengthKey,CSymmetryStrengthKey,SIZE(CSymmetryStrengthKey,DIM=1))

END SUBROUTINE RankSymmetryRelatedStructureFactor

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SimplexInitialisation(RSimplexVolume,RSimplexFoM,RIndependentVariableValues,&
     IIterationCount,RStandardDeviation,RMean,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind,jnd
  REAL(RKIND),DIMENSION(IIndependentVariables+1,IIndependentVariables),INTENT(OUT) :: &
       RSimplexVolume
  REAL(RKIND),DIMENSION(IIndependentVariables+1),INTENT(OUT) :: &
       RSimplexFoM
  REAL(RKIND) :: &
       SimplexFunction,RSimplexDummy
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(INOUT) :: &
       RIndependentVariableValues
  INTEGER(IKIND),INTENT(INOUT) :: &
       IIterationCount
  REAL(RKIND),INTENT(OUT) :: &
       RStandardDeviation,RMean
  REAL(RKIND) :: &
       RStandardError,RStandardTolerance

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"SimplexInitialisation(",my_rank,")"
  END IF
      
  CALL PerformDummySimulationToSetupSimplexValues(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation(", my_rank, ") error in PerformDummySimulationToSetupSimplexValues()"
     RETURN
  ENDIF
  
  CALL InitialiseWeightingCoefficients(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation(", my_rank, ") error in InitialiseWeightingCoefficients()"
     RETURN
  ENDIF

  CALL RefinementVariableSetup(RIndependentVariableValues,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation(", my_rank, ") error in RefinementVariableSetup()"
     RETURN
  ENDIF

  
  IF(IRefineModeSelectionArray(1).EQ.1) THEN
     
     CALL RankSymmetryRelatedStructureFactor(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexInitialisation(", my_rank, ") error in RankSymmetryRelatedStructureFactor()"
        RETURN
     ENDIF
     
     CALL StructureFactorRefinementSetup(RIndependentVariableValues,IIterationCount,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexInitialisation(", my_rank, ") error in StructureFactorRefinementSetup()"
        RETURN
     ENDIF
     
  END IF
  
  DEALLOCATE(&
       CUgmat,&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation (", my_rank, ") error in Deallocation()"
     RETURN
  ENDIF

!!$ RandomSequence

  IF(IContinueFLAG.EQ.0) THEN
     IF(my_rank.EQ.0) THEN
        CALL CreateRandomisedSimplex(RSimplexVolume,&
             RIndependentVariableValues,IErr)

        CALL MPI_BCAST(RSimplexVolume,(IIndependentVariables+1)*(IIndependentVariables),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
     ELSE
        CALL MPI_BCAST(RSimplexVolume,(IIndependentVariables+1)*(IIndependentVariables),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)

     END IF

     IPreviousPrintedIteration = -IPrint ! Ensures print out on first iteration

     DO ind = 1,(IIndependentVariables+1)
        
        IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"---------------------------------------------------------"
           PRINT*,"-------- Simplex",ind,"of",IIndependentVariables+1
           PRINT*,"---------------------------------------------------------"
        END IF
                

        RSimplexDummy = SimplexFunction(RSimplexVolume(ind,:),1,0,IErr)
        
        RStandardTolerance = RStandardError(RStandardDeviation,RMean,RSimplexDummy,IErr)
        
        RSimplexFoM(ind) =  RSimplexDummy
        
        IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"---------------------------------------------------------"
           PRINT*,"-------- Figure of Merit" ,RSimplexFoM(ind)        
           PRINT*,"---------------------------------------------------------"
        END IF
     END DO
     
  ELSE
     
     CALL RecoverSavedSimplex(RSimplexVolume,RSimplexFoM,RStandardDeviation,RMean,IIterationCount,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexInitialisation (", my_rank, ") error in RecoverSavedSimplex()"
        RETURN
     ENDIF
     
  END IF
       
END SUBROUTINE SimplexInitialisation

SUBROUTINE CreateRandomisedSimplex(RSimplexVolume,RIndependentVariableValues,IErr)

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
       RRandomSigns,RRandomNumbers
  REAL(RKIND),DIMENSION(IIndependentVariables+1,IIndependentVariables),INTENT(OUT) :: &
       RSimplexVolume
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(INOUT) :: &
       RIndependentVariableValues
  
  IF(IRefineModeSelectionArray(2).EQ.1) THEN
     DO ind = 1,(IIndependentVariables+1)
        ALLOCATE(&
             RRandomSigns(IAllowedVectors),&
             RRandomNumbers(IAllowedVectors),&
             STAT=IErr)
        
        
!!$           Randomise Atomic Displacements
        
        CALL RandomSequence(RRandomNumbers,IAllowedVectors,ind,IErr)
        CALL RandomSequence(RRandomSigns,IAllowedVectors,2*ind,IErr)
        WHERE (RRandomSigns.LT.HALF)
           RRandomSigns = ONE
        ELSEWHERE
           RRandomSigns = -ONE
        END WHERE
        
        RSimplexVolume(ind,:IAllowedVectors) = &
             RRandomNumbers*RRandomSigns*RSimplexLengthScale
        
        DEALLOCATE(&
             RRandomSigns,&
             RRandomNumbers)
        
        ALLOCATE(&
             RRandomSigns(IIndependentVariables-IAllowedVectors),&
             RRandomNumbers(IIndependentVariables-IAllowedVectors),&
             STAT=IErr)
        
!!$           Randomise Everything else
        
        CALL RandomSequence(RRandomNumbers,&
             IIndependentVariables-IAllowedVectors,ind,IErr)
        CALL RandomSequence(RRandomSigns,&
             IIndependentVariables-IAllowedVectors,2*ind,IErr)
        WHERE (RRandomSigns.LT.HALF)
           RRandomSigns = ONE
        ELSEWHERE
           RRandomSigns = -ONE
        END WHERE
        
        RSimplexVolume(ind,(IAllowedVectors+1):) = &
             RIndependentVariableValues((IAllowedVectors+1):)*&
             (1+(RRandomNumbers*RRandomSigns*RSimplexLengthScale))
        
        DEALLOCATE(&
             RRandomSigns,&
             RRandomNumbers)
        
     END DO
     
  ELSE
     ALLOCATE(&
          RRandomSigns(IIndependentVariables),&
          RRandomNumbers(IIndependentVariables),&
          STAT=IErr)
     
     DO ind = 1,(IIndependentVariables+1)
        
        CALL RandomSequence(RRandomNumbers,IIndependentVariables,ind,IErr)
        CALL RandomSequence(RRandomSigns,IIndependentVariables,2*ind,IErr)
        WHERE (RRandomSigns.LT.HALF)
           RRandomSigns = ONE
        ELSEWHERE
           RRandomSigns = -ONE
        END WHERE
        
        
        RSimplexVolume(ind,:) = &
             RIndependentVariableValues(:)*&
             (1+(RRandomNumbers*RRandomSigns*RSimplexLengthScale))
     END DO
     
     DEALLOCATE(&
          RRandomSigns,&
          RRandomNumbers)
     
  END IF
  

END SUBROUTINE CreateRandomisedSimplex

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE PerformDummySimulationToSetupSimplexValues(IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Calls setup subroutines to initialise Structure factors for 
!!$  % use in initialisation of simplex
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(",my_rank,")"
  END IF
  
  !-------------------------------------------------------------------- 
  !Setup Experimental Variables
  !--------------------------------------------------------------------

  CALL ExperimentalSetup (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error in ExperimentalSetup()"
     RETURN
  ENDIF
    
  !--------------------------------------------------------------------
  ! Setup Image
  !--------------------------------------------------------------------

  CALL ImageSetup( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error in ImageSetup()"
     RETURN
  ENDIF

 
  !--------------------------------------------------------------------
  ! MAIN section
  !--------------------------------------------------------------------
   
  IF(IAbsorbFLAG.NE.0) THEN ! Calculate Non-absorbative UGs
     
     IAbsorbFLAG = 0
     
     CALL StructureFactorSetup(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error in StructureFactorSetup()"
        RETURN
     ENDIF
     
     IAbsorbFLAG = 1
     
  END IF
  
  DEALLOCATE( &
       MNP,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation MNP"
     RETURN
  ENDIF
      
  DEALLOCATE( &
        SMNP, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation SMNP"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RFullAtomicFracCoordVec, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation RFullAtomicFracCoordVec"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       SFullAtomicNameVec,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation SFullAtomicNameVec"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IFullAnisotropicDWFTensor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation IFullAnisotropicDWFTensor"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IFullAtomNumber,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation IFullAtomNumber"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RFullIsotropicDebyeWallerFactor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation RFullIsotropicDebyeWallerFactor"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RFullPartialOccupancy,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation RFullPartialOccupancy"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RDWF,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation RDWF"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       ROcc,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation ROcc"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IAtoms,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation IAtoms"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IAnisoDWFT,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation IAnisoDWFT"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       Rhkl,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation Rhkl"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RgVecMatT,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMatT"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RgVecMag,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMag"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RGn,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation RGn" 
     RETURN
  ENDIF
       
  DEALLOCATE( &
       Rhklpositions,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation Rhklpositions"
     RETURN
  ENDIF  
       
  DEALLOCATE( &
       RMask,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation RMask"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IPixelLocations,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation IPixelLocations"
     RETURN
  ENDIF

  IDiffractionFLAG = 0

END SUBROUTINE PerformDummySimulationToSetupSimplexValues

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE InitialiseAtomicVectorMagnitudes(IVariableID,RCorrectedMovement,IErr)
  
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Creates pseudo random movements of atoms using allowed vectors
!!$  % to initialise the simplex, proposed movements which exit the unit
!!$  $ cell are corrected to bring the atom back in on the opposite side
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind,IVariableID
  REAL(RKIND) :: &
       RNegativeMovement,RPositiveMovement,RCorrectedMovement,RANDOMNUMBER
  RNegativeMovement = RSimplexLengthScale*(-1.0_RKIND)
  RPositiveMovement = RSimplexLengthScale

  IF(RANDOMNUMBER(IVariableID,IErr).LT.0.5_RKIND) THEN
     CALL OutofUnitCellCheck(IVariableID,RNegativeMovement,RCorrectedMovement,IErr)
  ELSE
     CALL OutofUnitCellCheck(IVariableID,RPositiveMovement,RCorrectedMovement,IErr)
  END IF

END SUBROUTINE InitialiseAtomicVectorMagnitudes

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE RandomSequence(RRandomSequence,IRandomSequenceLength,ISeedModifier,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Sets up a pseudo random sequence and selects a number
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,Ivalues(1:8), k,IRandomSequenceLength,ISeedModifier
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       seed
  REAL(RKIND),DIMENSION(IRandomSequenceLength) :: &
       RRandomSequence
  
  CALL DATE_AND_TIME(VALUES=Ivalues)

  IValues = IValues*ISeedModifier
!!$  CALL SYSTEM_CLOCK(
  IF (IRandomFLAG.EQ.0) THEN
     CALL RANDOM_SEED(size=k)
     allocate(seed(1:k))
     seed(:) = Ivalues(8)
     CALL RANDOM_SEED(put=seed)
  ELSE
     CALL RANDOM_SEED(size=k)
     ALLOCATE(seed(1:k))
     seed(:) = IFixedSeed*ISeedModifier
     CALL RANDOM_SEED(put=seed)
  END IF
   
  DEALLOCATE(seed)

  CALL RANDOM_NUMBER(RRandomSequence)
  
!!$  RANDOMSEQUENCE = RRandomNumberSequence(IRequestedNumber)
  
END SUBROUTINE  RANDOMSEQUENCE

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REAL(RKIND) FUNCTION RANDOMNUMBER(IRequestedNumber,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Sets up a pseudo random sequence and selects a number
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,values(1:8), k,IRequestedNumber
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       seed
  REAL(RKIND),DIMENSION(IRequestedNumber) :: &
       RRandomNumberSequence
  
  CALL DATE_AND_TIME(values=values)
  
  IF (IRandomFLAG.EQ.0) THEN
     CALL RANDOM_SEED(size=k)
     allocate(seed(1:k))
     seed(:) = values(8)
     CALL RANDOM_SEED(put=seed)
  ELSE
     CALL RANDOM_SEED(size=k)
     ALLOCATE(seed(1:k))
     seed(:) = IFixedSeed*IRequestedNumber
     CALL RANDOM_SEED(put=seed)
  END IF
   
  CALL RANDOM_NUMBER(RRandomNumberSequence)
  
  RANDOMNUMBER = RRandomNumberSequence(IRequestedNumber)
  
END FUNCTION RANDOMNUMBER

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE OutofUnitCellCheck(IVariableID,RProposedMovement,RCorrectedMovement,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Checks that vector movement applied by the simplex initialisation
!!$  % does not move an atom out fo the unit cell, and if it does
!!$  % the atom is moved back into the unit cell on the opposite side
!!$  % as if the atom had moved from one unit cell into the neighbouring one
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       ind,IErr,IVariableID,IAtomID,IVectorID
  REAL(RKIND),DIMENSION(THREEDIM) :: &
       RProposedAtomicCoordinate,RDummyMovement
  REAL(RKIND),INTENT(IN) :: &
       RProposedMovement
  REAL(RKIND),INTENT(OUT) :: &
       RCorrectedMovement
  
  IVectorID = IIterativeVariableUniqueIDs(IVariableID,3)
  
  IAtomID = IAllowedVectorIDs(IVectorID)
 
  RProposedAtomicCoordinate(:) = RAtomSiteFracCoordVec(IAtomID,:) + &
       RProposedMovement*RAllowedVectors(IVectorID,:)

  RDummyMovement = RProposedMovement

  IF(ANY(RProposedAtomicCoordinate.GT.ONE).OR.ANY(RProposedAtomicCoordinate.LT.ZERO)) THEN
     DO ind = 1,THREEDIM
        IF (RProposedAtomicCoordinate(ind).GT.ONE) THEN
           RDummyMovement(ind) = (ONE-RAtomSiteFracCoordVec(IAtomID,ind))/RAllowedVectors(IVectorID,ind)
        ELSEIF(RProposedAtomicCoordinate(ind).LT.ZERO) THEN
           RDummyMovement(ind) = (-RAtomSiteFracCoordVec(IAtomID,ind))/RAllowedVectors(IVectorID,ind)
        ELSE
           RDummyMovement(ind) = RProposedMovement
        END IF
     END DO
  END IF

  IF(RProposedMovement.LT.ZERO) THEN
     RCorrectedMovement = MAXVAL(RDummyMovement)
  ELSE
     RCorrectedMovement = MINVAL(RDummyMovement)
  END IF

END SUBROUTINE OutofUnitCellCheck

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ApplyNewStructureFactors(IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Subroutine to place iteratively determined Structure factors
!!$  % to Ug Matrix
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
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
  COMPLEX(CKIND),DIMENSION(nReflections,nReflections) :: &
       CUgMatDummy

!!$  Dummy Matrix to contain new iterative values
  
   CUgMatDummy = CZERO

!!$  Populate Ug Matrix with new iterative elements

  DO ind = 1,INoofUgs
     WHERE(ISymmetryRelations.EQ.ISymmetryStrengthKey(ind,2)) 
        CUgMatDummy = CSymmetryStrengthKey(ind)
     END WHERE
  END DO

  WHERE(ABS(CUgMatDummy).GT.TINY)
     CUgMat = CUgMatDummy
  END WHERE

!!$  CUgMat now contains the new values from the iterative process 
  
END SUBROUTINE ApplyNewStructureFactors

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CreateIdentityMatrix(IIdentityMatrix,ISize,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % This Subroutine creates an identity matrix of size
!!$  % ISize * ISize
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       IErr,ISize,ind
  INTEGER(IKIND),DIMENSION(ISize,ISize) :: &
       IIdentityMatrix

  IIdentityMatrix = 0

  DO ind = 1,ISize
     IIdentityMatrix(ind,ind) = 1
  END DO

END SUBROUTINE CreateIdentityMatrix

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE RecoverSavedSimplex(RSimplexVolume,RSimplexFoM,RStandardDeviation,RMean,IIterationCount,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % This Subroutine reads the fr-simplex.txt file from a previous
!!$  % refinement run, and recreates the simplex volume and tolerances
!!$  % allowing for the continuation of a previous refinement
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind,IIterationCount
  REAL(RKIND),DIMENSION(IIndependentVariables+1,IIndependentVariables) :: &
       RSimplexVolume
  REAL(RKIND),DIMENSION(IIndependentVariables+1) :: &
       RSimplexFoM
  REAL(RKIND) :: &
       RStandardDeviation,RMean
  CHARACTER*200 :: &
       CSizeofData,SFormatString,filename

  WRITE(filename,*) "fr-Simplex.txt"

  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
        FILE=TRIM(ADJUSTL(filename)))
  
  WRITE(CSizeofData,*) IIndependentVariables+1
  WRITE(SFormatString,*) "("//TRIM(ADJUSTL(CSizeofData))//"(1F6.3,1X),A1)"

  DO ind = 1,(IIndependentVariables+1)
     READ(IChOutSimplex,FMT=SFormatString) RSimplexVolume(ind,:),RSimplexFoM(ind)
  END DO
    
  READ(IChOutSimplex,FMT="(2(1F6.3,1X),I5.1,I5.1,A1)") RStandardDeviation,RMean,IStandardDeviationCalls,IIterationCount

  CLOSE(IChOutSimplex)

END SUBROUTINE RecoverSavedSimplex

SUBROUTINE DetermineNumberofRefinementVariablesPerType(INoofelementsforeachrefinementtype,IErr)

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
  INTEGER(IKIND),DIMENSION(IRefinementVariableTypes) :: &
       INoofelementsforeachrefinementtype

  INoofelementsforeachrefinementtype(1) = &
       IRefineModeSelectionArray(1)*INoofUgs*2
  INoofelementsforeachrefinementtype(2) = &
       IRefineModeSelectionArray(2)*IAllowedVectors
  INoofelementsforeachrefinementtype(3) = &
       IRefineModeSelectionArray(3)*SIZE(IAtomicSitesToRefine)
  INoofelementsforeachrefinementtype(4) = &
       IRefineModeSelectionArray(4)*SIZE(IAtomicSitesToRefine)
  INoofelementsforeachrefinementtype(5) = &
       IRefineModeSelectionArray(5)*SIZE(IAtomicSitesToRefine)*6
  INoofelementsforeachrefinementtype(6) = &
       IRefineModeSelectionArray(6)*3
  INoofelementsforeachrefinementtype(7) = &
       IRefineModeSelectionArray(7)*3
  INoofelementsforeachrefinementtype(8) = &
       IRefineModeSelectionArray(8)
  INoofelementsforeachrefinementtype(9) = &
       IRefineModeSelectionArray(9)

END SUBROUTINE DetermineNumberofRefinementVariablesPerType
