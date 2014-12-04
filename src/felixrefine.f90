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
       ind,IIterationCount
  REAL(RKIND) :: &
       StartTime, CurrentTime, Duration, TotalDurationEstimate,&
       RFigureOfMerit,SimplexFunction  
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RSimplexVolume
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       RSimplexFoM,RIndependentVariableValues
  REAL(RKIND) :: &
       RBCASTREAL

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
  
  ISoftwareMode = 2 ! felixrefinemode
  
  !Read from input files
  CALL ReadInput (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in ReadInput()"
     GOTO 9999
  ENDIF

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
!!$  IF(my_rank.EQ.0) THEN
     ALLOCATE( &
          RSimplexFoM(IIndependentVariables),&
          STAT=IErr)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixrefine (", my_rank, ") error in Allocation()"
        GOTO 9999
     ENDIF
!!$  END IF
    
  CALL SimplexInitialisation(RSimplexVolume,RSimplexFoM,RIndependentVariableValues,1,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in SimplexInitialisation()"
     GOTO 9999
  ENDIF


  !--------------------------------------------------------------------
  ! Apply Simplex Method
  !--------------------------------------------------------------------

!!$  IIterationCount = 0
!!$  IF(my_rank.EQ.0) THEN
!!$     PRINT*,"IIterationCount =",IIterationCount
!!$  END IF

  CALL NDimensionalDownhillSimplex(RSimplexVolume,RSimplexFoM,&
       IIndependentVariables+1,&
       IIndependentVariables,IIndependentVariables,&
       0.001d0,IIterationCount,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in NDimensionalDownhillSimplex()"
     GOTO 9999
  ENDIF

  
  
  IF(my_rank.EQ.0) THEN
     DO ind=1,(IIndependentVariables+1)
        PRINT*,RSimplexVolume(ind,:)
     END DO
     PRINT*,"IIterationCount =",IIterationCount
  END IF

  
!!$  DEALLOCATE(&
!!$       RSimplexFoM,&
!!$       STAT=IErr)  
!!$     IF( IErr.NE.0 ) THEN
!!$        PRINT*,"felixrefine (", my_rank, ") error in Deallocation()"
!!$        GOTO 9999
!!$     ENDIF
!!$     PRINT*,"--------------------DEALLOCATED----------------------------------"

  CALL MPI_BARRIER(MPI_COMM_WORLD,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error ", IErr, " in MPI_BARRIER()"
     STOP
  ENDIF

  PRINT*,"Im rank",my_rank

  STOP
!!$  DEALLOCATE(&
!!$       RSimplexVolume,&
!!$       STAT=IErr)  
!!$  IF( IErr.NE.0 ) THEN
!!$     PRINT*,"felixrefine (", my_rank, ") error in Deallocation()"
!!$     GOTO 9999
!!$  ENDIF

  !--------------------------------------------------------------------
  ! Apply Simplex Method
  !--------------------------------------------------------------------
  
!!$  CALL NDimensionalDownhillSimplex
!!$  IF( IErr.NE.0 ) THEN
!!$     PRINT*,"felixrefine (", my_rank, ") error in SetupSimplexVolume()"
!!$     GOTO 9999
!!$  ENDIF

!!$  RFigureOfMerit = SimplexFunction(IErr)

  !--------------------------------------------------------------------
  ! Deallocate Memory
  !--------------------------------------------------------------------

  
  IF(my_rank.EQ.0) THEN
     PRINT*,"--------------------DEALLOCATING----------------------------------"
  END IF

  PRINT*,my_rank,ALLOCATED(RImageExpi)

  DEALLOCATE( &
       RImageExpi,&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Deallocation()"
     GOTO 9999
  ENDIF

   IF(my_rank.EQ.0) THEN
     PRINT*,"--------------------DEALLOCATED----------------------------------"
  END IF
  
  PRINT*,my_rank,ALLOCATED(RImageExpi)
  !--------------------------------------------------------------------
  ! finish off
  !--------------------------------------------------------------------
  
  CALL cpu_time(CurrentTime)
  Duration=(CurrentTime-StartTime)
  IHours = FLOOR(Duration/3600.0D0)
  IMinutes = FLOOR(MOD(Duration,3600.0D0)/60.0D0)
  ISeconds = MOD(Duration,3600.0D0)-IMinutes*60.0D0
  IMilliSeconds = INT((Duration-(IHours*3600+IMinutes*60+ISeconds))*100,IKIND)


   IF(my_rank.EQ.0) THEN
     PRINT*,"--------------------TIME CALCULATED----------------------------------"
  END IF
  
  PRINT*, "Felixrefine(", my_rank, ") ", RStr, ", used time=", IHours, "hrs ",IMinutes,"mins ",ISeconds,"Seconds ",&
       IMilliSeconds,"Milliseconds"
  !--------------------------------------------------------------------
  ! Shut down MPI
  !--------------------------------------------------------------------
  
   IF(my_rank.EQ.0) THEN
     PRINT*,"--------------------TIME PRINTED----------------------------------"
  END IF

!!$  IF(my_rank.EQ.0) THEN
!!$     CALL MPI_ISEND(RBCASTREAL,1,MPI_DOUBLE_PRECISION,

  CALL MPI_BARRIER(MPI_COMM_WORLD,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error ", IErr, " in MPI_BARRIER()"
     STOP
  ENDIF

  
  IF(my_rank.EQ.0) THEN
     PRINT*,"--------------------BARRIERED----------------------------------"
  END IF

  
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
  INTEGER(IKIND),DIMENSION(7) :: &
       INoofelementsforeachrefinementtype

  INoofelementsforeachrefinementtype(1) = &
       IRefineModeSelectionArray(1)*INoofUgs
  INoofelementsforeachrefinementtype(2) = &
       IRefineModeSelectionArray(2)*SIZE(IAtomicSitesToRefine)*3
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

  ALLOCATE(&
       IIterativeVariableUniqueIDs(IIndependentVariables,5),&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Deallocation()"
     RETURN
  ENDIF

  IIterativeVariableUniqueIDs = 0
  ICalls = 0

  DO ind = 1,7 !Loop over all possible iterative variables
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
  INTEGER(IKIND),DIMENSION(7) :: &
       INoofelementsforeachrefinementtype

!!$  Calculate How Many of Each Variable Type There are

  INoofelementsforeachrefinementtype(1) = &
       IRefineModeSelectionArray(1)*INoofUgs
  INoofelementsforeachrefinementtype(2) = &
       IRefineModeSelectionArray(2)*SIZE(IAtomicSitesToRefine)*3
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
     IArrayToFill(IArrayIndex,3) = IAtomicSitesToRefine(INT(CEILING(REAL(IVariableNo/3.0D0,RKIND))))
     IArrayToFill(IArrayIndex,4) = &
          NINT(3.D0*(REAL(IVariableNo/3.0D0,RKIND)-CEILING(REAL(IVariableNo/3.0D0,RKIND)))+3.0D0)
     IArrayToFill(IArrayIndex,5) = 0

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
             RAtomSiteFracCoordVec(&
             IIterativeVariableUniqueIDs(ind,3),&
             IIterativeVariableUniqueIDs(ind,4))
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
       IErr
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(OUT) :: &
       RIndependentVariableValues
  INTEGER(IKIND),INTENT(IN) :: &
       IIterationCount

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"StructureFactorRefinementSetup(",my_rank,")"
  END IF

  IF(IRefineModeSelectionArray(1).EQ.1.AND.IIterationCount.EQ.1) THEN
     
     RIndependentVariableValues(:INoofUgs) = &
          REAL(CSymmetryStrengthKey(:INoofUgs),RKIND)
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

SUBROUTINE SimplexInitialisation(RSimplexVolume,RSimplexFoM,RIndependentVariableValues,IIterationCount,IErr)
  
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
  REAL(RKIND),DIMENSION(IIndependentVariables+1,IIndependentVariables),INTENT(OUT) :: &
       RSimplexVolume
  REAL(RKIND),DIMENSION(IIndependentVariables+1),INTENT(OUT) :: &
       RSimplexFoM
  REAL(RKIND) :: &
       SimplexFunction,RSimplexDummy
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(IN) :: &
       RIndependentVariableValues
  INTEGER(IKIND),INTENT(IN) :: &
       IIterationCount

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"SimplexInitialisation(",my_rank,")"
  END IF
      
  CALL PerformDummySimulationToSetupSimplexValues(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation(", my_rank, ") error in PerformDummySimulationToSetupSimplexValues()"
     RETURN
  ENDIF
  
  CALL RefinementVariableSetup(RIndependentVariableValues,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation(", my_rank, ") error in RefinementVariableSetup()"
     RETURN
  ENDIF

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

  DEALLOCATE(&
       CUgmat,&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadExperimentalImages (", my_rank, ") error in deAllocation()"
     RETURN
  ENDIF

  DEALLOCATE( &
       CUgMatPrime,&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadExperimentalImages (", my_rank, ") error in deAllocation()"
     RETURN
  ENDIF

  DEALLOCATE( &
       ISymmetryRelations,&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadExperimentalImages (", my_rank, ") error in deAllocation()"
     RETURN
  ENDIF

  DEALLOCATE( &
       ISymmetryStrengthKey,&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadExperimentalImages (", my_rank, ") error in deAllocation()"
     RETURN
  ENDIF

  DEALLOCATE( &
       CSymmetryStrengthKey)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
  
  DO ind = 1,(IIndependentVariables+1)
     RSimplexVolume(ind,:) = &
          RIndependentVariableValues
     IF(ind.GT.1) THEN
        RSimplexVolume(ind,ind-1) = &
             RIndependentVariableValues(ind-1)*1.1
     END IF
  END DO

  DO ind = 1,(IIndependentVariables+1)
          
     IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"---------------------------------------------------------"
        PRINT*,"-------- Simplex",ind,"of",IIndependentVariables+1
        PRINT*,"---------------------------------------------------------"
     END IF

     RSimplexDummy = SimplexFunction(RSimplexVolume(ind,:),1,IErr)
     
     RSimplexFoM(ind) =  RSimplexDummy
     
     IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"---------------------------------------------------------"
        PRINT*,"-------- Simplex",RSimplexVolume(ind,:),RSimplexFoM(ind)
        PRINT*,"---------------------------------------------------------"
     END IF
  END DO
     
  END SUBROUTINE SimplexInitialisation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE PerformDummySimulationToSetupSimplexValues(IErr)
  
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
 
  CALL StructureFactorSetup(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error in StructureFactorSetup()"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       RgMatMat,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation of RgMatMat"
     RETURN
  ENDIF
  
  DEALLOCATE(&
       RgMatMag,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation of RgMatMag"
     RETURN
  ENDIF
  
  DEALLOCATE(&
       RrVecMat,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation of RgMatMag"
     RETURN
  ENDIF
  
  DEALLOCATE( &
       MNP,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
      
  DEALLOCATE( &
        SMNP, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RFullAtomicFracCoordVec, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       SFullAtomicNameVec,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IFullAnisotropicDWFTensor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IFullAtomNumber,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RFullIsotropicDebyeWallerFactor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RFullPartialOccupancy,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RDWF,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       ROcc,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IAtoms,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IAnisoDWFT,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       Rhkl,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RgVecMatT,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RgVecMag,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RGn,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RSg,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       Rhklpositions,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RMask,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IPixelLocations,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"PerformDummySimulationToSetupSimplexValues(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF

  IDiffractionFLAG = 0

END SUBROUTINE PerformDummySimulationToSetupSimplexValues
