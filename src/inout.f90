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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: inout.f90,v 1.59 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! -----------------------------------------------------------------------
!Input: Read the input file
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE Input( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) IErr, ILine,ind,IPos,IPos1,IPos2
  REAL(KIND=RKIND) ROfIter
  CHARACTER*200 SImageMode,SElements
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Input()"
  END IF

  OPEN(UNIT= IChInp, ERR= 120, FILE= "Felix.inp",&
       STATUS= 'OLD')
  ILine= 1

  ! ----------------------------------------------------------------------
  ! introductory comment lines

  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ! ----------------------------------------------------------------------
  ! BLOCH method input

  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)') 
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ! ----------------------------------------------------------------------
  ! control flags

  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IWriteFLAG
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IWriteFLAG = ", IWriteFLAG
  END IF
  
  ILine= ILine+1
  READ(IChInp,FMT='(A)',ERR=20,END=30) SImageMode
  IPos = 0
  IF(SCAN(SImageMode,'0').NE.0) THEN
     IPos = IPos +1
  END IF
  IF(SCAN(SImageMode,'1').NE.0) THEN
     IPos = IPos +1
  END IF
  IF(SCAN(SImageMode,'2').NE.0) THEN
     IPos = IPos +1
  END IF
  
  SELECT CASE (IPos)
  CASE (1)
     IF(SCAN(SImageMode,'2').NE.0) THEN
        IImageFLAG = 3
     ELSE
        IF(SCAN(SImageMode,'1').NE.0) THEN
           IImageFLAG = 1
        ELSE
           IImageFlag = 0
        END IF
     END IF     
  CASE (2)     
     IF(SCAN(SImageMode,'2').NE.0) THEN
        IF(SCAN(SImageMode,'1').NE.0) THEN
           IImageFLAG = 5
        ELSE
           IImageFLAG = 4
        END IF
     ELSE
        IImageFLAG = 2
     END IF
     
  CASE (3)
     IImageFlag = 6
  END SELECT

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IImageFLAG = ", IImageFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IOutputFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IOutputFLAG = ", IOutputFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IBinorTextFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IBinorTextFLAG = ", IBinorTextFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IScatterFactorMethodFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IScatterFactorMethodFLAG = ", IScatterFactorMethodFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) ICentralBeamFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ICentralBeamFLAG = ", ICentralBeamFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMaskFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IMaskFLAG = ", IMaskFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IZolzFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IZolzFLAG = ", IZolzFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IAbsorbFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IAbsorbFLAG = ", IAbsorbFLAG
  ENDIF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IAnisoDebyeWallerFactorFlag
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IAnisoDebyeWallerFactorFlag = ", IAnisoDebyeWallerFactorFlag
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IBeamConvergenceFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IBeamConvergenceFLAG = ", IBeamConvergenceFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IPseudoCubicFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IPseudoCubicFLAG = ", IPseudoCubicFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IXDirectionFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IXDirectionFLAG = ", IXDirectionFLAG
  END IF

  ! ----------------------------------------------------------------------
  ! beam details
  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IPixelCount
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IPixelCount = ",IPixelCount
  END IF

  ! ----------------------------------------------------------------------
  ! beam selection criteria
  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMinReflectionPool
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IMinReflectionPool = ", IMinReflectionPool
  END IF
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMinStrongBeams
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IMinStrongBeams = ", IMinStrongBeams
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMinWeakBeams
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IMinWeakBeams = ", IMinWeakBeams
  END IF

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RBSBMax
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RBSBMax = ", RBSBMax
  END IF

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RBSPMax
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RBSPMax = ", RBSPMax
  END IF

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RConvergenceTolerance
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RConvergenceTolerance = ", RConvergenceTolerance
  END IF

  ! ----------------------------------------------------------------------
  ! crystal settings
  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RDebyeWallerConstant
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RDebyeWallerConstant = ", RDebyeWallerConstant
  END IF

  RMeanSquaredDisplacement= RDebyeWallerConstant/(8*PI**2)
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RAbsorptionPercentage
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RAbsorptionPercentage = ", RAbsorptionPercentage
  END IF

  ! ----------------------------------------------------------------------
  ! microscopy settings

  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RConvergenceAngle
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RConvergenceAngle = ", RConvergenceAngle
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IIncidentBeamDirectionX
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IIncidentBeamDirectionX = ", IIncidentBeamDirectionX
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IIncidentBeamDirectionY
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IIncidentBeamDirectionY = ", IIncidentBeamDirectionY
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IIncidentBeamDirectionZ
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IIncidentBeamDirectionZ = ", IIncidentBeamDirectionZ
  END IF

  RZDirC(1)= REAL(IIncidentBeamDirectionX,RKIND)
  RZDirC(2)= REAL(IIncidentBeamDirectionY,RKIND)
  RZDirC(3)= REAL(IIncidentBeamDirectionZ,RKIND)
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IXDirectionX
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IXDirectionX = ", IXDirectionX
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IXDirectionY
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IXDirectionY = ", IXDirectionY
  END IF
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IXDirectionZ
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IXDirectionZ = ", IXDirectionZ
  END IF

  RXDirC(1)= IXDirectionX
  RXDirC(2)= IXDirectionY
  RXDirC(3)= IXDirectionZ  

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) INormalDirectionX
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"INormalDirectionX = ", INormalDirectionX
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) INormalDirectionY
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"INormalDirectionY = ", INormalDirectionY
  END IF
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) INormalDirectionZ
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"INormalDirectionZ = ", INormalDirectionZ
  END IF

  RNormDirC(1)= INormalDirectionX
  RNormDirC(2)= INormalDirectionY
  RNormDirC(3)= INormalDirectionZ
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RAcceleratingVoltage

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RAcceleratingVoltage = ", RAcceleratingVoltage
  END IF

  ! ----------------------------------------------------------------------
  ! Title Space
  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ! ----------------------------------------------------------------------
  ! Image Output Options

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RInitialThickness
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RInitialThickness = ", RInitialThickness
  END IF

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RFinalThickness
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RFinalThickness = ", RFinalThickness
  END IF

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RDeltaThickness
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RDeltaThickness = ", RDeltaThickness
  END IF
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IReflectOut
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IReflectOut = ", IReflectOut
  END IF

  IF(ISoftwareMode.EQ.2) THEN

     !-----------------------------------------------------------------------
     ! FelixRefine Input
     
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     
     !-----------------------------------------------------------------------
     ! Refinement Specific Flags
     
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     
     ILine= ILine+1
     READ(IChInp,10,ERR=20,END=30) IImageOutputFLAG
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"IImageOutputFLAG = ", IImageOutputFLAG
     END IF
     
     ILine= ILine+1
     READ(IChInp,10,ERR=20,END=30) IDevFLAG
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"IDevFLAG = ", IDevFLAG
     END IF
     
     ILine= ILine+1
     READ(IChInp,10,ERR=20,END=30) IRefineModeFLAG
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"IRefineModeFLAG = ", IRefineModeFLAG
     END IF
     
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     
     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RInitialDebyeWallerFactor
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"RInitialDebyeWallerFactor = ", RInitialDebyeWallerFactor
     END IF
     
     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RFinalDebyeWallerFactor
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"RFinalDebyeWallerFactor = ", RFinalDebyeWallerFactor
     END IF
     
     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RDeltaDebyeWallerFactor
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"RDeltaDebyeWallerFactor = ", RDeltaDebyeWallerFactor
     END IF

     ILine= ILine+1
     READ(IChInp,FMT='(A)',ERR=20,END=30) SElements
     IPos1 = SCAN(SElements,'{')
     IPos2 = SCAN(SElements,'}')
     IPos = SCAN(SElements,'0')
     IElements = 1
     IF(IPos2.EQ.(IPos1+1).OR.IPos.EQ.(IPos1+1).OR.IPos2.EQ.0.OR.IPos1.EQ.0) THEN
        IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"No Elements have been specified, Felix will assume all elements are to be refined = "
        END IF
     ELSE
        DO 
           IPos = SCAN(SElements((IPos1):IPos2),',')
           IPos1 = IPos1+IPos
           IF(IPos.EQ.0) THEN
              EXIT
           ELSE
              IElements = IElements + 1
           END IF
              
        END DO
        
        ALLOCATE(&
             IElementList(IElements),&
             STAT=IErr)
        IF(IErr.NE.0) THEN
           PRINT*,"Input(",my_rank,") ERROR IN ALLOCATE OF IElementList"
           RETURN
        ENDIF
        
        
        IPos1 = SCAN(SElements,'{')
        IPos2 = SCAN(SElements,'}')

        DO ind = 1,IElements

           IPos = SCAN(SElements((IPos1+1):IPos2),',')
           IF(IPos.NE.0) THEN
              READ(SElements((IPos1+1):(IPos1+IPos-1)),FMT='(I3.1)') IElementList(ind) 
              IPos1 = IPos1+IPos
           ELSE
              READ(SElements((IPos1+1):(IPos2-1)),FMT='(I3.1)') IElementList(ind) 
           END IF
        END DO
     END IF

     
     !-----------------------------------------------------------------------
     ! Iterative Ug input
     
     
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     
     ILine= ILine+1
     READ(IChInp,10,ERR=20,END=30) INoofUgs
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"INoofUgs = ", INoofUgs
     END IF
     
     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RLowerBoundUgChange
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"RLowerBoundUgChange = ", RLowerBoundUgChange
     END IF
     
     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RUpperBoundUgChange
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"RUpperBoundUgChange = ", RUpperBoundUgChange
     END IF
     
     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RDeltaUgChange
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"RDeltaUgChange = ", RDeltaUgChange
     END IF
  END IF
  
10 FORMAT(27X,I15.1)
  ! 10	FORMAT("IMAXIteration= ",I15.1)
15 FORMAT(27X,F18.9)
  ! 15     FORMAT("IMAXIteration= ",F18.9)
  
  CLOSE(IChInp, ERR=130)
  
  ! check the parameters for validity
   
  RETURN
  
  !	error in OPEN detected
120 IF(my_rank.EQ.0.OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Input(): ERR in OPEN"
     PRINT*,""
     CALL WriteOutInputFile
  END IF
  GOTO 1000
  
  !	error in CLOSE detected
130 IF(my_rank.EQ.0.OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Input(): ERR in CLOSE"
  END IF
  GOTO 1000
  
  !	error in READ detected
20 IF(my_rank.EQ.0.OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Input(): ERR in READ at line", ILine
  END IF
  GOTO 1000
  
  !	EOF in READ occured prematurely
30 IF(my_rank.EQ.0.OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Input(): EOF in READ at line", ILine
  END IF
  
  ! dump the input help
  
1000 IF((my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.(IWriteFLAG.GE.10.AND.ISoftwareMode.LT.2)) THEN
     PRINT*,"# Input file for FelixSim/Draw version :VERSION: Build :BUILD:"
     PRINT*,"# ------------------------------------"
     PRINT*,""
     PRINT*,"# ------------------------------------"
     PRINT*,"# FelixSim input"
     PRINT*,""
     PRINT*,"# control flags"
     PRINT*,"IWriteFLAG                = 1"
     PRINT*,"IImageFLAG                = 1"
     PRINT*,"IOutputFLAG               = 0"
     PRINT*,"IBinorTextFLAG            = 0"  
     PRINT*,"IScatterFactorMethodFLAG  = 0"
     PRINT*,"ICentralBeamFLAG          = 1"
     PRINT*,"IMaskFLAG                 = 0"
     PRINT*,"IZolzFLAG                 = 1"
     PRINT*,"IAbsorbFLAG               = 1"
     PRINT*,"IAnisoDebyeWallerFlag     = 0"
     PRINT*,"IBeamConvergenceFLAG      = 1"
     PRINT*,"IPseudoCubicFLAG          = 0"
     PRINT*,"IXDirectionFLAG           = 1"
     PRINT*,""
     PRINT*,"# radius of the beam in pixels"
     PRINT*,"IPixelCount               = 64"
     PRINT*,""
     PRINT*,"# beam selection criteria"
     PRINT*,"IMinReflectionPool        = 600"
     PRINT*,"IMinStrongBeams           = 125"
     PRINT*,"IMinWeakBeams             = 20"
     PRINT*,"RBSBMax                   = 0.1"
     PRINT*,"RBSPMax                   = 0.1"
     PRINT*,"RConvergenceTolerance (%) = 1.0"
     PRINT*,""
     PRINT*,"# crystal settings"
     PRINT*,"RDebyeWallerConstant      = 0.4668"
     PRINT*,"RAbsorptionPer            = 2.9"
     PRINT*,""
     PRINT*,"# microscope settings"
     PRINT*,"RConvergenceAngle         = 6.0"
     PRINT*,"IIncidentBeamDirectionX   = 0"
     PRINT*,"IIncidentBeamDirectionY   = 1"
     PRINT*,"IIncidentBeamDirectionZ   = 1"
     PRINT*,"IXDirectionX              = 1"
     PRINT*,"IXDirectionY              = 0"
     PRINT*,"IXDirectionZ              = 0"
     PRINT*,"INormalDirectionX         = 0"
     PRINT*,"INormalDirectionY         = 1"
     PRINT*,"INormalDirectionZ         = 1"
     PRINT*,"RAcceleratingVoltage (kV) = 200.0"
     PRINT*,""
     PRINT*,"# Image Output Options"
     PRINT*,""
     PRINT*,"RInitialThickness        = 300.0"
     PRINT*,"RFinalThickness          = 1300.0"
     PRINT*,"RDeltaThickness          = 10.0"
     PRINT*,"IReflectOut              = 49"
     PRINT*,""
     
     PRINT*,"A Sample Input File Has been Written For you as Felix.inp.SimDraw_sample"
     PRINT*,"It must be renamed to Felix.inp before use"
  ELSE
     PRINT*,"# Input file for FelixRefine version :VERSION: Build :BUILD:"
     PRINT*,"# ------------------------------------"
     PRINT*,""
     PRINT*,"# ------------------------------------"
     PRINT*,"# FelixSim input"
     PRINT*,""
     PRINT*,"# control flags"
     PRINT*,"IWriteFLAG                = 1"
     PRINT*,"IImageFLAG                = 1"
     PRINT*,"IOutputFLAG               = 0"
     PRINT*,"IBinorTextFLAG            = 0"  
     PRINT*,"IScatterFactorMethodFLAG  = 0"
     PRINT*,"ICentralBeamFLAG          = 1"
     PRINT*,"IMaskFLAG                 = 0"
     PRINT*,"IZolzFLAG                 = 1"
     PRINT*,"IAbsorbFLAG               = 1"
     PRINT*,"IAnisoDebyeWallerFlag     = 0"
     PRINT*,"IBeamConvergenceFLAG      = 1"
     PRINT*,"IPseudoCubicFLAG          = 0"
     PRINT*,"IXDirectionFLAG           = 1"
     PRINT*,""
     PRINT*,"# radius of the beam in pixels"
     PRINT*,"IPixelCount               = 64"
     PRINT*,""
     PRINT*,"# beam selection criteria"
     PRINT*,"IMinReflectionPool        = 600"
     PRINT*,"IMinStrongBeams           = 125"
     PRINT*,"IMinWeakBeams             = 20"
     PRINT*,"RBSBMax                   = 0.1"
     PRINT*,"RBSPMax                   = 0.1"
     PRINT*,"RConvergenceTolerance (%) = 1.0"
     PRINT*,""
     PRINT*,"# crystal settings"
     PRINT*,"RDebyeWallerConstant      = 0.4668"
     PRINT*,"RAbsorptionPer            = 2.9"
     PRINT*,""
     PRINT*,"# microscope settings"
     PRINT*,"RConvergenceAngle         = 6.0"
     PRINT*,"IIncidentBeamDirectionX   = 0"
     PRINT*,"IIncidentBeamDirectionY   = 1"
     PRINT*,"IIncidentBeamDirectionZ   = 1"
     PRINT*,"IXDirectionX              = 1"
     PRINT*,"IXDirectionY              = 0"
     PRINT*,"IXDirectionZ              = 0"
     PRINT*,"INormalDirectionX         = 0"
     PRINT*,"INormalDirectionY         = 1"
     PRINT*,"INormalDirectionZ         = 1"
     PRINT*,"RAcceleratingVoltage (kV) = 200.0"
     PRINT*,""
     PRINT*,"# Image Output Options"
     PRINT*,""
     PRINT*,"RInitialThickness        = 300.0"
     PRINT*,"RFinalThickness          = 1300.0"
     PRINT*,"RDeltaThickness          = 10.0"
     PRINT*,"IReflectOut              = 49"
     PRINT*,""
     PRINT*,"# FelixRefine Input"
     PRINT*,""
     PRINT*,"#Refinement Specific Flags"
     PRINT*,"IImageOutputFLAG          = 1"
     PRINT*,"IDevFLAG                  = 0"
     PRINT*,"IRefineModeFLAG           = 0"
     PRINT*,""
     PRINT*,"# Debye Waller Factor Iteration"
     PRINT*,""
     PRINT*,"RInitialDebyeWallerFactor = 0.1"
     PRINT*,"RFinalDebyeWallerFactor = 1.0"
     PRINT*,"RDeltaDebyeWallerFactor = 0.1"
     PRINT*,"IElementsforDWFchange   = 0"
     PRINT*,""
     PRINT*,"# Ug Iteration"
     PRINT*,"INoofUgs                  = 1"
     PRINT*,"RLowerBoundUgChange       = 50.0"
     PRINT*,"RUpperBoundUgChange       = 50.0"
     PRINT*,"RDeltaUgChange            = 50.0"
     PRINT*,""
     

     PRINT*,"A Sample Input File Has been Written For you as Felix.inp.Refine_sample"
     PRINT*,"It must be renamed to Felix.inp before use"
  END IF
  IErr= 1
  RETURN
  
END SUBROUTINE Input
   

! -----------------------------------------------------------------------
!
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE InputScatteringFactors( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  CHARACTER*200 dummy
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       rdummy

  INTEGER IErr, ILine, ILength, ind
  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN

     PRINT*,"InputScatteringFactors()"
     
  END IF
  ILine= 0
  
!!$  OPEN(UNIT= IChInp, ERR= 120, FILE= "FelixDoyle.sca",&
!!$       STATUS= 'OLD')
  OPEN(UNIT= IChInp, ERR= 120, FILE= "Felix.sca",&
       STATUS= 'OLD')
!!$  
  DO
     READ(UNIT= IChInp, END=100, ERR=20, FMT='(a)') dummy
     ILine=ILine+1
  ENDDO
100 ILength=ILine
  
  IF(IWriteFLAG.GE.10) THEN
     
     PRINT*,"InputScatteringFactors(): ILength=", ILength
     
  END IF
  
  REWIND(UNIT=IChInp)
  IF(IWriteFLAG.GE.10) THEN
     
     PRINT*,"actual reading of data"
  END IF

  SELECT CASE (IScatterFactorMethodFLAG)

  CASE(0)

     ALLOCATE( &
          rdummy(1),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"InputScatteringFactors(): error in memory ALLOCATE()"
        RETURN
     ENDIF

     ALLOCATE(RScattFactors(ILength,12), STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"InputScatteringFactors(): error in memory ALLOCATE()"
        RETURN
     ENDIF

     DO ILine=1,ILength,1
        READ(UNIT= IChInp, FMT='(13(E15.11,1X))') &
             rdummy, RScattFactors(ILine,:)     
        !PRINT*,rdummy, RScattFactors(ILine,:)
     ENDDO

  CASE(1)

     ALLOCATE( &
          rdummy(13),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"InputScatteringFactors(): error in memory ALLOCATE()"
        RETURN
     ENDIF

     ALLOCATE(RScattFactors(ILength,8), STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"InputScatteringFactors(): error in memory ALLOCATE()"
        RETURN
     ENDIF
     
     DO ILine=1,ILength,1
        READ(UNIT= IChInp, FMT='(21(E15.11,1X))') &
             rdummy(:), RScattFactors(ILine,:)     
        !PRINT*,rdummy, RScattFactors(ILine,:)
     ENDDO

  CASE(2)

     ALLOCATE( &
          rdummy(21),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"InputScatteringFactors(): error in memory ALLOCATE()"
        RETURN
     ENDIF

     ALLOCATE(RScattFactors(ILength,8), STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"InputScatteringFactors(): error in memory ALLOCATE()"
        RETURN
     ENDIF
     
     DO ILine=1,ILength,1
        READ(UNIT= IChInp, FMT='(29(E15.11,1X))') &
             rdummy(:), RScattFactors(ILine,:)     
        !PRINT*,rdummy, RScattFactors(ILine,:)
     ENDDO

  END SELECT

  DEALLOCATE( &
       rdummy,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"InputScatteringFactors(): error in memory DEALLOCATE()"
     RETURN
  ENDIF  
  
  CLOSE(IChInp)
  
  RETURN
  
  !	error in OPEN detected
120 PRINT*,"InputScatteringFactors(): ERR in OPEN"
  GOTO 1000
  
  !	error in CLOSE detected
130 PRINT*,"InputScatteringFactors(): ERR in CLOSE"
  GOTO 1000
  
  !	error in READ detected
20 PRINT*,"InputScatteringFactors(): ERR in READ at line", ILine
  GOTO 1000
  
  !	EOF in READ occured prematurely
30 PRINT*,"InputScatteringFactors(): EOF in READ at line", ILine
  
1000 &
  PRINT*,"InputScatteringFactors expects a file such as"
  PRINT*,"--------------------------------------------------------------------"
  PRINT*,"22      0.737291729     9.96300042      0.355417565     9.96300042 ", &
       "0.496982599     0.072659586     0.016335034     0.360199179     1.42171303 ", &
       "     0.073173566     0.957232308      15.8512114"
  PRINT*,"23      0       0       0       0       0       0       0       0  ", &
       "     0       0       0       0"

  IErr= 1
  RETURN
END SUBROUTINE InputScatteringFactors
  
! -----------------------------------------------------------------------
!
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE ReadEigenSystemChunk( IAllocationChunk,IErr )
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara 
  USE IChannels

  USE MPI
  USE MyMPI
    
  IMPLICIT NONE

  CHARACTER*25 EIGENFORMAT
  INTEGER(IKIND) :: IRank, &
       Iindex, Ijndex, IWriteLine, IStrongBeamIndex,IInputBeams, &
       IAllocationChunk,InChunks,IindPrevious, IjndPrevious, &
       IErr, ILength, ind,IReadLine,INewLineFLAG,index
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: &
       REigenVectorRealTemp, REigenVectorImagTemp
  REAL(RKIND) :: &
       REigenValueRealTemp, REigenValueImagTemp

  PRINT*,"InputEigenSystem()"

  !-----------------------------------------------------------
  ! eigenvalues
  !-----------------------------------------------------------

  InBeams = 0
  ind = 0

  IF(IWriteFLAG.GE.10) THEN
     
     PRINT*,"actual reading of data"
  
  END IF
  IReadLine = 1
  INewLineFLAG = 1
  DO 
     READ(UNIT= IChInp, END=100, ERR=20, FMT='(7(I3.1,1X))',ADVANCE='NO') &
          IRank,Iindex, Ijndex, nReflections,&
          IInputBeams, IWriteLine, IStrongBeamIndex

     IF (IInputBeams.EQ.0) THEN
        IReadLine = IReadLine-1
        EXIT
     END IF
        
     InBeams(IReadLine) = IInputBeams
     ILACBEDStrongBeamList(IReadLine,IWriteLine) = IStrongBeamIndex
     IPixelLocation(IReadLine,1) = Iindex
     IPixelLocation(IReadLine,2) = Ijndex
      
     IF(INewLineFlag.EQ.1) THEN
        IPixelCountTotal = IPixelCountTotal + 1
        ALLOCATE( &
             REigenVectorRealTemp(IInputBeams),&
             REigenVectorImagTemp(IInputBeams),&
             STAT=IErr)
        IF(IErr.NE.0) THEN
           PRINT*,"ERROR IN ALLOCATE OF LACBED Temps"
           STOP
        ENDIF
        INewLineFlag=0
     ENDIF
     
     WRITE(EIGENFORMAT,*) IInputBeams

     READ(ICHInp,END=100,ERR=20,FMT="((1F13.10,1X),(1F13.10,1X),&
          "//TRIM(ADJUSTL(TRIM(EIGENFORMAT)))//"(1F13.10,1X), &
          "//TRIM(ADJUSTL(TRIM(EIGENFORMAT)))//"(1F13.10,1X))", &
          ADVANCE='YES') REigenValueRealTemp, REigenValueImagTemp,&     
          REigenVectorRealTemp,REigenVectorImagTemp
     CEigenVectorsChunk(IReadLine,IWriteLine,1:IInputBeams) = &
          REigenVectorRealTemp+REigenVectorImagTemp*CIMAGONE
     CEigenValuesChunk(IReadLine,IWriteLine) = &
          REigenValueRealTemp+ REigenValueImagTemp*CIMAGONE
     
     IF(IReadLine.EQ.IAllocationChunk.AND.IWriteLine.EQ.IInputBeams) EXIT
     IF(IWriteLine.EQ.IInputBeams) THEN
        IReadLine =  IReadLine + 1
        INewLineFLAG = 1
        DEALLOCATE(&
             REigenVectorRealTemp,&
             REigenVectorImagTemp)
     ENDIF
     
  ENDDO

  100 PRINT*,"Done with the reading already"

  RETURN
  
  !	error in OPEN detected
120 PRINT*,"InputEigenSystem(): ERR in OPEN"
  !GOTO 1000
  
  !	error in CLOSE detected
130 PRINT*,"InputEigenSystem(): ERR in CLOSE"
  !GOTO 1000
  
  !	error in READ detected
20 PRINT*,"InputEigenSystem(): ERR in READ at line", IReadLine
  !GOTO 1000
  
  !	EOF in READ occured prematurely
30 PRINT*,"InputEigenSystem(): EOF in READ at line", IReadLine
  
  IErr= 1
  RETURN
END SUBROUTINE ReadEigenSystemChunk
  
! --------------------------------------------------------------------
! OpenData
! --------------------------------------------------------------------


SUBROUTINE OpenData(IChOutWrite, prefix, surname, IErr)

  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  CHARACTER*27 surname, surnamelength
  CHARACTER*2 prefix,postfix
  INTEGER(KIND=IKIND) IChOutWrite, IErr

  CHARACTER*34 filename
  INTEGER index
  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN

     PRINT*,"OpenData()"

  END IF
  
  WRITE(surnamelength,*) LEN_TRIM(surname)

  WRITE(filename,"(A2,A2,A1,A"//TRIM(surnamelength)//",A4)") "F-",prefix,"-",surname,".txt"

  IF(IWriteFLAG.GE.10) THEN
     
     PRINT*,filename

  END IF
  IF (IWriteFLAG.GE.10) THEN
     SELECT CASE(IChOutWrite)
     CASE(IChOutWF)
        PRINT*, "OpenData: opening channel", IChOutWF, &
             "for WAVE FUNCTIONS (WF*.txt)"
     CASE(IChOutWI)
        PRINT*, "OpenData: opening channel", IChOutWI, &
             "for WAVE INTENSITIES (WI*.txt)"
     CASE(IChOutEV)
        PRINT*, "OpenData: opening channel", IChOutEV, &
             "for EIGENVALUES of UgMat (EV*.txt)"
     CASE(IChOutEX)
        PRINT*, "OpenData: opening channel", IChOutEX, &
             "for EIGENVECTORS of UgMat (EX*.txt)"
     CASE(IChOutUM)
        PRINT*, "OpenData: opening channel", IChOutUM, &
             "for UgMat (UM*.txt)"
     CASE DEFAULT
        PRINT*, "OpenData: opening UNKNOWN", IChOutWrite, &
             "channel "
     END SELECT
  END IF

  OPEN(UNIT= IChOutWrite, ERR= 10, STATUS= 'UNKNOWN', FILE=TRIM(filename))

  RETURN

  ! error in OPEN detected
10 PRINT*,"WriteDataC(): ERR in OPEN()"
  !PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenData
  
! --------------------------------------------------------------------
! OpenData
! --------------------------------------------------------------------


SUBROUTINE OpenImageForReadIn(IErr,filename)

  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(KIND=IKIND) IChOutWrite, IErr

  CHARACTER*34 filename

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN

     PRINT*,"OpenImageForReadIn()"

  END IF

  !filename = "Felix.img"

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,filename

  END IF

  OPEN(UNIT= IChInImage, ERR= 10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),FORM='UNFORMATTED',&
       ACCESS='DIRECT',IOSTAT=Ierr,RECL=IImageSizeXY(1)*8)

  RETURN

  ! error in OPEN detected
10 PRINT*,"WriteDataC(): ERR in OPEN()"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenImageForReadIn

SUBROUTINE ReadImageForRefinement(IErr)

  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN

     PRINT*,"ReadImageForRefinement()"

  END IF

  DO ind=1,IImageSizeXY(2)
     READ(IChInImage,rec=ind) RImageIn(ind,:)
  END DO
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadImageForRefinement (", my_rank, ") error in READ()",IErr
     RETURN
  ENDIF
  
END SUBROUTINE ReadImageForRefinement


! --------------------------------------------------------------------
! OpenData
! --------------------------------------------------------------------

SUBROUTINE OpenDataForAppend(IChOutWrite, prefix, surname, IErr)

  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  CHARACTER*27 surname, surnamelength
  CHARACTER*2 prefix,postfix
  INTEGER(KIND=IKIND) IChOutWrite, IErr

  CHARACTER*34 filename
  INTEGER index

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"OpenData()"

  END IF
  
  WRITE(surnamelength,*) LEN_TRIM(surname)
  
  WRITE(filename,"(A2,A2,A1,A"//TRIM(surnamelength)//",A4)") "F-",prefix,"-",surname,".txt"

  IF (IWriteFLAG.GE.10) THEN
     SELECT CASE(IChOutWrite)
     CASE(IChOutWF)
        PRINT*, "OpenData: opening channel", IChOutWF, &
             "for WAVE FUNCTIONS (WF*.txt)"
     CASE(IChOutWI)
        PRINT*, "OpenData: opening channel", IChOutWI, &
             "for WAVE INTENSITIES (WI*.txt)"
     CASE(IChOutEV)
        PRINT*, "OpenData: opening channel", IChOutEV, &
             "for EIGENVALUES of UgMat (EV*.txt)"
     CASE(IChOutEX)
        PRINT*, "OpenData: opening channel", IChOutEX, &
             "for EIGENVECTORS of UgMat (EX*.txt)"
     CASE(IChOutUM)
        PRINT*, "OpenData: opening channel", IChOutUM, &
             "for UgMat (UM*.txt)"
     CASE DEFAULT
        PRINT*, "OpenData: opening UNKNOWN", IChOutWrite, &
             "channel "
     END SELECT
  END IF
  
  OPEN(UNIT= IChOutWrite, ERR= 10, STATUS= 'UNKNOWN',&
       FILE=TRIM(filename),ACCESS='APPEND')

  RETURN

  ! error in OPEN detected
10 PRINT*,"WriteDataC(): ERR in OPEN()"
  !PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenDataForAppend

! --------------------------------------------------------------------
! Open Reflection Image
! --------------------------------------------------------------------

SUBROUTINE OpenReflectionImage(IChOutWrite, surname, IErr,IReflectWriting,IImageSizeX)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  USE IChannels

  CHARACTER*27 surname
  CHARACTER*20 prefix,postfix,h,k,l
  INTEGER(KIND=IKIND) IChOutWrite, IErr,IReflectWriting,IImageSizeX

  CHARACTER*50 filename
  CHARACTER*40 fileext
  INTEGER index

  IF ((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"OpenReflectionImage()"
  END IF
  

  SELECT CASE(IChOutWrite)
  CASE(MontageOut)
  CASE DEFAULT
     IF(IHKLSelectFLAG.EQ.0) THEn
        WRITE(h,*)  NINT(RHKL(IReflectWriting,1))
        WRITE(k,*)  NINT(RHKL(IReflectWriting,2))
        WRITE(l,*)  NINT(RHKL(IReflectWriting,3))
     ELSE
        
        WRITE(h,*)  NINT(RHKL(IOutPutReflections(IReflectWriting),1))
        WRITE(k,*)  NINT(RHKL(IOutPutReflections(IReflectWriting),2))
        WRITE(l,*)  NINT(RHKL(IOutPutReflections(IReflectWriting),3))
     END IF
  END SELECT

  IF (IWriteFLAG.GE.10) THEN
     PRINT*,filename
  END IF
  
  SELECT CASE (IBinorTextFLAG)
  CASE(0)
     WRITE(fileext,*) TRIM(ADJUSTL(".bin")) 
  CASE(1)
     WRITE(fileext,*) TRIM(ADJUSTL(".txt"))
  END SELECT
  
  SELECT CASE(IChOutWrite)
  CASE(IChOutWFImageReal)        
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"/F-WF-A_",&
          TRIM(ADJUSTL(h)),&
          TRIM(ADJUSTL(k)),&
          TRIM(ADJUSTL(l)),&
          TRIM(ADJUSTL(fileext))
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE FUNCTION REAL PART (WR*.txt)"
     END IF
  CASE(IChOutWFImagePhase)        
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"/F-WF-P_",&
          TRIM(ADJUSTL(h)),&
          TRIM(ADJUSTL(k)),&
          TRIM(ADJUSTL(l)),&
          TRIM(ADJUSTL(fileext))
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE FUNCTION PHASE PART (WP*.txt)"
     END IF
  CASE(IChOutWIImage)        
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"/F-WI_",&
          TRIM(ADJUSTL(h)),&
          TRIM(ADJUSTL(k)),&
          TRIM(ADJUSTL(l)),&
          TRIM(ADJUSTL(fileext))
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE INTENSITIES"
     END IF
  CASE(MontageOut)        
     WRITE(filename,*) "F-WI-",TRIM(ADJUSTL(surname)),&
          TRIM(ADJUSTL(fileext))
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE INTENSITIES"
     END IF
  CASE DEFAULT
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening UNKNOWN channel ", IChOutWrite
     END IF
  END SELECT
  
  
  SELECT CASE (IBinorTextFLAG)
     CASE(0)
        OPEN(UNIT=IChOutWrite, ERR=10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),FORM='UNFORMATTED',&
             ACCESS='DIRECT',IOSTAT=Ierr,RECL=IImageSizeX*8)
        
     CASE(1)
        OPEN(UNIT=IChOutWrite, ERR=10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)))
     END SELECT
     RETURN
     
  ! error in OPEN detected
10 PRINT*,"OpenReflectionImage(): ERR in OPEN()",IErr
  IErr= 1
  RETURN
  
END SUBROUTINE OpenReflectionImage

!-----------------------------------------------------------------
! Write Reflection Images
!-----------------------------------------------------------------

SUBROUTINE WriteReflectionImage( IChOutWrite, data, IErr,IImageSizeX,IImageSizeY)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  INTEGER(KIND=IKIND) IErr,IImageSizeY,IImageSizeX
  REAL(KIND=RKIND),DIMENSION(IImageSizeX,IImageSizeY) :: data
  CHARACTER*100 CSizeofData
  INTEGER ind, IChOutWrite
  CHARACTER*100 SFormatString

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"WriteReflectionImage"
  END IF

  SELECT CASE (IBinorTextFLAG)
     
  CASE(0)
     
     DO ind = 1,(IImageSizeY)
        WRITE(IChOutWrite,rec=ind) data(ind,:)
     END DO

  CASE(1)
     
     DO ind = 1,(2*IPixelCount)
        WRITE(CSizeofData,*) 2*IPixelCount
        WRITE(SFormatString,*) "("//TRIM(ADJUSTL(CSizeofData))//"(1F6.3,1X),A1)"
        WRITE(IChOutWrite,FMT=SFormatString,ERR=20) data(ind,:)
     END DO
     
  END SELECT

  RETURN
  ! error in WRITE detected
20 PRINT*,"WriteReflectionImage(): ERR in WRITE()",Ierr
  IErr= 1
  RETURN

  ! buffer too short
30 PRINT*,"WriteImageR_MPI(): EOF in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteReflectionImage

! --------------------------------------------------------------------
! WriteDataC

SUBROUTINE WriteDataC( IChOutWrite, ipos,jpos, data, size, step, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MyMPI

  INTEGER(KIND=IKIND) ipos,jpos, size, step, IErr
  COMPLEX(KIND=CKIND) data(size)

  INTEGER ind, IChOutWrite
  DO ind=1,size,step
     IF (ABS(data(ind)).GE. TINY) THEN
        WRITE(IChOutWrite,100) my_rank, ipos,jpos,ind, DBLE(data(ind)),AIMAG(data(ind)), &
             ARG(DBLE(data(ind)),AIMAG(data(ind)))
     ELSE
        WRITE(IChOutWrite,100) my_rank, ipos,jpos,ind, DBLE(data(ind)),AIMAG(data(ind)), &
             0.0D0
     ENDIF
     
100  FORMAT(4I4,3(G25.15))
  ENDDO
  
  RETURN
  
  ! error in WRITE detected
20 PRINT*,"WriteDataC(): ERR in WRITE()"
  IErr= 1
  RETURN
  
END SUBROUTINE WriteDataC

!---------------------------------------------------------------------
!Output WriteEigenSystem_MPI
!---------------------------------------------------------------------

SUBROUTINE WriteEigenSystem_MPI( IChOutWrite, &
     ipos,jpos,nReflect,nbeamout, CdataEVal,CdataEVec,ISbeamlist, rows, cols, step, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  INTEGER(KIND=IKIND) ipos,jpos, step,rows,cols ,IErr,nReflect,nbeamout
  COMPLEX(KIND=CKIND), DIMENSION(nbeamout,nbeamout) :: CdataEVec
  COMPLEX(CKIND), DIMENSION(nbeamout) :: CdataEVal
  INTEGER(IKIND), DIMENSION(nbeamout) :: ISbeamlist

  INTEGER my_status(MPI_STATUS_SIZE)

  CHARACTER*25 FORMATstring
  CHARACTER(IMAXCBuffer) DATAstring

  INTEGER ind, IChOutWrite
  
  IF ( 2*13*SIZE(CdataEVec)+2*13*SIZE(CdataEVal)+3*SIZE(ISbeamlist)+3*6*ADD_OUT_INFO .GT. IMAXCBuffer ) THEN
     IErr=1
     PRINT*, "WriteEigenSystem_MPI(", my_rank, ") error ", IErr, &
          ", output buffer size IMAXCBuffer=", IMAXCBuffer, &
          " smaller than SIZEOF(EVec+EVal+Sbeam)=", &
          SIZEOF(CdataEVec)+SIZEOF(CdataEVal)+SIZEOF(ISbeamlist), &
          SIZE(CdataEVec) + SIZE(CdataEVal) + SIZE(ISbeamlist), &
          2*13*SIZE(CdataEVec)+2*13*SIZE(CdataEVal)+3*SIZE(ISbeamlist)+3*6*ADD_OUT_INFO
     GOTO 20
  ENDIF
  
  DO ind=1,rows,step
     
     WRITE(FORMATstring,*) nbeamout
     WRITE(DATAstring,&
          "(7(I3.1,1X),"//TRIM(ADJUSTL(TRIM(FORMATstring)))//"(1F13.10,1X), &
          "//TRIM(ADJUSTL(TRIM(FORMATstring)))//"(1F13.10,1X),(1F13.10,1X),(1F13.10,1X),A1)") &
          my_rank, ipos,jpos,nReflect,nbeamout,ind, ISbeamlist(ind), &
          REAL(CdataEVal(ind)), AIMAG(CdataEVal(ind)),&
          REAL(CdataEVec(ind,:)), AIMAG(CdataEVec(ind,:)), &
          CHAR(10)
     
     CALL MPI_FILE_WRITE_SHARED(IChOutWrite, TRIM(DATAstring), LEN_TRIM(DATAstring), &
          MPI_CHARACTER, my_status, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteEigenSystem_MPI(", my_rank, ") error ", IErr, &
             " in MPI_FILE_WRITE_SHARED() for file handle ",IChOutWrite
        RETURN
     ENDIF
     
  ENDDO
  
  RETURN

  ! error in WRITE detected
20 PRINT*,"WriteEigenSystem_MPI(): ERR in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteEigenSystem_MPI
!---------------------------------------------------------------------
!Output WriteEigenSystem_MPI
!---------------------------------------------------------------------

SUBROUTINE WriteEigenSystemBinary_MPI( IChOutWrite, &
     ipos,jpos,nReflect,nbeamout, CdataEVal,CdataEVec,ISbeamlist, rows, cols, step, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  INTEGER(KIND=IKIND) ipos,jpos, step,rows,cols ,IErr,nReflect,nbeamout
  COMPLEX(KIND=CKIND), DIMENSION(nbeamout,nbeamout) :: CdataEVec
  COMPLEX(CKIND), DIMENSION(nbeamout) :: CdataEVal
  INTEGER(IKIND), DIMENSION(nbeamout) :: ISbeamlist

  INTEGER my_status(MPI_STATUS_SIZE)

  CHARACTER*25 FORMATstring
  CHARACTER(IMAXCBuffer) DATAstring
  
  INTEGER ind, IChOutWrite
  
  IF ( 2*13*SIZE(CdataEVec)+2*13*SIZE(CdataEVal)+3*SIZE(ISbeamlist)+3*6*ADD_OUT_INFO .GT. IMAXCBuffer ) THEN
     IErr=1
     PRINT*, "WriteEigenSystem_MPI(", my_rank, ") error ", IErr, &
          ", output buffer size IMAXCBuffer=", IMAXCBuffer, &
          " smaller than SIZEOF(EVec+EVal+Sbeam)=", &
          SIZEOF(CdataEVec)+SIZEOF(CdataEVal)+SIZEOF(ISbeamlist), &
          SIZE(CdataEVec) + SIZE(CdataEVal) + SIZE(ISbeamlist), &
          2*13*SIZE(CdataEVec)+2*13*SIZE(CdataEVal)+3*SIZE(ISbeamlist)+3*6*ADD_OUT_INFO
     GOTO 20
  ENDIF

  DO ind=1,rows,step

     CALL MPI_FILE_WRITE_SHARED(IChOutWrite, TRIM(DATAstring), LEN_TRIM(DATAstring), &
          MPI_CHARACTER, my_status, IErr)

     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteEigenSystem_MPI(", my_rank, ") error ", IErr, &
             " in MPI_FILE_WRITE_SHARED() for file handle ",IChOutWrite
        RETURN
     ENDIF

  ENDDO
  
  RETURN

  ! error in WRITE detected
20 PRINT*,"WriteEigenSystem_MPI(): ERR in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteEigenSystemBinary_MPI

! --------------------------------------------------------------------
! WriteDataR
!--------------------------------------------------------------------

SUBROUTINE WriteDataR( IChOutWrite, ipos,jpos, data, size, step, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MyMPI

  INTEGER(KIND=IKIND) ipos,jpos, size, step, IErr
  REAL(KIND=RKIND) data(size)

  INTEGER ind, IChOutWrite

  DO ind=1,size,step
     WRITE(IChOutWrite,100) my_rank, ipos,jpos,ind,data(ind)
100  FORMAT(4I4,G25.15)
  ENDDO
     
  RETURN

  ! error in WRITE detected
20 PRINT*,"WriteDataR(): ERR in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteDataR

! --------------------------------------------------------------------
! OpenData_MPI

SUBROUTINE OpenDataForAppend_MPI(IChOutWrite, prefix, surname, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE IChannels

  USE MPI
  USE MyMPI

  CHARACTER*27 surname
  CHARACTER*2 prefix,postfix
  INTEGER IChOutWrite, IErr

  CHARACTER*35 filename
  INTEGER index

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"OpenDataForAppend_MPI()"

  END IF

  WRITE(filename,'(A2,A2,A1,A12,A4)') "F-",prefix,"-",surname,".txt"
  
  CALL MPI_FILE_OPEN( MPI_COMM_WORLD, TRIM(filename), &
       MPI_MODE_CREATE+MPI_MODE_APPEND, MPI_INFO_NULL, IChOutWrite, IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"OpenDataForAppend_MPI(", my_rank, ") error ", IErr, &
          " in MPI_FILE_OPEN()"
     STOP
  ENDIF
  
  IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenData: opened channel", IChOutWrite, &
             "for ", filename
  END IF

  RETURN

  ! error in OPEN detected
10 PRINT*,"OpenDataForAppend_MPI(): ERR in OPEN()"
  !PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenDataForAppend_MPI

! --------------------------------------------------------------------
! OpenData_MPI

SUBROUTINE OpenData_MPI(IChOutWrite, prefix, surname, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE IChannels

  USE MPI
  USE MyMPI

  CHARACTER*27 surname
  CHARACTER*2 prefix,postfix
  INTEGER IChOutWrite, IErr

  CHARACTER*35 filename
  INTEGER index

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"OpenData_MPI()"
  END IF

  WRITE(filename,'(A2,A2,A1,A12,A4)') "F-",prefix,"-",surname,".txt"
  
  CALL MPI_FILE_OPEN( MPI_COMM_WORLD, TRIM(filename), &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, IChOutWrite, IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"OpenDataMPI(", my_rank, ") error ", IErr, &
          " in MPI_FILE_OPEN()"
     STOP
  ENDIF
  
  IF (IWriteFLAG.GE.10) THEN
     PRINT*, "OpenData: opened channel", IChOutWrite, &
          "for ", filename
  ENDIF
  RETURN
  
  ! error in OPEN detected
10 PRINT*,"OpenDatMPI(): ERR in OPEN()"
  !PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenData_MPI

! --------------------------------------------------------------------
! OpenImage_MPI

SUBROUTINE OpenImage_MPI(IChOutWrite, surname, IErr,IReflectWriting)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE IChannels

  USE MPI
  USE MyMPI

  CHARACTER*27 surname
  CHARACTER*2 prefix,postfix
  INTEGER IChOutWrite, IErr

  CHARACTER*100 filename,h,k,l
  INTEGER index,IReflectWriting

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"OpenImage_MPI()"
     
  END IF
  
  WRITE(h,*)  NINT(RHKL(IReflectWriting,1))
  WRITE(k,*)  NINT(RHKL(IReflectWriting,2))
  WRITE(l,*)  NINT(RHKL(IReflectWriting,3))

  WRITE(filename,*) TRIM(ADJUSTL(surname)),"/F-WI_",&
       TRIM(ADJUSTL(h)),&
       TRIM(ADJUSTL(k)),&
       TRIM(ADJUSTL(l)),&
       ".txt"
  
  CALL MPI_FILE_OPEN( MPI_COMM_WORLD, TRIM(filename), &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, IChOutWrite, IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"OpenImage_MPI(", my_rank, ") error ", IErr, &
          " in MPI_FILE_OPEN()"
     STOP
  ENDIF
  
  IF (IWriteFLAG.GE.10) THEN
     PRINT*, "OpenData: opened channel", IChOutWrite, &
          "for ", filename
  END IF
  
  RETURN
  
  ! error in OPEN detected
10 PRINT*,"OpenDatMPI(): ERR in OPEN()"
  !PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenImage_MPI

! --------------------------------------------------------------------
! WriteDataRMPI

SUBROUTINE WriteDataR_MPI( IChOutWrite, ipos,jpos, data, size, step, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  INTEGER(KIND=IKIND) ipos,jpos, size, step, IErr
  INTEGER my_status(MPI_STATUS_SIZE)
  REAL(KIND=RKIND) data(size)

  CHARACTER(IMAXRBuffer) outstring
  INTEGER ind, IChOutWrite

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"WriteDataR_MPI()"
  END IF
  ! line by line

  IF ( SIZEOF(data)+ADD_OUT_INFO .GT. IMAXRBuffer ) THEN
     IErr=1
     PRINT*, "WriteDataR_MPI(", my_rank, ") error ", IErr, &
          ", output buffer size IMAXRBuffer=", IMAXRBuffer, &
          " smaller than SIZEOF(data)=", SIZEOF(data)
     GOTO 20
  ENDIF

  WRITE(outstring,*,ERR=20) my_rank,ipos,jpos,size,(data(ind),ind=1,size), CHAR(10)

  CALL MPI_FILE_WRITE_SHARED(IChOutWrite, TRIM(outstring), LEN_TRIM(outstring), &
       MPI_CHARACTER, my_status, IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteDataRMPI(", my_rank, ") error ", IErr, &
          " in MPI_FILE_WRITE_SHARED()"
     STOP
  ENDIF
  
  GOTO 9

  ! element by element
  DO ind=1,size

     CALL MPI_FILE_WRITE_SHARED(IChOutWrite, TRIM(outstring), LEN_TRIM(outstring), &
          MPI_CHARACTER, my_status, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteDataR_MPI(", my_rank, ") error ", IErr, &
             " in MPI_FILE_WRITE_SHARED() for file handle ",IChOutWrite
        RETURN
     ENDIF
  END DO

9 RETURN

  ! error in WRITE detected
20 PRINT*,"WriteDataR_MPI(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! buffer too short
30 PRINT*,"WriteDataR_MPI(): EOF in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteDataR_MPI! --------------------------------------------------------------------
! WriteDataRMPI

SUBROUTINE WriteImageR_MPI( IChOutWrite, data, IErr,ILocalPixelCountMin,ILocalPixelCountMax)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  INTEGER(KIND=IKIND) IErr,ILocalPixelCountMin,ILocalPixelCountMax
  INTEGER my_status(MPI_STATUS_SIZE)
  REAL(KIND=RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: data

  CHARACTER(IMAXRBuffer) outstring,CSizeofData
  INTEGER ind, IChOutWrite
  CHARACTER*100 SFormatString
  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"WriteImageR_MPI"
  END IF

  IF ( 6*SIZE(data)+ADD_OUT_INFO .GT. IMAXRBuffer ) THEN
     IErr=1
     PRINT*, "WriteImageR_MPI(", my_rank, ") error ", IErr, &
          ", output buffer size IMAXRBuffer=", IMAXRBuffer, &
          " smaller than SIZEOF(data)=", SIZE(data)
  ENDIF

  DO ind = 1,(2*IPixelCount)
     WRITE(CSizeofData,*) 2*IPixelCount
     WRITE(SFormatString,*) "("//TRIM(ADJUSTL(CSizeofData))//"(1F6.3,1X),A1)"
     WRITE(outstring,FMT=SFormatString,ERR=20) data(ind,:), CHAR(10)

     CALL MPI_FILE_WRITE_SHARED(IChOutWrite,TRIM(outstring), &
          LEN_TRIM(outstring),MPI_CHARACTER, my_status, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteImageRMPI(", my_rank, ") error ", IErr, &
             " in MPI_FILE_WRITE_SHARED()"
        STOP
     ENDIF
  END DO
  
  RETURN
  ! error in WRITE detected
20 PRINT*,"WriteImageR_MPI(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! buffer too short
30 PRINT*,"WriteImageR_MPI(): EOF in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteImageR_MPI

! --------------------------------------------------------------------
! WriteDataCMPI

SUBROUTINE WriteDataC_MPI( IChOutWrite, ipos,jpos, Cdata, Isize, step, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara
  USE CPara
  USE SPara

  USE MPI
  USE MyMPI

  INTEGER(KIND=IKIND) ipos,jpos, Isize, step, IErr
  INTEGER my_status(MPI_STATUS_SIZE)
  COMPLEX(CKIND), DIMENSION(ISize):: &
       Cdata
  CHARACTER(IMAXCBuffer) outstring
  INTEGER ind, IChOutWrite,ITotalSize
  CHARACTER*150 SFormatstring,Sisize
  
  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"WriteDataC_MPI()"
  END IF
  
  WRITE(sisize,*) isize
  
  ! line by line

  IF ( (2*13*SIZE(Cdata)+3*6*ADD_OUT_INFO ).GT. IMAXCBuffer ) THEN
     IErr=1
     PRINT*, "WriteDataC_MPI(", my_rank, ") error ", IErr, &
          ", output buffer size IMAXCBuffer=", IMAXCBuffer, &
          " smaller than SIZEOF(data)=", SIZE(Cdata)
     GOTO 20
  ENDIF

  WRITE(SFormatString,*) &
       "(4(I7.1,1X),"//TRIM(ADJUSTL(TRIM(Sisize)))//"(1F13.10,1X), &
       "//TRIM(ADJUSTL(TRIM(sisize)))//"(1F13.10,1X),A1)"
  WRITE(outstring,FMT=SFormatString, ERR=20) my_rank,ipos,jpos,Isize,&
       REAL(Cdata,RKIND),AIMAG(Cdata), CHAR(10)
  CALL MPI_FILE_WRITE_SHARED(IChOutWrite, TRIM(outstring), LEN_TRIM(outstring), &
       MPI_CHARACTER, my_status, IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteDataC_MPI(", my_rank, ") error ", IErr, &
          " in MPI_FILE_WRITE_SHARED() for file handle ",IChOutWrite
     RETURN
  ENDIF
  RETURN

  ! error in WRITE detected
20 PRINT*,"WriteDataC_MPI(): ERR in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteDataC_MPI

SUBROUTINE WriteOutInputFile
  
  USE MyNumbers
  
  USE IConst
  USE RConst
  
  USE IPara
  USE RPara
  USE CPara
  USE SPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  IF(ISoftwareMode.LT.2) THEN
     
     OPEN(UNIT= IChInp,FILE= "Felix.inp.SimDraw_sample",&
       STATUS= 'UNKNOWN')
  
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Input file for FelixSim/Draw/Refine version :VERSION: Build :BUILD:")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# ------------------------------------")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# ------------------------------------")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# FelixSim input")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# control flags")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IWriteFLAG                = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IImageFLAG                = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IOutputFLAG               = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IBinorTextFLAG            = 0") 
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IScatterFactorMethodFLAG  = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("ICentralBeamFLAG          = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMaskFLAG                 = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IZolzFLAG                 = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IAbsorbFLAG               = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IAnisoDebyeWallerFlag     = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IBeamConvergenceFLAG      = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IPseudoCubicFLAG          = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionFLAG           = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# radius of the beam in pixels")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IPixelCount               = 64")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# beam selection criteria")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMinReflectionPool        = 600")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMinStrongBeams           = 125")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMinWeakBeams             = 20")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RBSBMax                   = 0.1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RBSPMax                   = 0.1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RConvergenceTolerance (%) = 1.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# crystal settings")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RDebyeWallerConstant      = 0.4668")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RAbsorptionPer            = 2.9")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# microscope settings")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RConvergenceAngle         = 6.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IIncidentBeamDirectionX   = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IIncidentBeamDirectionY   = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IIncidentBeamDirectionZ   = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionX              = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionY              = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionZ              = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INormalDirectionX         = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INormalDirectionY         = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INormalDirectionZ         = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RAcceleratingVoltage (kV) = 200.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Image Output Options")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RInitialThickness        = 300.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RFinalThickness          = 1300.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RDeltaThickness          = 10.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IReflectOut              = 49")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
  ELSE
     
     OPEN(UNIT= IChInp,FILE= "Felix.inp.Refine_sample",&
       STATUS= 'UNKNOWN')
  
  
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Input file for FelixSim/Draw/Refine version :VERSION: Build :BUILD:")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# ------------------------------------")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# ------------------------------------")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# FelixSim input")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# control flags")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IWriteFLAG                = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IImageFLAG                = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IOutputFLAG               = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IBinorTextFLAG            = 0") 
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IScatterFactorMethodFLAG  = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("ICentralBeamFLAG          = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMaskFLAG                 = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IZolzFLAG                 = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IAbsorbFLAG               = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IAnisoDebyeWallerFlag     = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IBeamConvergenceFLAG      = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IPseudoCubicFLAG          = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionFLAG           = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# radius of the beam in pixels")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IPixelCount               = 64")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# beam selection criteria")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMinReflectionPool        = 600")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMinStrongBeams           = 125")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMinWeakBeams             = 20")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RBSBMax                   = 0.1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RBSPMax                   = 0.1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RConvergenceTolerance (%) = 1.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# crystal settings")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RDebyeWallerConstant      = 0.4668")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RAbsorptionPer            = 2.9")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# microscope settings")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RConvergenceAngle         = 6.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IIncidentBeamDirectionX   = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IIncidentBeamDirectionY   = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IIncidentBeamDirectionZ   = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionX              = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionY              = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionZ              = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INormalDirectionX         = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INormalDirectionY         = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INormalDirectionZ         = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RAcceleratingVoltage (kV) = 200.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Image Output Options")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RInitialThickness        = 300.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RFinalThickness          = 1300.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RDeltaThickness          = 10.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IReflectOut              = 49")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# FelixRefine Input")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("#Refinement Specific Flags")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IImageOutputFLAG          = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IDevFLAG                  = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IRefineModeFLAG           = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Debye Waller Factor Iteration")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RInitialDebyeWallerFactor = 0.1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RFinalDebyeWallerFactor = 1.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RDeltaDebyeWallerFactor = 0.1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IElementsforDWFchange = {0}")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Ug Iteration")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INoofUgs                  = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RLowerBoundUgChange       = 50.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RUpperBoundUgChange       = 50.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RDeltaUgChange            = 50.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
  CLOSE(UNIT=IChInp)
END IF

END SUBROUTINE WriteOutInputFile

SUBROUTINE ReadInHKLs(IErr)

  USE MyNumbers
  USE IChannels

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara
  USE CPara
  USE SPara

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: ILine,ILength,IErr,h,k,l,ind,IPos1,IPos2
  CHARACTER*200 :: dummy1,dummy2

  
  IF((my_rank.EQ.0.AND.IWriteFLAG.GE.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ReadInHKLs(",my_rank,")"
  END IF

  OPEN(Unit = IChInp,FILE="Felix.hkl",&
       STATUS='OLD',ERR=10)

  IF((my_rank.EQ.0.AND.IWriteFLAG.GE.2).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ReadInHKLs(",my_rank,") Felix Has Detected .hkl and is entering Selected hkl mode"
  END IF
   
  ILine = 0

  IHKLSelectFLAG=1

  
  
  IF((my_rank.EQ.0.AND.IWriteFLAG.GE.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ReadInHKLs(",my_rank,") Selected HKL mode Engaged"
  END IF

  
  DO
     READ(UNIT= IChInp, END=100, FMT='(a)') dummy1
     ILine=ILine+1
  ENDDO

  100 ILength=ILine

  ALLOCATE(&
       RInputHKLs(ILength,THREEDIM),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadInHKLs(): error in memory ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE(&
       IOutputReflections(ILength),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadInHKLs(): error in memory ALLOCATE()"
     RETURN
  ENDIF

  REWIND(UNIT=IChInp)

  ILine = 0

  DO ind = 1,ILength

     ! READ HKL in as String

     READ(UNIT= IChInp,FMT='(A)' ) dummy1

     
     !Scan String for h

     IPos1 = SCAN(dummy1,'[')
     IPos2 = SCAN(dummy1,',')
     dummy2 = dummy1((IPos1+1):(IPos2-1))
     READ(dummy2,'(I20)') h

     !Scan String for k
     
     IPos1 = SCAN(dummy1((IPos2+1):),',') + IPos2
     dummy2 = dummy1((IPos2+1):(IPos1-1))
     READ(dummy2,'(I20)') k

     !Scan String for l     

     IPos2 = SCAN(dummy1((IPos1+1):),']') + IPos1
     dummy2 = dummy1((IPos1+1):(IPos2-1))
     READ(dummy2,'(I20)') l

     ! Convert to REALs

     ILine=ILine+1
     RInputHKLs(ILine,1) = REAL(h,RKIND)
     RInputHKLs(ILine,2) = REAL(k,RKIND)
     RInputHKLs(ILine,3) = REAL(l,RKIND)
     
     IF((my_rank.EQ.0.AND.IWriteFLAG.GE.3).OR.IWriteFLAG.GE.10) THEN
        
        PRINT*,RInputHKLs(ILine,:)
     END IF
  END DO

  IReflectOut = ILength

  RETURN

  10 IHKLSelectFLAG=0
  
  IF((my_rank.EQ.0.AND.IWriteFLAG.GE.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ReadInHKLs(",my_rank,") Did not find .hkl file continuing in normal mode"
  END IF
  RETURN
END SUBROUTINE ReadInHKLs
  
  
