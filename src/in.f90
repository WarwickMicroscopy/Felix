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
! $Id: inout.f90,v 1.59 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! -----------------------------------------------------------------------
!ReadInputParameters: Read the input file
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE ReadInpFile( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI
  USE WriteToScreen
  
  IMPLICIT NONE


  INTEGER(IKIND) :: &
       IErr, ILine,ind,IPos,IPos1,IPos2
  REAL(RKIND) :: &
       ROfIter
  CHARACTER*200 ::&
       SImageMode,SElements,SRefineMode,SStringFromNumber,SRefineYESNO,&
       SAtomicSites,SFormatString,SLengthofNumberString

  OPEN(UNIT= IChInp, ERR= 120, FILE= "felix.inp",&
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

!!$If the user wants to print out the normal sim messages 
!!$for debugging, they need to select a WriteFlag of above 100
  IF ((IWriteFLAG.GE.100).AND.(ISoftwareMode.EQ.2)) THEN
     IRefineSwitch = 3
  ELSE
     IRefineSwitch = 2
  END IF

  CALL Message ("ReadInpFile",IMust,IErr)
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable="IWriteFLAG",IVariable=IWriteFLAG)
  
  !IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
  !   PRINT*,"IWriteFLAG = ", IWriteFLAG
  !END IF
  
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

  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IImageFLAG",IVariable=IImageFLAG)

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IOutputFLAG
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IOutputFLAG",IVariable=IOutputFLAG)
 
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IBinorTextFLAG
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IBinorTextFLAG",IVariable=IBinorTextFLAG)

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IScatterFactorMethodFLAG
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IScatterFactorMethodFLAG",IVariable=IScatterFactorMethodFLAG)

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) ICentralBeamFLAG
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="ICentralBeamFLAG",IVariable=ICentralBeamFLAG)


  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMaskFLAG
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IMaskFLAG",IVariable=IMaskFLAG)

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IZolzFLAG
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IZolzFLAG",IVariable=IZolzFLAG)

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IAbsorbFLAG
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IAbsorbFLAG",IVariable=IAbsorbFLAG)

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IAnisoDebyeWallerFactorFlag
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IAnisoDebyeWallerFactorFlag",IVariable=IAnisoDebyeWallerFactorFlag)

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IBeamConvergenceFLAG
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IBeamConvergenceFLAG",IVariable=IBeamConvergenceFLAG)

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IPseudoCubicFLAG
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IPseudoCubicFLAG",IVariable=IPseudoCubicFLAG)

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IXDirectionFLAG
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IXDirectionFLAG",IVariable=IXDirectionFLAG)


  ! ----------------------------------------------------------------------
  ! beam details
  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IPixelCount
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IPixelCount",IVariable=IPixelCount)


  ! ----------------------------------------------------------------------
  ! beam selection criteria
  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMinReflectionPool
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IMinReflectionPool",IVariable=IMinReflectionPool)
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMinStrongBeams
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IMinStrongBeams",IVariable=IMinStrongBeams)

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMinWeakBeams
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IMinWeakBeams",IVariable=IMinWeakBeams)

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RBSBMax
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RBSBMax",RVariable=RBSBMax)

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RBSPMax
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RBSPMax",RVariable=RBSPMax)


  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RConvergenceTolerance
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RConvergenceTolerance",RVariable=RConvergenceTolerance)
  ! ----------------------------------------------------------------------
  ! crystal settings
  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RDebyeWallerConstant
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RDebyeWallerConstant",RVariable=RDebyeWallerConstant)

  RMeanSquaredDisplacement= RDebyeWallerConstant/(8*PI**2)
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RAbsorptionPercentage
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RAbsorptionPercentage",RVariable=RAbsorptionPercentage)


  ! ----------------------------------------------------------------------
  ! microscopy settings

  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RConvergenceAngle
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="ROuterConvergenceAngle",RVariable=RConvergenceAngle)

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RInnerConvergenceAngle
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RInnerConvergenceAngle",RVariable=RInnerConvergenceAngle)

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IIncidentBeamDirectionX
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IIncidentBeamDirectionX",IVariable=IIncidentBeamDirectionX)


  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IIncidentBeamDirectionY
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IIncidentBeamDirectionY",IVariable=IIncidentBeamDirectionY)


  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IIncidentBeamDirectionZ
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IIncidentBeamDirectionZ",IVariable=IIncidentBeamDirectionZ)


  RZDirC(1)= REAL(IIncidentBeamDirectionX,RKIND)
  RZDirC(2)= REAL(IIncidentBeamDirectionY,RKIND)
  RZDirC(3)= REAL(IIncidentBeamDirectionZ,RKIND)
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IXDirectionX
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IXDirectionX",IVariable=IXDirectionX)


  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IXDirectionY
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IXDirectionY",IVariable=IXDirectionY)
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IXDirectionZ
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IXDirectionZ",IVariable=IXDirectionZ)

  RXDirC(1)= IXDirectionX
  RXDirC(2)= IXDirectionY
  RXDirC(3)= IXDirectionZ  

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) INormalDirectionX
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="INormalDirectionX",IVariable=INormalDirectionX)

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) INormalDirectionY
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="INormalDirectionY",IVariable=INormalDirectionY)
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) INormalDirectionZ
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="INormalDirectionZ",IVariable=INormalDirectionZ)

  RNormDirC(1)= INormalDirectionX
  RNormDirC(2)= INormalDirectionY
  RNormDirC(3)= INormalDirectionZ
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RAcceleratingVoltage
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RAcceleratingVolatage",RVariable=RAcceleratingVoltage)

  ! ----------------------------------------------------------------------
  ! Title Space
  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ! ----------------------------------------------------------------------
  ! Image Output Options

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RInitialThickness
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RInitialThickness",RVariable=RInitialThickness)

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RFinalThickness
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RFinalThickness",RVariable=RFinalThickness)

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RDeltaThickness
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RDeltaThickness",RVariable=RDeltaThickness)
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IReflectOut
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IReflectOut",IVariable=IReflectOut)


  IF(ISoftwareMode.EQ.2) THEN

     !-----------------------------------------------------------------------
     ! felixrefine Input
     
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     
     !-----------------------------------------------------------------------
     ! Refinement Specific Flags
     
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     
     ILine= ILine+1

     READ(IChInp,FMT='(A)',ERR=20,END=30) SRefineMode

     IRefineModeSelectionArray = 0

     DO ind = 1,7
        WRITE(SStringFromNumber,'(I1)') ind-1
        IF(SCAN(SRefineMode,TRIM(ADJUSTL(SStringFromNumber))).NE.0) THEN
           IRefineModeSelectionArray(ind) = 1
        END IF
     END DO
     
     IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        IF(IWriteFLAG.GE.4) THEN
           DO ind = 1,7
              SRefineYESNO = 'NO'
              SELECT CASE (ind)
              CASE(1)
                 IF(IRefineModeSelectionArray(ind).EQ.1) SRefineYESNO = 'YES'
                 PRINT*,"Refine Structure Factors ",SRefineYESNO
              CASE(2)
                 IF(IRefineModeSelectionArray(ind).EQ.1) SRefineYESNO = 'YES'
                 PRINT*,"Refine Atomic Coordinates ",SRefineYESNO
              CASE(3)
                 IF(IRefineModeSelectionArray(ind).EQ.1) SRefineYESNO = 'YES'
                 PRINT*,"Refine Atomic Site Occupancies ",SRefineYESNO
              CASE(4)
                 IF(IRefineModeSelectionArray(ind).EQ.1) SRefineYESNO = 'YES'
                 PRINT*,"Refine Isotropic Debye Waller Factors ",SRefineYESNO
              CASE(5)
                 IF(IRefineModeSelectionArray(ind).EQ.1) SRefineYESNO = 'YES'
                 PRINT*,"Refine Anisotropic Debye Waller Factors ",SRefineYESNO
              CASE(6)
                 IF(IRefineModeSelectionArray(ind).EQ.1) SRefineYESNO = 'YES'
                 PRINT*,"Refine Lattice Lengths ",SRefineYESNO
              CASE(7)
                 IF(IRefineModeSelectionArray(ind).EQ.1) SRefineYESNO = 'YES'
                 PRINT*,"Refine Lattice Angles ",SRefineYESNO
              END SELECT
           END DO
        ELSE
           PRINT*,"IRefineModeSelectionArray = ",IRefineModeSelectionArray
        END IF
     END IF    
     
     ILine= ILine+1
     READ(IChInp,10,ERR=20,END=30) IWeightingFLAG
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IWeightingFLAG",IVariable=IWeightingFLAG)
     
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     
     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RInitialDebyeWallerFactor
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RInitialDebyeWallerFactor",RVariable =  RInitialDebyeWallerFactor)

     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RFinalDebyeWallerFactor
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RFinalDebyeWallerFactor",RVariable =  RFinalDebyeWallerFactor)
     
     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RDeltaDebyeWallerFactor
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RDeltaDebyeWallerFactor",RVariable =  RDeltaDebyeWallerFactor)

     ILine= ILine+1
     READ(IChInp,FMT='(A)',ERR=20,END=30) SElements
     IPos1 = SCAN(SElements,'{')
     IPos2 = SCAN(SElements,'}')
     IPos = SCAN(SElements,'0')
     IElements = 1

     !not sure what type of message this is?
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
           PRINT*,"ReadInpFile(",my_rank,") ERROR IN ALLOCATE OF IElementList"
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
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     
     ILine= ILine+1
     READ(IChInp,10,ERR=20,END=30) INoofUgs
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="INoofUgs",IVariable = INoofUgs)
     
     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RLowerBoundUgChange
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RLowerBoundUgChange",RVariable = RLowerBoundUgChange)
     
     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RUpperBoundUgChange
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RUpperBoundUgChange",RVariable = RUpperBoundUgChange)

     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RDeltaUgChange
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RDeltaUgChange",RVariable = RDeltaUgChange)

     
     !-----------------------------------------------------------------------
     ! Iterative Structural input

     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     
     READ(IChInp,FMT='(A)',ERR=20,END=30) SAtomicSites

     CALL DetermineRefineableAtomicSites(SAtomicSites,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadInpFile(): error in DetermineRefineableAtomicSites()"
        RETURN
     ENDIF

     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     
     ILine= ILine+1
     READ(IChInp,10,ERR=20,END=30) IPrint
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IPrint",IVariable = IPrint)
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     
     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RSimplexLengthScale
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RSimplexLengthScale",RVariable = RSimplexLengthScale)
     RSimplexLengthScale = RSimplexLengthScale/100.0

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
     PRINT*,"ReadInpFile(): ERR in OPEN"
     PRINT*,""
     CALL WriteOutInputFile (IErr)
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
     CALL WriteOutInputFile (IErr)
  END IF
  GOTO 1000
  
  !	EOF in READ occured prematurely
30 IF(my_rank.EQ.0.OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Input(): EOF in READ at line", ILine
     CALL WriteOutInputFile (IErr)
  END IF
  
  ! dump the input help
  
1000 IF((my_rank.EQ.0.AND.ISoftwareMode.LT.2).OR.(IWriteFLAG.GE.10.AND.ISoftwareMode.LT.2)) THEN
     PRINT*,"# Input file for felixsim/draw version :VERSION: Build :BUILD:"
     PRINT*,"# ------------------------------------"
     PRINT*,""
     PRINT*,"# ------------------------------------"
     PRINT*,"# felixsim input"
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
     PRINT*,"ROuterConvergenceAngle    = 6.0"
     PRINT*,"RInnerConvergenceAngle    = 3.0"
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
     
     PRINT*,"A Sample Input File Has been Written For you as felix.inp.simdraw_sample"
     PRINT*,"It must be renamed to felix.inp before use"
  ELSE
     PRINT*,"# Input file for felixrefine version :VERSION: Build :BUILD:"
     PRINT*,"# ------------------------------------"
     PRINT*,""
     PRINT*,"# ------------------------------------"
     PRINT*,"# felixsim input"
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
     PRINT*,"ROuterConvergenceAngle    = 6.0"
     PRINT*,"RInnerConvergenceAngle    = 3.0"
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
     PRINT*,"# felixrefine Input"
     PRINT*,""
     PRINT*,"#Refinement Specific Flags"
     PRINT*,"IRefineModeFLAG           = 0"
     PRINT*,"IWeightingFLAG           = 0"
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
     PRINT*,"# Structural Refinement"
     PRINT*,""
     PRINT*,"IAtomicSites              = (1,2,6)"
     PRINT*,""
     PRINT*,"# Refinement Output"
     PRINT*,""
     PRINT*,"IPrint                    = 10"
     PRINT*,""
     PRINT*,"IAtomicSites              = (1,2,3,4,5,6,7)"
     PRINT*,""
     PRINT*,"# Refinement Output"
     PRINT*,""
     PRINT*,"IPrint                    = 20"
     PRINT*,""
     PRINT*,"# Simplex Initialisation"
     PRINT*,""
     PRINT*,"RSimplexLengthScale       = 5.0"
     PRINT*,""     

     PRINT*,"A Sample Input File Has been Written For you as felix.inp.refine_sample"
     PRINT*,"It must be renamed to felix.inp before use"
  END IF
  IErr= 1
  RETURN
!changed name from input to readinputparameters  
END SUBROUTINE ReadInpFile
   

! -----------------------------------------------------------------------
!
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE ReadScaFile( IErr )

  USE MyNumbers
  USE WriteToScreen
  
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

     CALL Message ("ReadScaFile",IMust,IErr)

  ILine= 0
  
!!$  OPEN(UNIT= IChInp, ERR= 120, FILE= "FelixDoyle.sca",&
!!$       STATUS= 'OLD')
  OPEN(UNIT= IChInp, ERR= 120, FILE= "felix.sca",&
       STATUS= 'OLD')
!!$  
  DO
     READ(UNIT= IChInp, END=100, ERR=20, FMT='(a)') dummy
     ILine=ILine+1
  ENDDO
100 ILength=ILine
  
  IF(IWriteFLAG.GE.10) THEN
     
     PRINT*,"ReadScaFile(): ILength=", ILength
     
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
        PRINT*,"ReadScaFile(): error in memory ALLOCATE()"
        RETURN
     ENDIF

     ALLOCATE(RScattFactors(ILength,12), STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadScaFile(): error in memory ALLOCATE()"
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
        PRINT*,"ReadScaFile(): error in memory ALLOCATE()"
        RETURN
     ENDIF

     ALLOCATE(RScattFactors(ILength,8), STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadScaFile(): error in memory ALLOCATE()"
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
        PRINT*,"ReadScaFile(): error in memory ALLOCATE()"
        RETURN
     ENDIF

     ALLOCATE(RScattFactors(ILength,8), STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadScaFile(): error in memory ALLOCATE()"
        RETURN
     ENDIF
     
     DO ILine=1,ILength,1
        READ(UNIT= IChInp, FMT='(29(E15.11,1X))') &
             rdummy(:), RScattFactors(ILine,:)     
        PRINT*,RScattFactors(ILine,:)
     ENDDO

  END SELECT

  DEALLOCATE( &
       rdummy,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadScaFile(): error in memory DEALLOCATE()"
     RETURN
  ENDIF  
  
  CLOSE(IChInp)
  
  RETURN
  
  !	error in OPEN detected
120 PRINT*,"ReadScaFile(): ERR in OPEN"
  GOTO 1000
  
  !	error in CLOSE detected
130 PRINT*,"ReadScaFile(): ERR in CLOSE"
  GOTO 1000
  
  !	error in READ detected
20 PRINT*,"ReadScaFile(): ERR in READ at line", ILine
  GOTO 1000
  
  !	EOF in READ occured prematurely
30 PRINT*,"ReadScaFile(): EOF in READ at line", ILine
  
!bug need to change layout
1000 &
  PRINT*,"ReadScaFile expects a file such as"
  PRINT*,"--------------------------------------------------------------------"
  PRINT*,"22      0.737291729     9.96300042      0.355417565     9.96300042 ", &
       "0.496982599     0.072659586     0.016335034     0.360199179     1.42171303 ", &
       "     0.073173566     0.957232308      15.8512114"
  PRINT*,"23      0       0       0       0       0       0       0       0  ", &
       "     0       0       0       0"

  IErr= 1
  RETURN
END SUBROUTINE ReadScaFile
  
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
          &"//TRIM(ADJUSTL(TRIM(EIGENFORMAT)))//"(1F13.10,1X), &
          &"//TRIM(ADJUSTL(TRIM(EIGENFORMAT)))//"(1F13.10,1X))", &
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
  
!Write out the sample input file, when none provided
SUBROUTINE WriteOutInputFile (IErr)
  
  USE WriteToScreen
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

  INTEGER(IKIND):: IErr

  IF(ISoftwareMode.LT.2) THEN

     CALL Message("WriteOutInputFile",IMust,IErr)
     
     OPEN(UNIT= IChInp,FILE= "felix.inp.simdraw_sample",&
       STATUS= 'UNKNOWN')
  
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Input file for felixsim/draw/refine version :VERSION: Build :BUILD:")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# ------------------------------------")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# ------------------------------------")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# felixsim input")
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
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("ROuterConvergenceAngle    = 6.0") 
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RInnerConvergenceAngle    = 0.0")
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
     
     OPEN(UNIT= IChInp,FILE= "felix.inp.refine_sample",&
       STATUS= 'UNKNOWN')
  
  
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Input file for felixsim/draw/refine version :VERSION: Build :BUILD:")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# ------------------------------------")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# ------------------------------------")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# felixsim input")
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
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("ROuterConvergenceAngle    = 6.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RInnerConvergenceAngle    = 0.0")
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
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# felixrefine Input")
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
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Structural Refinement")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IAtomicSites              = (1,2,6)")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Refinement Output")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IPrint                    = 10")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Simplex Initialisation")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RSimplexLengthScale       = 5.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
  CLOSE(UNIT=IChInp)
END IF

END SUBROUTINE WriteOutInputFile

SUBROUTINE ReadHklFile(IErr)

  USE WriteToScreen
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

  CALL Message ("ReadHklFile",IMust,IErr)  
 ! IF((my_rank.EQ.0.AND.IWriteFLAG.GE.0).OR.IWriteFLAG.GE.10) THEN
 !    PRINT*,"ReadHklFile(",my_rank,")"
 ! END IF

  OPEN(Unit = IChInp,FILE="felix.hkl",&
       STATUS='OLD',ERR=10)

  CALL Message ("ReadHklFile",IInfo,IErr, &
       MessageString =" Felix has detected .hkl and is entering selected hkl mode") 
  !IF((my_rank.EQ.0.AND.IWriteFLAG.GE.2).OR.IWriteFLAG.GE.10) THEN
  !   PRINT*,"ReadHklFile(",my_rank,") felix Has Detected .hkl and is entering Selected hkl mode"
  !END IF
   
  ILine = 0

  IHKLSelectFLAG=1
  
  CALL Message ("ReadHklFile",IInfo,IErr,MessageString = " Selected hkl mode engaged") 
 ! IF((my_rank.EQ.0.AND.IWriteFLAG.GE.0).OR.IWriteFLAG.GE.10) THEN
 !    PRINT*,"ReadHklFile(",my_rank,") Selected hkl Mode Engaged"
 ! END IF

  
  DO
     READ(UNIT= IChInp, END=100, FMT='(a)') dummy1
     ILine=ILine+1
  ENDDO

  100 ILength=ILine

  ALLOCATE(&
       RInputHKLs(ILength,THREEDIM),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadHklFile(): error in memory ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE(&
       IOutputReflections(ILength),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadHklFile(): error in memory ALLOCATE()"
     RETURN
  ENDIF

  REWIND(UNIT=IChInp)

  ILine = 0

CALL Message ("ReadHklFile",IMoreInfo,IErr,MessageString = " Printing selected hkl's") 

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
     
     !Why Zero?
     IF((my_rank.EQ.0.AND.IWriteFLAG.GE.3).OR.IWriteFLAG.GE.10) THEN
        
        PRINT*,RInputHKLs(ILine,:)
     END IF
  END DO

  IReflectOut = ILength

  RETURN

  10 IHKLSelectFLAG=0
  
  CALL Message ("ReadHklFile",IInfo,IErr,MessageString = "Did not find .hkl file continuing in normal mode") 

  RETURN
END SUBROUTINE ReadHklFile
  
  
SUBROUTINE DetermineRefineableAtomicSites(SAtomicSites,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE  
  
  INTEGER(IKIND) :: &
       IPos,IPos1,IPos2,IErr,ind
  CHARACTER*200 :: &
       SAtomicSites,SFormatString,SLengthofNumberString

  IPos1 = SCAN(SAtomicSites,'(')
  IPos2 = SCAN(SAtomicSites,')')
  IF(((IPos2-IPos1).EQ.1).OR.(IPos1.EQ.0).OR.(IPos2.EQ.0)) THEN
     IF(IRefineModeSelectionArray(2).EQ.1) THEN
        PRINT*,"You Have Not Specfied Atomic Sites to Refine" 
        IErr = 1
        RETURN
     END IF
  END IF

  IF ((IPos2-IPos1).GT.1.AND.SCAN(SAtomicSites,',').EQ.0) THEN
     ALLOCATE(&
          IAtomicSitesToRefine(1),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadInpFile(): error in memory ALLOCATE()"
        RETURN
     ENDIF

     
     PRINT*,SIZE(IAtomicSitesToRefine)
     WRITE(SLengthofNumberString,*) LEN(SAtomicSites((IPos1+1):(IPos2-1))) 
     WRITE(SFormatString,*) "(I"//TRIM(ADJUSTL(SLengthofNumberString))//")"
     READ(SAtomicSites((IPos1+1):(IPos2-1)),FMT=SFormatString) IAtomicSitesToRefine(1)
  ELSE
     IPos = 1
     DO 
        IF(SCAN(SAtomicSites(IPos1:IPos2),',').NE.0) THEN
           IPos1 = IPos1 + LEN(SAtomicSites(IPos1:(IPos1+SCAN(SAtomicSites(IPos1:IPos2),','))))
           IPos = IPos+1
        END IF
        IF (IPos2-IPos1.LE.1) EXIT
     END DO
     
     ALLOCATE(&
          IAtomicSitesToRefine(IPos),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadInpFile(): error in memory ALLOCATE()"
        RETURN
     ENDIF
     
     IPos1 = SCAN(SAtomicSites,'(')
     DO ind = 1,SIZE(IAtomicSitesToRefine,DIM=1)
        IF(SCAN(SAtomicSites((IPos1+1):IPos2),',').NE.0) THEN
           IPos = SCAN(SAtomicSites((IPos1+1):IPos2),',')-1
           WRITE(SLengthofNumberString,*) LEN(SAtomicSites((IPos1+1):(IPos1+IPos))) 
           WRITE(SFormatString,*) "(I"//TRIM(ADJUSTL(SLengthofNumberString))//")"
           READ(SAtomicSites((IPos1+1):(IPos1+IPos)),FMT=SFormatString) IAtomicSitesToRefine(ind)
           IPos1 = IPos1 + IPos + 1 
        ELSE
           WRITE(SLengthofNumberString,*) LEN(SAtomicSites((IPos1+1):(IPos2-1))) 
           WRITE(SFormatString,*) "(I"//TRIM(ADJUSTL(SLengthofNumberString))//")"
           READ(SAtomicSites((IPos1+1):(IPos2-1)),FMT=SFormatString) IAtomicSitesToRefine(ind)
        END IF
     END DO
     
  END IF
  
END SUBROUTINE DetermineRefineableAtomicSites
