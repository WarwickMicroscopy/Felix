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
! $Id: out.f90,v 1.59 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------------------------------------------------
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
  USE InputModules
  
  IMPLICIT NONE


  INTEGER(IKIND) :: &
       IErr, ILine,ind,IPos,IPos1,IPos2
  REAL(RKIND) :: &
       ROfIter
  CHARACTER*200 ::&
       SImageMode,SElements,SRefineMode,SStringFromNumber,SRefineYESNO,&
       SAtomicSites,SFormatString,SLengthofNumberString
  CHARACTER*200 :: &
       SDirectionX,SIncidentBeamDirection,SNormalDirectionX
  

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
  READ(IChInp,10,ERR=20,END=30) IScatterFactorMethodFLAG
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IScatterFactorMethodFLAG",IVariable=IScatterFactorMethodFLAG)

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
  READ(IChInp,FMT='(27X,A)',END=30) SIncidentBeamDirection
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IIncidentBeamDirection", &
       MessageString=ADJUSTL(TRIM(SIncidentBeamDirection)))
  CALL ThreeDimVectorReadIn(SIncidentBeamDirection,'[',']',RZDirC)
  CALL Message ("ReadInpFile",IInfo+IDebug,IErr,MessageVariable ="RZDirC",RVector=RZDirC)

  ILine= ILine+1
  READ(IChInp,FMT='(27X,A)',END=30) SDirectionX
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IDirection", &
       MessageString=ADJUSTL(TRIM(SDirectionX)))
  CALL ThreeDimVectorReadIn(SDirectionX,'[',']',RXDirC)
  CALL Message ("ReadInpFile",IInfo+IDebug,IErr,MessageVariable ="RXDirC",RVector=RXDirC)

  ILine= ILine+1
  READ(IChInp,FMT='(27X,A)',ERR=20,END=30) SNormalDirectionX
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="INormalDirection", &
       MessageString=ADJUSTL(TRIM(SNormalDirectionX)))
  CALL ThreeDimVectorReadIn(SNormalDirectionX,'[',']',RNormDirC)
  CALL Message ("ReadInpFile",IInfo+IDebug,IErr,MessageVariable ="RNormDirC",RVector=RNormDirC)

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RAcceleratingVoltage
  CALL Message ("ReadInpFile",IInfo+IDebug,IErr,MessageVariable ="RAcceleratingVolatage",RVariable=RAcceleratingVoltage)

  ILine=ILine+1
  READ(IChInp,15,ERR=20,END=30) RAcceptanceAngle
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RAcceptanceAngle",RVariable=RAcceptanceAngle)

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
     SRefineMode = SRefineMode((SCAN(SRefineMode,"=")+1):)
     IRefineModeSelectionArray = 0

     DO ind = 1,IRefinementVariableTypes
!!$        PRINT*,CAlphabet(ind),SRefineMode,SCAN(TRIM(ADJUSTL(SRefineMode)),TRIM(ADJUSTL(CAlphabet(ind))))
        IF(SCAN(TRIM(ADJUSTL(SRefineMode)),TRIM(ADJUSTL(CAlphabet(ind)))).NE.0) THEN
           IRefineModeSelectionArray(ind) = 1
        END IF
     END DO
     
     IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        IF(IWriteFLAG.GE.4) THEN
           DO ind = 1,IRefinementVariableTypes
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
              CASE(8)
                 IF(IRefineModeSelectionArray(ind).EQ.1) SRefineYESNO = 'YES'
                 PRINT*,"Refine Convergence Angle ",SRefineYESNO
              CASE(9)
                 IF(IRefineModeSelectionArray(ind).EQ.1) SRefineYESNO = 'YES'
                 PRINT*,"Refine Total Absorption ",SRefineYESNO
              CASE(10)
                 IF(IRefineModeSelectionArray(ind).EQ.1) SRefineYESNO = 'YES'
                 PRINT*,"Refine Accelerating Voltage ",SRefineYESNO
              CASE(11)
                 IF(IRefineModeSelectionArray(ind).EQ.1) SRefineYESNO = 'YES'
                 PRINT*,"Refine Residual Sum of Squares Scaling Factor ",SRefineYESNO
              END SELECT
           END DO
        ELSE
           PRINT*,"IRefineModeSelectionArray = ",IRefineModeSelectionArray
        END IF
     END IF
    
     !Check if user has requested Ug refinement and anything else which isnt possible
        
     IF(IRefineModeSelectionArray(1).EQ.1.AND.SUM(IRefineModeSelectionArray).GT.1) THEN         
!!$        CALL Message ("ReadInpFile",IMust,IErr,MessageVariable ="Structure factors must be refined seperately")
        IF(my_rank.EQ.0) THEN
           PRINT*,"Structure factors must be refined seperately"
        END IF
       IErr = 1
        RETURN
     END IF


     ILine= ILine+1
     READ(IChInp,10,ERR=20,END=30) IWeightingFLAG
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IWeightingFLAG",IVariable=IWeightingFLAG)
     
     ILine= ILine+1
     READ(IChInp,10,ERR=20,END=30) IContinueFLAG
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IContinueFLAG",IVariable=IContinueFLAG)
     
     ILine= ILine+1
     READ(IChInp,10,ERR=20,END=30) ICorrelationFLAG
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="ICorrelationFLAG",IVariable=ICorrelationFLAG)
     
     ILine= ILine+1
     READ(IChInp,10,ERR=20,END=30) IImageProcessingFLAG
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IImageProcessingFLAG",IVariable=IImageProcessingFLAG)
     
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
     ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

     ILine= ILine+1
     READ(IChInp,10,ERR=20,END=30) INoofUgs
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="INoofUgs",IVariable = INoofUgs)
     
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

     ILine= ILine+1
     READ(IChInp,15,ERR=20,END=30) RExitCriteria
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="RExitCriteria",RVariable = RExitCriteria)
  END IF

10 FORMAT(27X,I15.1)
  ! 10	FORMAT("IMAXIteration= ",I15.1)
15 FORMAT(27X,F18.9)
  ! 15     FORMAT("IMAXIteration= ",F18.9)
  
  CLOSE(IChInp, ERR=130)
  
  ! check the parameters for validity
   
  RETURN
  
  !	error in OPEN detected
120 IF(my_rank.EQ.0) THEN
     PRINT*,"ReadInpFile(): ERR in OPEN"
     PRINT*,""
     CALL WriteOutInputFile (IErr)
  END IF
  
  !	error in CLOSE detected
130 IF(my_rank.EQ.0) THEN
     PRINT*,"Input(): ERR in CLOSE"
  END IF
  
  !	error in READ detected
20 IF(my_rank.EQ.0) THEN
     PRINT*,"Input(): ERR in READ at line", ILine
     CALL WriteOutInputFile (IErr)
  END IF
  
  !	EOF in READ occured prematurely
30 IF(my_rank.EQ.0) THEN
     PRINT*,"Input(): EOF in READ at line", ILine
     CALL WriteOutInputFile (IErr)
  END IF
  
  ! dump the input help

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

  INTEGER(IKIND) :: &
       ILine,ILength,IErr,h,k,l,ind,IPos1,IPos2
  CHARACTER*200 :: &
       dummy1,dummy2,SHKLString

  CALL Message ("ReadHklFile",IMust,IErr)  

  OPEN(Unit = IChInp,FILE="felix.hkl",&
       STATUS='OLD',ERR=10)

  CALL Message ("ReadHklFile",IInfo,IErr, &
       MessageString =" Felix has detected .hkl and is entering selected hkl mode") 
   
  ILine = 0

  IHKLSelectFLAG=1
  
  CALL Message ("ReadHklFile",IInfo,IErr,MessageString = " Selected hkl mode engaged") 
  
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
!!$     IF((my_rank.EQ.0.AND.IWriteFLAG.GE.3).OR.IWriteFLAG.GE.10) THEN
     WRITE(SHKLString,'(3(I4.1,1X))') NINT(RInputHKLs(ILine,:))
     CALL Message ("ReadHklFile",IInfo,IErr,MessageVariable = "Input HKL is",MessageString = SHKLString) 

!!$     END IF
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
  INTEGER(IKIND) :: &
       INegError = 0
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
     
     IF(MINVAL(RImageIn).LT.ZERO.AND.(my_rank.EQ.0)) THEN
        PRINT*,"Warning! There are negative values in your experimental images"
        INegError = INegError + 1
     END IF

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

  IF (INegError.NE.0) THEN
     IErr = 1
     PRINT*,"No. of Images with Negative Values",INegError
  END IF


END SUBROUTINE ReadExperimentalImages

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
