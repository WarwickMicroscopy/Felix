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


  INTEGER(IKIND) :: IErr, ILine,ind,IPos,IPos1,IPos2
  REAL(RKIND) :: ROfIter
  CHARACTER*200 :: SImageMode,SElements,SRefineMode,SStringFromNumber,SRefineYESNO,&
       SAtomicSites,SFormatString,SLengthofNumberString
  CHARACTER*200 :: SDirectionX,SIncidentBeamDirection,SNormalDirectionX
  

  OPEN(UNIT= IChInp, ERR= 120, FILE= "felix.inp",STATUS= 'OLD')
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

  !NB Mean Square Displacement= RDebyeWallerConstant/(8*PI**2)
  
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


  ! RZDirC,RXDirC,RNormDirC vectors are reciprocal lattice vectors that define the beam direction, x-axis of the
  ! diffraction pattern and the surface normal respectively
  ILine= ILine+1
  READ(IChInp,FMT='(27X,A)',END=30) SIncidentBeamDirection
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IIncidentBeamDirection", &
       MessageString=ADJUSTL(TRIM(SIncidentBeamDirection)))
  CALL ThreeDimVectorReadIn(SIncidentBeamDirection,'[',']',RZDirC)
  CALL Message ("ReadInpFile",IInfo+IDebug,IErr,MessageVariable ="RZDirC",RVector=RZDirC)

  ILine= ILine+1
  READ(IChInp,FMT='(27X,A)',END=30) SDirectionX

  !Read in the user specified X-direction, if 'A' felix automatically selects the closest g-vector
  IF (INDEX(SDirectionX,'A').NE.0) THEN
     IXDirectionFLAG=0
  ELSE
     IXDirectionFLAG=1
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IDirection", &
          MessageString=ADJUSTL(TRIM(SDirectionX)))
     CALL ThreeDimVectorReadIn(SDirectionX,'[',']',RXDirC)
     CALL Message ("ReadInpFile",IInfo+IDebug,IErr,MessageVariable ="RXDirC",RVector=RXDirC)
  END IF

  ILine= ILine+1
  READ(IChInp,FMT='(27X,A)',ERR=20,END=30) SNormalDirectionX
     CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="INormalDirection", &
          MessageString=ADJUSTL(TRIM(SNormalDirectionX)))
     CALL ThreeDimVectorReadIn(SNormalDirectionX,'[',']',RNormDirC)
     CALL Message ("ReadInpFile",IInfo+IDebug,IErr,MessageVariable ="RNormDirC",RVector=RNormDirC)
  END IF

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
  READ(IChInp,10,ERR=20,END=30) INoOfLacbedPatterns
  CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="INoOfLacbedPatterns",IVariable=INoOfLacbedPatterns)


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
     IRefineMode = 0

     DO ind = 1,IRefinementVariableTypes
        IF(SCAN(TRIM(ADJUSTL(SRefineMode)),TRIM(ADJUSTL(CAlphabet(ind)))).NE.0) THEN
           IRefineMode(ind) = 1
        END IF
     END DO
     
     IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
       IF(IRefineMode(1).EQ.1) PRINT*, "A:Refining Structure Factors by Simplex"
       IF(IRefineMode(2).EQ.1) PRINT*, "B:Refining Atomic Coordinates"
       IF(IRefineMode(3).EQ.1) PRINT*, "C:Refining Occupancies "
       IF(IRefineMode(4).EQ.1) PRINT*, "D:Refining Isotropic Debye Waller Factors"
       IF(IRefineMode(5).EQ.1) PRINT*, "E:Refining Anisotropic Debye Waller Factors "
       IF(IRefineMode(6).EQ.1) PRINT*, "F:Refining Lattice Lengths "
       IF(IRefineMode(7).EQ.1) PRINT*, "G:Refining Lattice Angles "
       IF(IRefineMode(8).EQ.1) PRINT*, "H:Refining Convergence Angle "
       IF(IRefineMode(9).EQ.1) PRINT*, "I:Refining Absorption"
       IF(IRefineMode(10).EQ.1) PRINT*,"J:Refining Accelerating Voltage "
       IF(IRefineMode(11).EQ.1) PRINT*,"K:Refining Scale Factor "
       IF(IRefineMode(12).EQ.1) PRINT*,"L:Refining Structure Factors by bisection"
     END IF
    
     !Check if user has requested Ug refinement and anything else which isnt possible
        
     IF((IRefineMode(1).EQ.1 .OR. IRefineMode(12).EQ.1).AND.SUM(IRefineMode).GT.1) THEN         
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
  IErr= 1
  RETURN
!changed name from input to readinputparameters  
END SUBROUTINE ReadInpFile

! -----------------------------------------------------------------------

SUBROUTINE ReadScaFile( IErr )
!completely pointless subroutine that just calls another subroutine
  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE


  INTEGER IErr, ILine

     CALL Message("ReadScaFile",IMust,IErr)
     CALL Message("ReadScaFile",IInfo,IErr,MessageString="Reading in Scattering Factors")
     CALL ScatteringFactors(IScatterFactorMethodFLAG,IErr)
  
END SUBROUTINE ReadScaFile
  
! -----------------------------------------------------------------------

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
  CHARACTER*200 :: dummy1,dummy2,SHKLString

  CALL Message ("ReadHklFile",IMust,IErr)  

  OPEN(Unit = IChInp,FILE="felix.hkl",STATUS='OLD',ERR=10)

  CALL Message ("ReadHklFile",IInfo,IErr,MessageString ="Using hkl list") 
   
  ILine = 0

  IHKLSelectFLAG=1
  
  DO
     READ(UNIT= IChInp, END=100, FMT='(a)') dummy1
     ILine=ILine+1
  ENDDO

  100 ILength=ILine

  ALLOCATE(RInputHKLs(ILength,ITHREE),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadHklFile(",my_rank,")error allocating RInputHKLs"
     RETURN
  ENDIF
  ALLOCATE(IOutputReflections(ILength),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadHklFile(",my_rank,")error allocating IOutputReflections"
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

     IF((my_rank.EQ.0.AND.IWriteFLAG.GE.3).OR.IWriteFLAG.GE.10) THEN
        WRITE(SHKLString,'(3(I4.1,1X))') NINT(RInputHKLs(ILine,:))
        CALL Message ("ReadHklFile",IInfo,IErr,MessageVariable = "Input HKL is",MessageString = SHKLString) 
     END IF
  END DO

  INoOfLacbedPatterns = ILength

  RETURN

  10 IHKLSelectFLAG=0
  
  CALL Message ("ReadHklFile",IInfo,IErr,MessageString = "Did not find .hkl file continuing in normal mode") 

  RETURN
END SUBROUTINE ReadHklFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
SUBROUTINE DetermineRefineableAtomicSites(SAtomicSites,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE  
  
  INTEGER(IKIND) :: IPos,IPos1,IPos2,IErr,ind
  CHARACTER*200 :: SAtomicSites,SFormatString,SLengthofNumberString

  IPos1 = SCAN(SAtomicSites,'(')
  IPos2 = SCAN(SAtomicSites,')')
  IF(((IPos2-IPos1).EQ.1).OR.(IPos1.EQ.0).OR.(IPos2.EQ.0)) THEN
     IF(IRefineMode(2).EQ.1) THEN
        PRINT*,"You Have Not Specfied Atomic Sites to Refine" 
        IErr = 1
        RETURN
     END IF
  END IF

  IF ((IPos2-IPos1).GT.1.AND.SCAN(SAtomicSites,',').EQ.0) THEN
     ALLOCATE(IAtomicSitesToRefine(1),STAT=IErr)
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
     
     ALLOCATE(IAtomicSitesToRefine(IPos),STAT=IErr)
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
!XX PRINT*, "Refining atoms",IAtomicSitesToRefine!XX     
  
END SUBROUTINE DetermineRefineableAtomicSites

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
  
  INTEGER(IKIND) :: ind,jnd,IErr
  INTEGER(IKIND) :: INegError = 0
  CHARACTER*34 :: filename
  CHARACTER*200 :: SPrintString

  DO ind = 1,INoOfLacbedPatterns
     
     WRITE(filename,"(A6,I3.3,A4)") "felix.",ind,".img"
     OPEN(UNIT= IChInImage, ERR= 10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),FORM='UNFORMATTED',&
       ACCESS='DIRECT',IOSTAT=Ierr,RECL=2*IPixelCount*8)
	   
    ALLOCATE(RImageIn(2*IPixelCount,2*IPixelCount), STAT=IErr)  
    IF( IErr.NE.0 ) THEN
      PRINT*,"ReadExperimentalImages(",my_rank,")error allocating RImageIn"
      RETURN
    ENDIF
    
     DO jnd=1,2*IPixelCount
        READ(IChInImage,rec=jnd,ERR=10) RImageIn(jnd,:)
     END DO
     
     IF(MINVAL(RImageIn).LT.ZERO.AND.(my_rank.EQ.0)) THEN
        PRINT*,"Warning! There are negative values in your experimental images"
        INegError = INegError + 1
     END IF

     RImageExpi(:,:,ind) = RImageIn
     DEALLOCATE(RImageIn, STAT=IErr)  
     
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

  IF (my_rank.EQ.0) THEN
    IF (IErr.EQ.0) THEN
     WRITE(SPrintString,FMT='(I3,A40)') INoOfLacbedPatterns," experimental images successfully loaded"
     PRINT*,TRIM(ADJUSTL(SPrintString))
    END IF
  END IF

  RETURN

10 IErr=1
    PRINT*,"ReadExperimentalImages(", my_rank, ") error in READ() at record=", &
         jnd, " of file ", TRIM(filename), " with ind=", ind
    RETURN

END SUBROUTINE ReadExperimentalImages

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
