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
  ! There are six introductory comment lines which are ignored
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)') 
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ! -----IWriteFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IWriteFLAG
  ! -----IImageFLAG-----------------------------------------------------------------
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
  ! -----IScatterFactorMethodFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IScatterFactorMethodFLAG
  IF(my_rank.EQ.0.AND.IWriteFLAG.EQ.7) PRINT*,"IScatterFactorMethodFLAG=",IScatterFactorMethodFLAG
  ! -----IMaskFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IMaskFLAG
  !!!CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IMaskFLAG",IVariable=IMaskFLAG)
  IF(my_rank.EQ.0.AND.IWriteFLAG.EQ.7) PRINT*,"IMaskFLAG=",IMaskFLAG
  ! -----IHolzFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IHolzFLAG
  !!!CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IHolzFLAG",IVariable=IHolzFLAG)
  IF(my_rank.EQ.0.AND.IWriteFLAG.EQ.7) PRINT*,"IHolzFLAG=",IHolzFLAG
  ! -----IAbsorbFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IAbsorbFLAG
  !!!CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IAbsorbFLAG",IVariable=IAbsorbFLAG)
  IF(my_rank.EQ.0.AND.IWriteFLAG.EQ.7) PRINT*,"IAbsorbFLAG=",IAbsorbFLAG
  ! -----IAnisoDebyeWallerFactorFlag-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IAnisoDebyeWallerFactorFlag
  !!!CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IAnisoDebyeWallerFactorFlag",IVariable=IAnisoDebyeWallerFactorFlag)
  IF(my_rank.EQ.0.AND.IWriteFLAG.EQ.7) PRINT*,"IAnisoDebyeWallerFactorFlag=",IAnisoDebyeWallerFactorFlag
  ! -----IByteSize-----------------------------------------------------------------depends on system being used, 8 for csc,4 for tinis
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IByteSize
  !!!CALL Message ("ReadInpFile",IInfo,IErr,MessageVariable ="IByteSize",IVariable=IByteSize)
  IF(my_rank.EQ.0.AND.IWriteFLAG.EQ.7) PRINT*,"IByteSize=",IByteSize

  ! ----Two comment lines----------------------# radius of the beam in pixels
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ! -----IPixelCount-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IPixelCount
  IF(my_rank.EQ.0.AND.IWriteFLAG.EQ.7) PRINT*,"IPixelCount=",IPixelCount

  ! ----Two comment lines----------------------# beam selection criteria
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ! -----IMinReflectionPool-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IMinReflectionPool
  IF(my_rank.EQ.0.AND.IWriteFLAG.EQ.7) PRINT*,"IMinReflectionPool=",IMinReflectionPool
  ! -----IMinStrongBeams-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IMinStrongBeams
  IF(my_rank.EQ.0.AND.IWriteFLAG.EQ.7) PRINT*,"IMinStrongBeams=",IMinStrongBeams
  ! -----IMinWeakBeams-----------------------------------------------------------------
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMinWeakBeams
  IF(my_rank.EQ.0.AND.IWriteFLAG.EQ.7) PRINT*,"IMinWeakBeams=",IMinWeakBeams

  ! ----Two comment lines----------------------# crystal settings
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ! -----RDebyeWallerConstant-----------------------------------------------------------------default, if not specified in .cif
  ILine= ILine+1; READ(IChInp,15,ERR=20,END=30) RDebyeWallerConstant
  ! -----RAbsorptionPercentage-----------------------------------------------------------------for proportional model of absorption
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RAbsorptionPercentage

  ! ----Two comment lines----------------------# microscope settings
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ! -----RConvergenceAngle-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,15,ERR=20,END=30) RConvergenceAngle
   ! -----SIncidentBeamDirection-----------------------------------------------------------------IIncidentBeamDirection
  ILine= ILine+1; READ(IChInp,FMT='(27X,A)',END=30) SIncidentBeamDirection
  CALL ThreeDimVectorReadIn(SIncidentBeamDirection,'[',']',RZDirC)
  ! -----SDirectionX-----------------------------------------------------------------IXDirection
  ILine= ILine+1; READ(IChInp,FMT='(27X,A)',END=30) SDirectionX
  CALL ThreeDimVectorReadIn(SDirectionX,'[',']',RXDirC)
  ! -----SNormalDirectionX-----------------------------------------------------------------INormalDirection
  ILine= ILine+1; READ(IChInp,FMT='(27X,A)',ERR=20,END=30) SNormalDirectionX
  CALL ThreeDimVectorReadIn(SNormalDirectionX,'[',']',RNormDirC)
  ! -----RAcceleratingVoltage-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,15,ERR=20,END=30) RAcceleratingVoltage
  ! -----RAcceptanceAngle-----------------------------------------------------------------
  ILine=ILine+1; READ(IChInp,15,ERR=20,END=30) RAcceptanceAngle

  ! ----Two comment lines----------------------# specimen thickness
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ! -----RInitialThickness-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,15,ERR=20,END=30) RInitialThickness
  ! -----RFinalThickness-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,15,ERR=20,END=30) RFinalThickness
  ! -----RDeltaThickness-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,15,ERR=20,END=30) RDeltaThickness
  ! -----INoOfLacbedPatterns-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) INoOfLacbedPatterns

  ! ----Two comment lines----------------------# Refinement Specific Flags
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ! -----IRefineModeFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,FMT='(A)',ERR=20,END=30) SRefineMode
  IF(SCAN(TRIM(ADJUSTL(SRefineMode)),TRIM(ADJUSTL(CAlphabet(19)))).NE.0) THEN
    ISimFLAG=1!Simulation only
  ELSE
    ISimFLAG=0
    SRefineMode = SRefineMode((SCAN(SRefineMode,"=")+1):)
    IRefineMode = 0
    DO ind = 1,IRefinementVariableTypes
      IF(SCAN(TRIM(ADJUSTL(SRefineMode)),TRIM(ADJUSTL(CAlphabet(ind)))).NE.0) THEN
        IRefineMode(ind) = 1
      END IF
    END DO
    IF(my_rank.EQ.0) THEN
      IF(IRefineMode(1).EQ.1) PRINT*, "A:Refining Structure Factors"
      IF(IRefineMode(2).EQ.1) PRINT*, "B:Refining Atomic Coordinates"
      IF(IRefineMode(3).EQ.1) PRINT*, "C:Refining Occupancies "
      IF(IRefineMode(4).EQ.1) PRINT*, "D:Refining Isotropic Debye Waller Factors"
      IF(IRefineMode(5).EQ.1) PRINT*, "E:Refining Anisotropic Debye Waller Factors "
      IF(IRefineMode(6).EQ.1) PRINT*, "F:Refining Lattice Lengths "
      IF(IRefineMode(7).EQ.1) PRINT*, "G:Refining Lattice Angles "
      IF(IRefineMode(8).EQ.1) PRINT*, "H:Refining Convergence Angle"
      IF(IRefineMode(9).EQ.1) PRINT*, "I:Refining Absorption"
      IF(IRefineMode(10).EQ.1) PRINT*,"J:Refining Accelerating Voltage "
      !IF(IRefineMode(11).EQ.1) PRINT*,"K:Refinement by parabola"
      !IF(IRefineMode(12).EQ.1) PRINT*,"L:Refining Structure Factors by bisection"
      IF(ISimFlag.EQ.1) PRINT*,"S:Simulation mode"
    END IF
    !Check if user has requested Ug refinement and anything else which isnt possible
    IF((IRefineMode(1).EQ.1).AND.SUM(IRefineMode).GT.1) THEN         
      IF(my_rank.EQ.0) THEN
        PRINT*,"Structure factors must be refined separately"
      END IF
      IErr = 1
      RETURN
    END IF
  END IF
  ! -----IWeightingFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IWeightingFLAG
  ! -----IMethodFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IMethodFLAG
  IF (my_rank.EQ.0) THEN
      IF(IMethodFLAG.EQ.1) PRINT*, "Refining by Simplex"
      IF(IMethodFLAG.EQ.2) PRINT*, "Refining by Bisection"
      IF(IMethodFLAG.EQ.3) PRINT*, "Refining by Parabola"
  END IF  
  ! -----ICorrelationFLAG: 0=phase,1=sumSq,2=NormalisedCC,3=masked
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) ICorrelationFLAG
  ! -----IImageProcessingFLAG: 0=no processing,1=sqrt,2=log,4=Gaussian blur
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IImageProcessingFLAG
  ! -----RBlurRadius-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,15,ERR=20,END=30) RBlurRadius
  ! -----INoofUgs-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) INoofUgs
  ! -----SAtomicSites-----------------------------------------------------------------
  ILine=ILine+1; READ(IChInp,FMT='(A)',ERR=20,END=30) SAtomicSites
  CALL DetermineRefineableAtomicSites(SAtomicSites,IErr)
  IF( IErr.NE.0 ) THEN
    PRINT*,"ReadInpFile(): error in DetermineRefineableAtomicSites()"
    RETURN
  END IF
  ! -----IPrint-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IPrint
  ! -----RSimplexLengthScale-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,15,ERR=20,END=30) RSimplexLengthScale
  RSimplexLengthScale = RSimplexLengthScale/100.0
  ! -----RExitCriteria-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,15,ERR=20,END=30) RExitCriteria

10 FORMAT(27X,I15.1)
15 FORMAT(27X,F18.9)

  CLOSE(IChInp, ERR=130)

  RETURN

  !	error in OPEN detected
120 IF(my_rank.EQ.0) THEN
     PRINT*,"ReadInpFile(): ERR in OPEN"
     CALL WriteOutInputFile (IErr)
  END IF
  IErr=1
  RETURN
  
  !	error in CLOSE detected
130 IF(my_rank.EQ.0) THEN
     PRINT*,"Input(): ERR in CLOSE"
  END IF
  IErr=1
  RETURN
  
  !	error in READ detected
20 IF(my_rank.EQ.0) THEN
     PRINT*,"Input(): ERRor in READ at line", ILine
     CALL WriteOutInputFile (IErr)
  END IF
  IErr=1
  RETURN
  
  !	EOF in READ occured prematurely
30 IF(my_rank.EQ.0) THEN
     PRINT*,"Input(): EOF in READ at line", ILine
     CALL WriteOutInputFile (IErr)
  END IF
  IErr=1
  RETURN

END SUBROUTINE ReadInpFile

! -----------------------------------------------------------------------

SUBROUTINE ReadScaFile( IErr )
  !This is now redundant
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

  INTEGER(IKIND) :: ILine,IErr,h,k,l,ind,IPos1,IPos2
  CHARACTER*200 :: dummy1,dummy2,SPrintString

  OPEN(Unit = IChInp,FILE="felix.hkl",STATUS='OLD',ERR=10)
  ILine = 0
  IHKLSelectFLAG=1

  !count the number of lines in felix.hkl: this is the number of reflections to output, INoOfLacbedPatterns
  DO
     READ(UNIT= IChInp, END=100, FMT='(a)') dummy1
     ILine=ILine+1
  ENDDO
100 INoOfLacbedPatterns=ILine
  IF (IWriteFLAG.EQ.7.AND.my_rank.EQ.0) THEN
     WRITE(SPrintString,FMT='(I3,A28)') INoOfLacbedPatterns," experimental images to load"
     PRINT*,TRIM(ADJUSTL(SPrintString))
  END IF

  ALLOCATE(RInputHKLs(INoOfLacbedPatterns,ITHREE),STAT=IErr)
  ALLOCATE(IOutputReflections(INoOfLacbedPatterns),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadHklFile(",my_rank,")error in allocations"
     RETURN
  END IF

  !read in the hkls
  REWIND(UNIT=IChInp)
  ILine = 0 
  DO ind = 1,INoOfLacbedPatterns
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

     IF(my_rank.EQ.0.AND.IWriteFLAG.EQ.6) THEN
        WRITE(SPrintString,'(3(I4.1,1X))') NINT(RInputHKLs(ILine,:))
        PRINT*,TRIM(ADJUSTL(SPrintString))
     END IF
  END DO

  RETURN

10 IHKLSelectFLAG=0
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
  CHARACTER*200 :: SAtomicSites,SFormatString,SLengthofNumberString,SPrintString

  IPos1 = SCAN(SAtomicSites,'(')
  IPos2 = SCAN(SAtomicSites,')')
  IF(((IPos2-IPos1).EQ.1).OR.(IPos1.EQ.0).OR.(IPos2.EQ.0)) THEN
     IF(IRefineMode(2).EQ.1) THEN
        IF (my_rank.EQ.0) PRINT*,"You Have Not Specfied Atomic Sites to Refine" 
        IErr = 1
        RETURN
     END IF
  END IF

  IF ((IPos2-IPos1).GT.1.AND.SCAN(SAtomicSites,',').EQ.0) THEN
     ALLOCATE(IAtomicSitesToRefine(1),STAT=IErr)
     IF(IErr.NE.0.AND.(my_rank.EQ.0)) PRINT*,&
       "DetermineRefineableAtomicSites: error allocating IAtomicSitesToRefine"
     IF (my_rank.EQ.0) PRINT*,SIZE(IAtomicSitesToRefine)
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
     IF(IErr.NE.0.AND.(my_rank.EQ.0)) PRINT*,&
       "DetermineRefineableAtomicSites: error allocating IAtomicSitesToRefine"
       
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
  WRITE(SPrintString,FMT='(A15,10(I2,1X))') "Refining atoms ",IAtomicSitesToRefine
  IF(my_rank.EQ.0) PRINT*,TRIM(ADJUSTL(SPrintString))
  
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
  CHARACTER*10 :: h,k,l


  !for IByteSize: 2bytes=64-bit input file (NB tinis specifies in bytes, not bits)
  !NB when this subroutine is working get rid of the pointless variable RImageIn
  ALLOCATE(RImageIn(2*IPixelCount,2*IPixelCount), STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadExperimentalImages(",my_rank,")error allocating RImageIn"
     RETURN
  ENDIF
  RImageIn=ZERO

  DO ind = 1,INoOfLacbedPatterns
     WRITE(h,'(I3.1)')  NINT(RInputHKLs(ind,1))
     WRITE(k,'(I3.1)')  NINT(RInputHKLs(ind,2))
     WRITE(l,'(I3.1)')  NINT(RInputHKLs(ind,3))
     WRITE(filename,*) TRIM(ADJUSTL(h)),TRIM(ADJUSTL(k)),TRIM(ADJUSTL(l)),".img"
     IF (IWriteFLAG.EQ.7.AND.my_rank.EQ.0) THEN
        PRINT*,filename
     END IF
     !WRITE(filename,"(A6,I3.3,A4)") "felix.",ind,".img"  !old version with format felix.000.img
     OPEN(UNIT= IChInImage, ERR= 10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),FORM='UNFORMATTED',&
          ACCESS='DIRECT',IOSTAT=Ierr,RECL=2*IPixelCount*IByteSize)
     DO jnd=1,2*IPixelCount
        READ(IChInImage,rec=jnd,ERR=10) RImageIn(jnd,:)
     END DO
     IF(MINVAL(RImageIn).LT.ZERO.AND.(my_rank.EQ.0)) THEN
        PRINT*,"Warning! There are negative values in your experimental images"
        INegError = INegError + 1
     END IF
     RImageExpi(:,:,ind) = RImageIn
     CLOSE(IChInImage,IOSTAT=IErr) 
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadExperimentalImages (", my_rank, ") error in CLOSE()"
        RETURN
     END IF
  END DO
  DEALLOCATE(RImageIn,STAT=IErr)

  IF (INegError.NE.0) THEN
     IErr = 1
     PRINT*,"No. of Images with Negative Values",INegError
  END IF

  IF (my_rank.EQ.0.AND.IErr.EQ.0) THEN
     WRITE(SPrintString,FMT='(I3,A40)') INoOfLacbedPatterns," experimental images successfully loaded"
     PRINT*,TRIM(ADJUSTL(SPrintString))
  END IF

  RETURN

10 IErr=1
  PRINT*,"ReadExperimentalImages(", my_rank, ")error reading ",TRIM(ADJUSTL(filename)),", line ",jnd
  RETURN
END SUBROUTINE ReadExperimentalImages

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
