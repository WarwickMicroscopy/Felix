!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-17, all rights reserved
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
!  Felix is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  Felix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with Felix.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! All procedures conatained in this file:
! ReadInput( )
! ReadInpFile( )
! ReadHklFile( )
! DetermineRefineableAtomicSites( )
! ReadExperimentalImages( )
! ThreeDimVectorReadIn( )
! WriteOutInputFile( )
! WriteToScreenandFile( )



!>
!! Procedure-description: Calls the various functions which read in all the
!! required data/parameters - read felix.inp, .sca, .cif, .hkl
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE ReadInput(IErr)

  USE MyNumbers
  USE IConst
  USE IPara

  USE MPI
  USE MyMPI
  USE alert_function_mod
  
  IMPLICIT NONE

  INTEGER(IKIND):: IErr

  !Calling all functions which read felix.inp, felix.sca and felix.cif
  !(felix.hkl too depending on user preference)
  !ensure all input files are in working directory
  
  !felix.inp
  CALL ReadInpFile(IErr)
  IF(LALERT(IErr,"ReadInput","ReadInpFile")) RETURN  

  !felix.hkl
  CALL ReadHklFile(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Error:ReadInput(", my_rank, ") error",IErr, &
          "in ReadHklFile()"
     RETURN
  ENDIF

  !felix.sca
  CALL ScatteringFactors(IScatterFactorMethodFLAG,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Error:ReadInput(", my_rank, ") error",IErr, &
          " in ReadScaFile()"
     RETURN
  ENDIF

  !felix.cif
  CALL ReadCif(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Error:ReadInput(", my_rank, ") error",IErr, &
          "in ReadCif()"
     RETURN
  ENDIF

END SUBROUTINE ReadInput

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description:
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE ReadInpFile( IErr )

  USE MyNumbers

  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE message_mod
  USE MPI
  USE MyMPI
  USE alert_function_mod

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
  CALL message ( LL, dbg7, "IScatterFactorMethodFLAG=",IScatterFactorMethodFLAG )
  ! -----IMaskFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IMaskFLAG
  CALL message ( LL, dbg3, "IMaskFLAG=",IMaskFLAG)
  ! -----IHolzFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IHolzFLAG
  CALL message ( LL, dbg3, "IHolzFLAG=",IHolzFLAG)
  ! -----IAbsorbFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IAbsorbFLAG
  CALL message ( LL, dbg3, "IAbsorbFLAG=",IAbsorbFLAG)
  ! -----IAnisoDebyeWallerFactorFlag-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IAnisoDebyeWallerFactorFlag
  CALL message ( LL, dbg3, "IAnisoDebyeWallerFactorFlag=",IAnisoDebyeWallerFactorFlag)
  ! -----IByteSize-----------------------------------------------------------------depends on system being used, 8 for csc,4 for tinis
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IByteSize
  CALL message ( LL, dbg3, "IByteSize=",IByteSize)

  ! ----Two comment lines----------------------# radius of the beam in pixels
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ! -----IPixelCount-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IPixelCount
  CALL message ( LL, dbg3, "IPixelCount=",IPixelCount)

  ! ----Two comment lines----------------------# beam selection criteria
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ! -----IMinReflectionPool-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IMinReflectionPool
  CALL message ( LL, dbg3, "IMinReflectionPool=",IMinReflectionPool)
  ! -----IMinStrongBeams-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IMinStrongBeams
  CALL message ( LL, dbg3, "IMinStrongBeams=",IMinStrongBeams)
  ! -----IMinWeakBeams-----------------------------------------------------------------
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMinWeakBeams
  CALL message ( LL, dbg3, "IMinWeakBeams=",IMinWeakBeams)

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
    IF(IRefineMode(1).EQ.1) CALL message( LS, "Mode A selected: Refining Structure Factors" )
    IF(IRefineMode(2).EQ.1) CALL message( LS, "Mode B selected: Refining Atomic Coordinates" )
    IF(IRefineMode(3).EQ.1) CALL message( LS, "Mode C selected: Refining Occupancies ")
    IF(IRefineMode(4).EQ.1) CALL message( LS, "Mode D selected: Refining Isotropic Debye Waller Factors")
    IF(IRefineMode(5).EQ.1) CALL message( LS, "Mode E selected: Refining Anisotropic Debye Waller Factors ")
    IF(IRefineMode(6).EQ.1) CALL message( LS, "Mode F selected: Refining Lattice Lengths ")
    IF(IRefineMode(7).EQ.1) CALL message( LS, "Mode G selected: Refining Lattice Angles ")
    IF(IRefineMode(8).EQ.1) CALL message( LS, "Mode H selected: Refining Convergence Angle")
    IF(IRefineMode(9).EQ.1) CALL message( LS, "Mode I selected: Refining Absorption")
    IF(IRefineMode(10).EQ.1) CALL message( LS, "Mode J selected: Refining Accelerating Voltage ")
    IF(ISimFlag.EQ.1) CALL message( LS, "Mode S selected: Simulation mode")
    !Check if user has requested Ug refinement and anything else which isn't possible
    IF((IRefineMode(1).EQ.1).AND.SUM(IRefineMode).GT.1) THEN         
      IF(my_rank.EQ.0) THEN
        PRINT*,"Error:Structure factors must be refined separately"
      END IF
      IErr = 1
      RETURN
    END IF
  END IF
  ! -----IWeightingFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IWeightingFLAG
  ! -----IMethodFLAG-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IMethodFLAG
  IF(IMethodFLAG.EQ.1) CALL message( LS, "Mode selected: Refining by simplex")
  IF(IMethodFLAG.EQ.2) CALL message( LS, "Mode selected: Refining by maximum gradient")
  IF(IMethodFLAG.EQ.3) CALL message( LS, "Mode selected: Refining by pairwise maximum gradient")
 
  ! -----ICorrelationFLAG: 0=phase,1=sumSq,2=NormalisedCC,3=masked
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) ICorrelationFLAG
  ! -----IImageProcessingFLAG: 0=no processing,1=sqrt,2=log
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) IImageProcessingFLAG
  ! -----RBlurRadius-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,15,ERR=20,END=30) RBlurRadius
  ! -----INoofUgs-----------------------------------------------------------------
  ILine= ILine+1; READ(IChInp,10,ERR=20,END=30) INoofUgs
  ! -----SAtomicSites-----------------------------------------------------------------
  ILine=ILine+1; READ(IChInp,FMT='(A)',ERR=20,END=30) SAtomicSites
  CALL DetermineRefineableAtomicSites(SAtomicSites,IErr)
  IF( IErr.NE.0 ) THEN
    PRINT*,"Error:ReadInpFile(): error in DetermineRefineableAtomicSites()"
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
     PRINT*,"Error:ReadInpFile(): ERR in OPEN"
     CALL WriteOutInputFile (IErr)
  END IF
  IErr=1
  RETURN
  
  !	error in CLOSE detected
130 IF(my_rank.EQ.0) THEN
     PRINT*,"Error:Input(): ERR in CLOSE"
  END IF
  IErr=1
  RETURN
  
  !	error in READ detected
20 IErr=20;IF(LALERT(IErr,"ReadInpFile","READ on felix.inp")) RETURN
!20 IF(my_rank.EQ.0) THEN
!     PRINT*,"Error:Input(): ERRor in READ at line", ILine
!     CALL WriteOutInputFile (IErr)
!  END IF
!  IErr=1
!  RETURN
  
  !	EOF in READ occured prematurely
30 IF(my_rank.EQ.0) THEN
     PRINT*,"Error:Input(): EOF in READ at line", ILine
     CALL WriteOutInputFile (IErr)
  END IF
  IErr=1
  RETURN

END SUBROUTINE ReadInpFile

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!>
!! Procedure-description:
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE ReadHklFile(IErr)

  USE MyNumbers
  USE IChannels

  USE IConst
  USE RConst

  USE IPara
  USE RPara
  USE CPara
  USE SPara
  USE message_mod
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

  CALL message ( LL, dbg7, "Number of experimental images to load = ", INoOfLacbedPatterns)

  ALLOCATE(RInputHKLs(INoOfLacbedPatterns,ITHREE),STAT=IErr)
  ALLOCATE(IOutputReflections(INoOfLacbedPatterns),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Error:ReadHklFile(",my_rank,")error in allocations"
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

     CALL message ( LL, dbg7, "RInputHKLs", NINT(RInputHKLs(ILine,:)) )

  END DO

  RETURN

10 IHKLSelectFLAG=0
  RETURN

END SUBROUTINE ReadHklFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description:
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE DetermineRefineableAtomicSites(SAtomicSites,IErr)

  USE MyNumbers

  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE message_mod
  USE MPI
  USE MyMPI

  IMPLICIT NONE  

  INTEGER(IKIND) :: IPos,IPos1,IPos2,IErr,ind
  CHARACTER*200 :: SAtomicSites,SFormatString,SLengthofNumberString

  IPos1 = SCAN(SAtomicSites,'(')
  IPos2 = SCAN(SAtomicSites,')')
  IF(((IPos2-IPos1).EQ.1).OR.(IPos1.EQ.0).OR.(IPos2.EQ.0)) THEN
     IF(IRefineMode(2).EQ.1) THEN
        IF (my_rank.EQ.0) PRINT*,"Error:You Have Not Specfied Atomic Sites to Refine" 
        IErr = 1
        RETURN
     END IF
  END IF

  IF ((IPos2-IPos1).GT.1.AND.SCAN(SAtomicSites,',').EQ.0) THEN
     ALLOCATE(IAtomicSitesToRefine(1),STAT=IErr)
     IF(IErr.NE.0.AND.(my_rank.EQ.0)) PRINT*,"Error:",&
       "DetermineRefineableAtomicSites: error allocating IAtomicSitesToRefine"
     CALL message (LM, "SIZE(IAtomicSitesToRefine) = ",SIZE(IAtomicSitesToRefine) )
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
     IF(IErr.NE.0.AND.(my_rank.EQ.0)) PRINT*,"Error:",&
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
  CALL message (LM, "Refining atoms ",IAtomicSitesToRefine )
  
END SUBROUTINE DetermineRefineableAtomicSites

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description:
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE ReadExperimentalImages(IErr)

  USE MyNumbers

  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels
  USE message_mod
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: ind,jnd,IErr
  INTEGER(IKIND) :: INegError = 0
  CHARACTER*34 :: filename
  CHARACTER*200 :: SPrintString
  CHARACTER*10 :: intstring

  !for IByteSize: 2bytes=64-bit input file (NB tinis specifies in bytes, not bits)
  !NB when this subroutine is working get rid of the pointless variable RImageIn
  ALLOCATE(RImageIn(2*IPixelCount,2*IPixelCount), STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Error:ReadExperimentalImages(",my_rank,")error allocating RImageIn"
     RETURN
  ENDIF
  RImageIn=ZERO

  DO ind = 1,INoOfLacbedPatterns  
    ! An image expected for each LacbedPattern
    ! Write corresponding filenames including chemical formula
    WRITE(filename,*) TRIM(ADJUSTL(chemicalformula)),"_"
    DO jnd = 1,3
    WRITE(intstring,'(I3.1)')  NINT(RInputHKLs(ind,jnd))
    IF (NINT(RInputHKLs(ind,jnd))>= 0) THEN    
      filename = TRIM(filename) // '+'
    END IF
    filename = TRIM(filename) // TRIM(ADJUSTL(IntString))
    END DO
    filename = TRIM(filename) // '.img'
    
    CALL message( LL, dbg7, "filename = ", filename)
    OPEN(UNIT= IChInImage, ERR= 10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),FORM='UNFORMATTED',&
          ACCESS='DIRECT',IOSTAT=Ierr,RECL=2*IPixelCount*IByteSize)
    DO jnd=1,2*IPixelCount
      READ(IChInImage,rec=jnd,ERR=10) RImageIn(jnd,:)
    END DO
    RImageExpi(:,:,ind) = RImageIn
    CLOSE(IChInImage,IOSTAT=IErr) 
    IF( IErr.NE.0 ) THEN
      PRINT*,"Error:ReadExperimentalImages (", my_rank, ") error in CLOSE()"
      RETURN
    END IF
  END DO
  DEALLOCATE(RImageIn,STAT=IErr)

  CALL message ( LM, "Number of experimental images successfully loaded = ", INoOfLacbedPatterns )

  RETURN

10 IErr=1
  PRINT*,"Error:Error Message:"
  PRINT*,"Error:ReadExperimentalImages(", my_rank, ")error reading ",TRIM(ADJUSTL(filename)),", line ",jnd
  PRINT*,"Error:From felix.cif ChemicalFormula = '","ChemicalFormula","'"
  PRINT*,"Error:Expect input .img files in a format like ChemicalFormula_-2-2+0.img"
  PRINT*,"Error:e.g. GaAs_+0+0-8.img"

  RETURN

END SUBROUTINE ReadExperimentalImages

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description: Function that reads in a 3D vector from a file,
!! which has to have the form [I1,I2,I3] where I1,I2,I3 are Integers of any length.
!! Used to read-in integer string vectors from felix.inp into real vectors 
!!
!! Major-Authors: Alexander Hubert (2015)
!!
SUBROUTINE ThreeDimVectorReadIn(SUnformattedVector,SOpenBracketDummy, &
     SCloseBracketDummy,RFormattedVector)

  ! e.g. CALL ThreeDimVectorReadIn(SIncidentBeamDirection,'[',']',RZDirC)
  ! Current format has interesting bracket inputs

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI 
 
  IMPLICIT NONE

  CHARACTER(*), INTENT(IN) :: &
       SUnformattedVector,SOpenBracketDummy,SCloseBracketDummy
  REAL(RKIND),INTENT(OUT),DIMENSION(3) :: &
       RFormattedVector

  CHARACTER*1 :: &
       SComma=',',SOpenBracket,SCloseBracket
  CHARACTER*100 :: &
       SFormattedVectorX,SFormattedVectorY,SFormattedVectorZ
 
  LOGICAL :: &
       LBACK=.TRUE.
  
  INTEGER(IKIND) :: &
       IOpenBracketPosition, ICloseBracketPosition, IFirstCommaPosition, &
       ILastCommaPosition
  
  ! Trim and adjustL bracket to ensure one character string
  SOpenBracket=TRIM(ADJUSTL(SOpenBracketDummy))
  SCloseBracket=TRIM(ADJUSTL(SCloseBracketDummy))
  
  ! Read in string to for the incident beam direction, convert these
  ! to x,y,z coordinate integers with string manipulation
  
  IOpenBracketPosition=INDEX(SUnformattedVector,SOpenBracket)
  ICloseBracketPosition=INDEX(SUnformattedVector,SCloseBracket)
  IFirstCommaPosition=INDEX(SUnformattedVector,SComma)
  ILastCommaPosition=INDEX(SUnformattedVector,SComma,LBACK)

  ! Separate the Unformatted Vector String into its three X,Y,Z components
  SFormattedVectorX=TRIM(ADJUSTL(SUnformattedVector(IOpenBracketPosition+1:IFirstCommaPosition-1)))
  SFormattedVectorY=TRIM(ADJUSTL(SUnformattedVector(IFirstCommaPosition+1:ILastCommaPosition-1)))
  SFormattedVectorZ=TRIM(ADJUSTL(SUnformattedVector(ILastCommaPosition+1:ICloseBracketPosition-1)))

  ! Read in each string component into its associated position in real variable RFormattedVector
  READ(SFormattedVectorX,*) RFormattedVector(1)
  READ(SFormattedVectorY,*) RFormattedVector(2)
  READ(SFormattedVectorZ,*) RFormattedVector(3)

END SUBROUTINE ThreeDimVectorReadIn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description: Write out the sample input file, when none provided
!!
!! Major-Authors: Keith Evans (2014)
!! 
SUBROUTINE WriteOutInputFile (IErr)
  
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
     
     OPEN(UNIT= IChInp,FILE= "felix.inp.sample",&
       STATUS= 'UNKNOWN')
     !todo - look into this subroutine
     CALL WriteToScreenandFile(ADJUSTL("# Input file for felixsim/draw/refine version :VERSION: Build :BUILD:"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# ------------------------------------"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# ------------------------------------"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# control flags"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IWriteFLAG                = 3"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IImageFLAG                = 1"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IScatterFactorMethodFLAG  = 0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IMaskFLAG                 = 1"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IHolzFLAG                 = 0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IAbsorbFLAG               = 1"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IAnisoDebyeWallerFlag     = 0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# radius of the beam in pixels"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IPixelCount               = 64"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# beam selection criteria"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IMinReflectionPool        = 100"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IMinStrongBeams           = 20"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IMinWeakBeams             = 0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# crystal settings"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RDebyeWallerConstant      = 0.4668"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RAbsorptionPer            = 2.9"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# microscope settings"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("ROuterConvergenceAngle    = 6.0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IIncidentBeamDirection    = [0,1,1]"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IXDirection               = [1,0,0]"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("INormalDirection          = [0,1,1]"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RAcceleratingVoltage (kV) = 200.0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RAcceptanceAngle          = 0.0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# Image Output Options"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RInitialThickness        = 400.0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RFinalThickness          = 700.0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RDeltaThickness          = 10.0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("INoOfLacbedPatterns              = 7"),IErr)

     IF(ISoftwareMode.EQ.2) THEN
        CALL WriteToScreenandFile(ADJUSTL(""),IErr)
        CALL WriteToScreenandFile(ADJUSTL("#Refinement Specific Flags"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("IRefineModeFLAG          = B"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("IWeightingFLAG           = 0"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("IContinueFLAG            = 0"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("ICorrelationFLAG         = 0"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("IImageProcessingFLAG     = 0"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("RBlurRadius              = 1.45"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("INoofUgs                 = 10"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("IAtomicSites             = (1,2,6)"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("IPrint                   = 10"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("RSimplexLengthScale      = 5.0"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("RExitCriteria            = 0.0001"),IErr)
     END IF
        CLOSE(UNIT=IChInp)
        
END SUBROUTINE WriteOutInputFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description: Prints to screen and file
!!
!! Major-Authors: 
!!  
SUBROUTINE WriteToScreenandFile(SStringtoWrite,IErr)
  !This is a pointless subroutine
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
  CHARACTER(*) :: SStringtoWrite

  PRINT*,TRIM(ADJUSTL(SStringtoWrite))
  WRITE(UNIT=IChInp,FMT='(A)') TRIM(ADJUSTL(SStringtoWrite))

 END SUBROUTINE WriteToScreenandFile

