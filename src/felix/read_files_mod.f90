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
! ReadInpFile( )
! ReadHklFile( )
! ReadExperimentalImages( )
! DetermineRefineableAtomicSites( )
! ThreeDimVectorReadIn( )

MODULE read_files_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ReadInpFile, ReadHklFile, ReadExperimentalImages

  CONTAINS

  !>
  !! Procedure-description:
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE ReadInpFile ( IErr )

    !?? this should contain clear input information
    !?? github and other places should referance here

    USE MyNumbers
    USE message_mod

    ! global outputs, read from .inp
    USE IPARA, ONLY : IWriteFLAG, IImageFLAG, IScatterFactorMethodFLAG, IMaskFLAG, IHolzFLAG, &
          IAbsorbFLAG, IAnisoDebyeWallerFactorFlag, IByteSize, IMinReflectionPool, &
          IMinStrongBeams, IMinWeakBeams, INoOfLacbedPatterns, ISimFLAG, IRefineMode, &
          IWeightingFLAG, IMethodFLAG, ICorrelationFLAG, IImageProcessingFLAG, &
          INoofUgs, IPrint, IPixelCount
    USE RPARA, ONLY : RDebyeWallerConstant, RAbsorptionPercentage, RConvergenceAngle, &
          RZDirC, RXDirC, RNormDirC, RAcceleratingVoltage, RAcceptanceAngle, &
          RInitialThickness, RFinalThickness, RDeltaThickness, RBlurRadius, &
          RSimplexLengthScale, RExitCriteria

    ! global inputs
    USE IChannels, ONLY : IChInp
    USE IPARA, ONLY : IRefinementVariableTypes
    USE SConst, ONLY : SAlphabet

    IMPLICIT NONE


    INTEGER(IKIND),INTENT(OUT) :: IErr
    INTEGER(IKIND) :: ILine, ind, IPos
    REAL(RKIND) :: ROfIter
    CHARACTER(200) :: SImageMode, SElements, SRefineMode, SStringFromNumber, SRefineYESNO, &
          SAtomicSites, SFormatString, SLengthofNumberString, &
          SDirectionX, SIncidentBeamDirection, SNormalDirectionX, SPrintString

    OPEN(UNIT= IChInp, IOSTAT=IErr, FILE= "felix.inp",STATUS= 'OLD')
    IF(l_alert(IErr,"ReadInpFile","OPEN() felix.inp")) RETURN
    ILine= 1

    ! There are six introductory comment lines which are ignored
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)') 
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

    !--------------------------------------------------------------------
    ! control flags
    !--------------------------------------------------------------------

    ! IWriteFLAG
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IWriteFLAG

    ! IImageFLAG
    ILine= ILine+1; READ(IChInp,FMT='(A)',ERR=20,END=30) SImageMode
    ! IImageFLAG handling - select which type(s) of output images to be produce
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

    ! IScatterFactorMethodFLAG
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IScatterFactorMethodFLAG
    CALL message ( LXL, dbg7, "IScatterFactorMethodFLAG=",IScatterFactorMethodFLAG )
    ! IMaskFLAG
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IMaskFLAG
    CALL message ( LXL, dbg3, "IMaskFLAG=",IMaskFLAG)
    ! IHolzFLAG
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IHolzFLAG
    CALL message ( LXL, dbg3, "IHolzFLAG=",IHolzFLAG)
    ! IAbsorbFLAG
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IAbsorbFLAG
    CALL message ( LXL, dbg3, "IAbsorbFLAG=",IAbsorbFLAG)
    ! IAnisoDebyeWallerFactorFlag
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IAnisoDebyeWallerFactorFlag
    CALL message ( LXL, dbg3, "IAnisoDebyeWallerFactorFlag=",IAnisoDebyeWallerFactorFlag)
    ! IByteSize
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IByteSize
    CALL message ( LXL, dbg3, "IByteSize=",IByteSize) !?? depends on system, 8 for csc, 4 tinis

    !--------------------------------------------------------------------
    ! radius of the beam in pixels
    !--------------------------------------------------------------------

    ! Two comment lines
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ! IPixelCount
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IPixelCount
    CALL message ( LXL, dbg3, "IPixelCount=",IPixelCount)

    !--------------------------------------------------------------------
    ! beam selection criteria
    !--------------------------------------------------------------------

    ! Two comment lines 
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ! IMinReflectionPool
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IMinReflectionPool
    CALL message ( LXL, dbg3, "IMinReflectionPool=",IMinReflectionPool)
    ! IMinStrongBeams
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IMinStrongBeams
    CALL message ( LXL, dbg3, "IMinStrongBeams=",IMinStrongBeams)
    ! IMinWeakBeams
    ILine= ILine+1
    READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IMinWeakBeams
    CALL message ( LXL, dbg3, "IMinWeakBeams=",IMinWeakBeams)

    !--------------------------------------------------------------------
    ! crystal settings
    !--------------------------------------------------------------------

    ! Two comment lines          
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ! RDebyeWallerConstant          !?? default, if not specified in .cif
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RDebyeWallerConstant
    ! RAbsorptionPercentage         !?? for proportional model of absorption
    ILine= ILine+1
    READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RAbsorptionPercentage

    !--------------------------------------------------------------------
    ! microscope settings
    !--------------------------------------------------------------------

    ! Two comment lines
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ! RConvergenceAngle
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RConvergenceAngle

    ! SIncidentBeamDirection -> IIncidentBeamDirection
    ILine= ILine+1; READ(IChInp,FMT='(27X,A)',END=30) SIncidentBeamDirection
    CALL ThreeDimVectorReadIn(SIncidentBeamDirection,'[',']',RZDirC)
    ! SDirectionX -> IXDirection
    ILine= ILine+1; READ(IChInp,FMT='(27X,A)',END=30) SDirectionX
    CALL ThreeDimVectorReadIn(SDirectionX,'[',']',RXDirC)
    ! SNormalDirectionX -> INormalDirection
    ILine= ILine+1; READ(IChInp,FMT='(27X,A)',ERR=20,END=30) SNormalDirectionX
    CALL ThreeDimVectorReadIn(SNormalDirectionX,'[',']',RNormDirC)

    ! RAcceleratingVoltage
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RAcceleratingVoltage
    ! RAcceptanceAngle
    ILine=ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RAcceptanceAngle

    !--------------------------------------------------------------------
    ! Image Output Options
    !--------------------------------------------------------------------

    !?? update input files and section here to 'specimin thickness'

    ! Two comment lines
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ! RInitialThickness
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RInitialThickness
    ! RFinalThickness
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RFinalThickness
    ! RDeltaThickness
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RDeltaThickness
    ! INoOfLacbedPatterns             !?? update IReflectOut in felix.inp files
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) INoOfLacbedPatterns 

    !--------------------------------------------------------------------
    ! Refinement Specific Flags
    !--------------------------------------------------------------------

    ! Two comment lines
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')  
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

    ! IRefineModeFLAG
    ILine= ILine+1; READ(IChInp,FMT='(A)',ERR=20,END=30) SRefineMode
    IF(SCAN(TRIM(ADJUSTL(SRefineMode)),TRIM(ADJUSTL(SAlphabet(19)))).NE.0) THEN
      ISimFLAG=1 ! Simulation only
    ELSE
      ISimFLAG=0
      SRefineMode = SRefineMode((SCAN(SRefineMode,"=")+1):)
      IRefineMode = 0
      DO ind = 1,IRefinementVariableTypes
        IF(SCAN(TRIM(ADJUSTL(SRefineMode)),TRIM(ADJUSTL(SAlphabet(ind)))).NE.0) THEN
          IRefineMode(ind) = 1
        END IF
      END DO
      ! Check which refinement modes have been selected
      IF(IRefineMode(1) .EQ.1) CALL message( LS, "Refining Structure Factors")
      IF(IRefineMode(2) .EQ.1) CALL message( LS, "Refining Atomic Coordinates")
      IF(IRefineMode(3) .EQ.1) CALL message( LS, "Refining Occupancies ")
      IF(IRefineMode(4) .EQ.1) CALL message( LS, "Refining Isotropic Debye Waller Factors")
      IF(IRefineMode(5) .EQ.1) CALL message( LS, "Refining Anisotropic Debye Waller Factors ")
      IF(IRefineMode(6) .EQ.1) CALL message( LS, "Refining Lattice Lengths ")
      IF(IRefineMode(7) .EQ.1) CALL message( LS, "Refining Lattice Angles ")
      IF(IRefineMode(8) .EQ.1) CALL message( LS, "Refining Convergence Angle")
      IF(IRefineMode(9) .EQ.1) CALL message( LS, "Refining Absorption")
      IF(IRefineMode(10).EQ.1) CALL message( LS, "Refining Accelerating Voltage ")
      !?? ISimFlag should always be zero here
      IF(ISimFlag.EQ.1) CALL message( LS, "Simulation only")
      ! Error Check - user cannot request Ug refinement and anything else
      IF((IRefineMode(1).EQ.1).AND.SUM(IRefineMode).GT.1) THEN
        IErr = 1; IF(l_alert(IErr,"ReadInpFile",&
              "Structure factors must be refined separately")) RETURN
      END IF
    END IF

    ! IWeightingFLAG          !?? is this still used JR
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IWeightingFLAG
    ! IMethodFLAG
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IMethodFLAG
    IF(IMethodFLAG.EQ.1) CALL message( LS, "Refining by simplex")
    IF(IMethodFLAG.EQ.2) CALL message( LS, "Refining by maximum gradient")
    IF(IMethodFLAG.EQ.3) CALL message( LS, "Refining by pairwise maximum gradient")
   
    ! ICorrelationFLAG: 0=phase,1=sumSq,2=NormalisedCC,3=masked
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) ICorrelationFLAG
    ! IImageProcessingFLAG: 0=no processing,1=sqrt,2=log
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IImageProcessingFLAG
    ! RBlurRadius
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RBlurRadius
    ! INoofUgs
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) INoofUgs
    ! SAtomicSites !?? what global does this later set?
    ILine=ILine+1; READ(IChInp,FMT='(A)',ERR=20,END=30) SAtomicSites
    CALL DetermineRefineableAtomicSites(SAtomicSites,IErr) !?? what does this do? JR
    IF(l_alert(IErr,"ReadInpFile","DetermineRefineableAtomicSites()")) RETURN
    ! IPrint
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IPrint
    ! RSimplexLengthScale
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RSimplexLengthScale
    RSimplexLengthScale = RSimplexLengthScale/100.0
    ! RExitCriteria
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RExitCriteria

    !--------------------------------------------------------------------
    ! finish reading, close felix.inp
    !--------------------------------------------------------------------
    
    CLOSE(IChInp, IOSTAT=IErr)
    IF(l_alert(IErr,"ReadInpFile","CLOSE() felix.inp")) RETURN
    RETURN

    !--------------------------------------------------------------------
    ! GOTO error handling from READ()
    !--------------------------------------------------------------------

    !	error in READ() detected  
  20 CONTINUE
    IErr=1
    WRITE(SPrintString,*) ILine
    IF(l_alert(IErr,"ReadInpFile",&
          "READ() felix.inp line number ="//TRIM(SPrintString))) RETURN
    
    !	EOF in READ() occured prematurely
  30 CONTINUE
    IErr=1
    WRITE(SPrintString,*) ILine
    IF(l_alert( IErr,"ReadInpFile",&
          "READ() felix.inp, premature end of file, line number =" // &
          TRIM(SPrintString) )) RETURN

    !?? JR Old write example felix.inp to screen removed
    !?? JR could direct to github wiki page, it currently has clear information
    !?? JR should be able to direct user to felix.inp in sample folders
    !?? JR here should contain clear information about inputs and github can referance/link this

  END SUBROUTINE ReadInpFile

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




  !>
  !! Procedure-description:
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE ReadHklFile(IErr)

    USE MyNumbers
    USE message_mod

    ! global outputs
    USE RPARA, ONLY : RInputHKLs
    USE IPARA, ONLY : INoOfLacbedPatterns, IHKLSelectFLAG, &
          IOutputReflections !?? allocated here

    ! global inputs
    USE IChannels, ONLY : IChInp

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(OUT) :: IErr
    INTEGER(IKIND) :: ILine, h, k, l, ind, IPos1, IPos2
    CHARACTER*200 :: dummy1, dummy2

    OPEN(Unit = IChInp,FILE="felix.hkl",STATUS='OLD',ERR=10)
    ILine = 0
    IHKLSelectFLAG=1

    ! count the number of lines in felix.hkl:
    ! this is the number of reflections to output, INoOfLacbedPatterns
    DO
      READ(UNIT= IChInp, END=100, FMT='(a)') dummy1
      ILine = ILine+1
    ENDDO  !?? JR can we change this? IOSTAT EOF should be identifiable, IF EOF -> EXIT
  100 INoOfLacbedPatterns = ILine
    CALL message ( LXL, dbg7, "Number of experimental images to load = ", INoOfLacbedPatterns)

    ALLOCATE(RInputHKLs(INoOfLacbedPatterns,ITHREE),STAT=IErr)
    IF(l_alert(IErr,"ReadHklFile","allocate RInputHKLs")) RETURN
    ALLOCATE(IOutputReflections(INoOfLacbedPatterns),STAT=IErr) !?? should we allocate here JR
    IF(l_alert(IErr,"ReadHklFile","allocate IOutputReflections")) RETURN

    ! read in the hkls
    REWIND(UNIT=IChInp) ! goes to beggining of felix.hkl file
    ILine = 0 
    DO ind = 1,INoOfLacbedPatterns
      ! READ HKL in as String
      READ(UNIT= IChInp,FMT='(A)' ) dummy1
      ! Scan String for h
      IPos1 = SCAN(dummy1,'[')
      IPos2 = SCAN(dummy1,',')
      dummy2 = dummy1((IPos1+1):(IPos2-1))
      READ(dummy2,'(I20)') h
      ! Scan String for k   
      IPos1 = SCAN(dummy1((IPos2+1):),',') + IPos2
      dummy2 = dummy1((IPos2+1):(IPos1-1))
      READ(dummy2,'(I20)') k
      ! Scan String for l     
      IPos2 = SCAN(dummy1((IPos1+1):),']') + IPos1
      dummy2 = dummy1((IPos1+1):(IPos2-1))
      READ(dummy2,'(I20)') l
      ! Convert to REALs
      ILine=ILine+1
      RInputHKLs(ILine,1) = REAL(h,RKIND)
      RInputHKLs(ILine,2) = REAL(k,RKIND)
      RInputHKLs(ILine,3) = REAL(l,RKIND)   
      CALL message ( LXL, dbg7, "RInputHKLs", NINT(RInputHKLs(ILine,:)) )
    END DO

    RETURN

  10 CONTINUE
    IHKLSelectFLAG=0 !?? JR hkl not found, elaborate what is used instead
    CALL message( LL, "felix.hkl not found, continuing")
    RETURN

  END SUBROUTINE ReadHklFile

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description:
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE ReadExperimentalImages(IErr)

    USE MyNumbers
    USE message_mod

    ! global outputs
    USE RPARA, ONLY : RImageExpi

    ! global inputs
    USE IChannels, ONLY : IChInImage
    USE SPARA, ONLY : SChemicalFormula
    USE IPARA, ONLY : INoOfLacbedPatterns, IPixelCount, IByteSize
    USE RPARA, ONLY : RInputHKLs

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(OUT) :: IErr
    INTEGER(IKIND) :: ind, jnd, INegError = 0
    CHARACTER :: filename*34, SPrintString*200, IntString*10

    ! for IByteSize: 2bytes=64-bit input file (NB tinis specifies in bytes, not bits)

    DO ind = 1,INoOfLacbedPatterns  !?? JR should standardise this filename writing throughout
      ! An image expected for each LacbedPattern
      ! Write corresponding filenames including chemical formula
      WRITE(filename,*) TRIM(ADJUSTL(SChemicalFormula)),"_"
      DO jnd = 1,3
      WRITE(intstring,'(I3.1)')  NINT(RInputHKLs(ind,jnd))
      IF (NINT(RInputHKLs(ind,jnd))>= 0) THEN
        filename = TRIM(filename) // '+'
      END IF
      filename = TRIM(filename) // TRIM(ADJUSTL(IntString))
      END DO
      filename = TRIM(filename) // '.img'

      CALL message(LL, dbg7, "filename = ", filename)
      OPEN(UNIT= IChInImage, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)), &
            FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=2*IPixelCount*IByteSize)
      IF(l_alert(IErr,"ReadExperimentalImages",&
            "OPEN() an experimental image, filename ="//TRIM(ADJUSTL(filename)))) RETURN
      !?? input image filenames should follow format like 'GaAs_+0+0-8.img'
      !?? Chemical formula should match _chemical_formula_structural from felix.cif
      DO jnd=1,2*IPixelCount
        READ(IChInImage,rec=jnd,IOSTAT=IErr) RImageExpi(jnd,:,ind)
        IF(l_alert(IErr,"ReadExperimentalImages",&
              "OPEN() an experimental image, filename ="//TRIM(ADJUSTL(filename)))) RETURN
      END DO
      CLOSE(IChInImage,IOSTAT=IErr)
      IF(l_alert(IErr,"ReadExperimentalImages","CLOSE() an experimental input image")) RETURN
    END DO

    CALL message(LM,"Number of experimental images successfully loaded =",INoOfLacbedPatterns)

    RETURN

  END SUBROUTINE ReadExperimentalImages

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description:
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE DetermineRefineableAtomicSites(SAtomicSites,IErr)

    USE MyNumbers
    USE message_mod

    ! global outputs
    USE IPARA, ONLY : IAtomsToRefine
    
    ! global inputs
    USE IPARA, ONLY : IRefineMode 

    IMPLICIT NONE  

    CHARACTER(200), INTENT(IN) :: SAtomicSites
    INTEGER(IKIND),INTENT(OUT) :: IErr
    INTEGER(IKIND) :: IPos, IPos1, IPos2, ind
    CHARACTER(200) :: SFormatString, SLengthofNumberString

    IPos1 = SCAN(SAtomicSites,'(')
    IPos2 = SCAN(SAtomicSites,')')
    ! error check
    IF(((IPos2-IPos1).EQ.1).OR.(IPos1.EQ.0).OR.(IPos2.EQ.0)) THEN 
      IF(IRefineMode(2).EQ.1) IErr = 1
      IF(l_alert(IErr,"DetermineRefineableAtomicSites",&
              "You Have Not Specfied Atomic Sites to Refine")) RETURN
    END IF

    IF ((IPos2-IPos1).GT.1.AND.SCAN(SAtomicSites,',').EQ.0) THEN

      ALLOCATE(IAtomsToRefine(1),STAT=IErr)
      IF(l_alert(IErr,"DetermineRefineableAtomicSites","allocate IAtomsToRefine")) RETURN
      CALL message (LM, "SIZE(IAtomsToRefine) = ",SIZE(IAtomsToRefine) )
      WRITE(SLengthofNumberString,*) LEN(SAtomicSites((IPos1+1):(IPos2-1))) 
      WRITE(SFormatString,*) "(I"//TRIM(ADJUSTL(SLengthofNumberString))//")"
      READ(SAtomicSites((IPos1+1):(IPos2-1)),FMT=SFormatString) IAtomsToRefine(1)
    ELSE
      IPos = 1
      DO 
        IF(SCAN(SAtomicSites(IPos1:IPos2),',').NE.0) THEN
          IPos1 = IPos1 + LEN(SAtomicSites(IPos1:(IPos1+SCAN(SAtomicSites(IPos1:IPos2),','))))
          IPos = IPos+1
        END IF
        IF (IPos2-IPos1.LE.1) EXIT
      END DO

      ALLOCATE(IAtomsToRefine(IPos),STAT=IErr)
      IF(l_alert(IErr,"DetermineRefineableAtomicSites","allocate IAtomsToRefine")) RETURN
       
      IPos1 = SCAN(SAtomicSites,'(')
      DO ind = 1,SIZE(IAtomsToRefine,DIM=1)
        IF(SCAN(SAtomicSites((IPos1+1):IPos2),',').NE.0) THEN
          IPos = SCAN(SAtomicSites((IPos1+1):IPos2),',')-1
          WRITE(SLengthofNumberString,*) LEN(SAtomicSites((IPos1+1):(IPos1+IPos))) 
          WRITE(SFormatString,*) "(I"//TRIM(ADJUSTL(SLengthofNumberString))//")"
          READ(SAtomicSites((IPos1+1):(IPos1+IPos)),FMT=SFormatString) IAtomsToRefine(ind)
          IPos1 = IPos1 + IPos + 1 
        ELSE
          WRITE(SLengthofNumberString,*) LEN(SAtomicSites((IPos1+1):(IPos2-1))) 
          WRITE(SFormatString,*) "(I"//TRIM(ADJUSTL(SLengthofNumberString))//")"
          READ(SAtomicSites((IPos1+1):(IPos2-1)),FMT=SFormatString) IAtomsToRefine(ind)
        END IF
      END DO
    END IF
    CALL message (LM, "Refining atoms ", IAtomsToRefine )
    
  END SUBROUTINE DetermineRefineableAtomicSites

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

    !?? interesting exercise - look into making pure JR

    ! e.g. CALL ThreeDimVectorReadIn(SIncidentBeamDirection,'[',']',RZDirC)
    ! Current format has interesting bracket inputs

    USE MyNumbers
   
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
    SFormattedVectorX = & 
          TRIM(ADJUSTL(SUnformattedVector(IOpenBracketPosition+1:IFirstCommaPosition-1)))
    SFormattedVectorY = &
          TRIM(ADJUSTL(SUnformattedVector(IFirstCommaPosition+1:ILastCommaPosition-1)))
    SFormattedVectorZ = &
          TRIM(ADJUSTL(SUnformattedVector(ILastCommaPosition+1:ICloseBracketPosition-1)))

    ! Read each string component into its associated position in real variable RFormattedVector
    READ(SFormattedVectorX,*) RFormattedVector(1)
    READ(SFormattedVectorY,*) RFormattedVector(2)
    READ(SFormattedVectorZ,*) RFormattedVector(3)

  END SUBROUTINE ThreeDimVectorReadIn

END MODULE read_files_mod
