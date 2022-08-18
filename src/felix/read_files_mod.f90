!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: :VERSION: RB_coord / 1.15 /
! Date:    :DATE: 16-01-2019
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build:   :BUILD: Mode F: test different lattice types" 
! Author:  :AUTHOR: r.beanland
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

    !?? this should contain clear input information and be well documented
    !?? github and other places should referance here

    USE MyNumbers
    USE message_mod

    ! global outputs, read from .inp
    USE IPARA, ONLY : IWriteFLAG, IImageFLAG, IScatterFactorMethodFLAG, IHolzFLAG, &
          IAbsorbFLAG, IByteSize, INhkl, &
          IMinStrongBeams, IMinWeakBeams, ISimFLAG, IRefineMode, &
          IWeightingFLAG, IRefineMethodFLAG, ICorrelationFLAG, IImageProcessingFLAG, &
          INoofUgs, IPrint, IPixelCount
    USE RPARA, ONLY : RDebyeWallerConstant, RAbsorptionPercentage, RConvergenceAngle, &
          RZDirC, RXDirC, RNormDirC, RAcceleratingVoltage, RAcceptanceAngle, &
          RInitialThickness, RFinalThickness, RDeltaThickness, RBlurRadius, &
          RSimplexLengthScale, RExitCriteria,RPrecision
    USE SPARA, ONLY : SPrintString
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
          SDirectionX, SIncidentBeamDirection, SNormalDirectionX

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
    ! IImageFLAG - select which type(s) of images to output 
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IImageFLAG
    ! IScatterFactorMethodFLAG
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IScatterFactorMethodFLAG
    CALL message ( LXL, dbg7, "IScatterFactorMethodFLAG=",IScatterFactorMethodFLAG )
    ! IHolzFLAG
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IHolzFLAG
    CALL message ( LXL, dbg3, "IHolzFLAG=",IHolzFLAG)
    ! IAbsorbFLAG
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IAbsorbFLAG
    CALL message ( LXL, dbg3, "IAbsorbFLAG=",IAbsorbFLAG)
    ! IByteSize
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IByteSize
    CALL message ( LXL, dbg3, "IByteSize=",IByteSize) ! depends on system, 8 for csc, 2 tinis

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
    ! INhkl
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) INhkl
    CALL message ( LXL, dbg3, "INhkl=",INhkl)
    ! IMinStrongBeams
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IMinStrongBeams
    CALL message ( LXL, dbg3, "IMinStrongBeams=",IMinStrongBeams)
    ! IMinWeakBeams
    ILine= ILine+1
    READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IMinWeakBeams
    CALL message ( LXL, dbg3, "IMinWeakBeams=",IMinWeakBeams)
    WRITE(SPrintString,FMT='(A19,I4,A19,I4,A13)') "Reflection pool of ",INhkl," with a minimum of ",IMinStrongBeams," strong beams"
    CALL message ( LS, SPrintString)

    !--------------------------------------------------------------------
    ! crystal settings
    !--------------------------------------------------------------------

    ! Two comment lines          
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ! RDebyeWallerConstant          ! default, if not specified in .cif
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RDebyeWallerConstant
    ! RAbsorptionPercentage         ! for proportional model of absorption
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RAbsorptionPercentage

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
    ! Two comment lines
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
    ! RInitialThickness
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RInitialThickness
    ! RFinalThickness
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RFinalThickness
    ! RDeltaThickness
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RDeltaThickness
    ! RPrecision - used for error calculations
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RPrecision 

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
      IF(IRefineMode(1) .EQ.1) CALL message( LS, "Refining Structure Factors, A")
      IF(IRefineMode(2) .EQ.1) CALL message( LS, "Refining Atomic Coordinates, B")
      IF(IRefineMode(3) .EQ.1) CALL message( LS, "Refining Occupancies, C")
      IF(IRefineMode(4) .EQ.1) CALL message( LS, "Refining Isotropic Debye Waller Factors, D")
      IF(IRefineMode(5) .EQ.1) CALL message( LS, "Refining Anisotropic Debye Waller Factors, E")
      IF(IRefineMode(6) .EQ.1) CALL message( LS, "Refining Lattice Parameters, F")
      IF(IRefineMode(7) .EQ.1) CALL message( LS, "Refining Lattice Angles, G")
      IF(IRefineMode(8) .EQ.1) CALL message( LS, "Refining Convergence Angle, H")
      IF(IRefineMode(9) .EQ.1) CALL message( LS, "Refining Accelerating Voltage, I")
      ! Error Check - user cannot request Ug refinement and anything else
      IF((IRefineMode(1).EQ.1).AND.SUM(IRefineMode).GT.1) THEN
        IErr = 1; IF(l_alert(IErr,"ReadInpFile",&
              "Structure factors must be refined separately")) RETURN
      END IF
    END IF

    ! IWeightingFLAG         
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IWeightingFLAG
    ! IRefineMethodFLAG
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IRefineMethodFLAG
    IF(ISimFLAG==0) THEN
      IF(IRefineMethodFLAG.EQ.1) CALL message( LS, "Refining by simplex")
      IF(IRefineMethodFLAG.EQ.2) CALL message( LS, "Refining by downhill (2-point) gradient")
      IF(IRefineMethodFLAG.EQ.3) CALL message( LS, "Refining by maximum (3-point) gradient")
      IF(IRefineMethodFLAG.EQ.4) CALL message( LS, "Refining by pairwise (2x2-point) gradient")
    END IF
   
    ! ICorrelationFLAG: 0=phase,1=sumSq,2=NormalisedCC,3=masked
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) ICorrelationFLAG
    ! IImageProcessingFLAG: 0=no processing,1=sqrt,2=log
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) IImageProcessingFLAG
    ! RBlurRadius
    ILine= ILine+1; READ(IChInp,'(27X,F18.9)',ERR=20,END=30) RBlurRadius
    ! INoofUgs
    ILine= ILine+1; READ(IChInp,'(27X,I15.1)',ERR=20,END=30) INoofUgs
    ! SAtomicSites
    ILine=ILine+1; READ(IChInp,FMT='(A)',ERR=20,END=30) SAtomicSites
    CALL DetermineRefineableAtomicSites(SAtomicSites,IErr)
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
    USE IPARA, ONLY : INoOfLacbedPatterns,IOutputReflections ! allocated here
    ! global inputs
    USE IChannels, ONLY : IChInp

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(OUT) :: IErr
    INTEGER(IKIND) :: ILine, h, k, l, ind, IPos1, IPos2
    CHARACTER(200) :: dummy1, dummy2

    OPEN(Unit = IChInp,FILE="felix.hkl",STATUS='OLD',ERR=10)
    ILine = 0

    ! count the number of lines in felix.hkl:
    ! this is the number of reflections to output, INoOfLacbedPatterns
    DO
      READ(UNIT= IChInp, END=5, FMT='(a)') dummy1
      ILine = ILine+1
    ENDDO
5   INoOfLacbedPatterns = ILine
    CALL message ( LXL, dbg7, "Number of experimental images to load = ", INoOfLacbedPatterns)

    ALLOCATE(RInputHKLs(INoOfLacbedPatterns,ITHREE),STAT=IErr)
    IF(l_alert(IErr,"ReadHklFile","allocate RInputHKLs")) RETURN
    ALLOCATE(IOutputReflections(INoOfLacbedPatterns),STAT=IErr) !?? should we allocate here JR
    IF(l_alert(IErr,"ReadHklFile","allocate IOutputReflections")) RETURN

    ! read in the hkls
    REWIND(UNIT=IChInp) ! goes to beginning of felix.hkl file
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
    IErr=1
    CALL message( LL, "felix.hkl not found")
    RETURN

  END SUBROUTINE ReadHklFile

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: This subroutine reads in the input experimental images.
  !! It inquires whether there are .img or .dm3 files in LR_NxN/, DM3/ or directly in the sample folder.
  !! It then reads in the images files expecting them to match the necessary LACBED patterns.
  !! With the .dm3 files, it writes out the equivalent .img files to be used next time.
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016), Jacob Richardson (2017)
  !!
  SUBROUTINE ReadExperimentalImages(IErr)

    USE MyNumbers
    USE message_mod

    USE read_dm3_mod

    ! global outputs
    USE RPARA, ONLY : RImageExpi

    ! global inputs
    USE IChannels, ONLY : IChInImage, IChOutWIImage
    USE SPARA, ONLY : SChemicalFormula
    USE IPARA, ONLY : INoOfLacbedPatterns,IPixelCount,IByteSize,ILN
    USE RPARA, ONLY : RInputHKLs
    USE SPARA, ONLY : SPrintString
    
    IMPLICIT NONE

    INTEGER(IKIND),INTENT(OUT) :: IErr
    INTEGER(IKIND) :: ind, jnd, INegError = 0, IPixelArray(2), IFileTypeID, IFound
    INTEGER(8) :: IFileSize
    CHARACTER(100) :: SFilename,SPath,SImageExtension,SFilePath
    LOGICAL :: LFileExist
    REAL(4),ALLOCATABLE :: RImage4ByteFloatDM3(:,:)

    ! for IByteSize: 2bytes=64-bit input file (NB tinis specifies in bytes, not bits)
    !?? JR Is it 8bytes=64-bit and not 2bytes=64-bit?
    ! iteratively INQUIRE each possible location for +0+0+0 .bin or .dm3 image
    IFound=0!flag to say if we found a readable a file
    DO IFileTypeID=1,4
      SELECT CASE(IFileTypeID)
        CASE(1) ! .img in LR_NxN/
          WRITE(SPath,'(A,I0,A,I0,A)') 'LR_',2*IPixelCount,'x',2*IPixelCount,'/'
          WRITE(SImageExtension,'(A)') '.img'
          ! NB pixel size read from felix.inp and this is expected to match pixels in foldername
        CASE(2) ! .dm3 in DM3/
          SPath='DM3/'
          WRITE(SImageExtension,'(A)') '.dm3'
        CASE(3) ! .img directly in sample directory
          SPath=''
          WRITE(SImageExtension,'(A)') '.img'
        CASE(4) ! .dm3 directly in sample directory
          SPath=''
          WRITE(SImageExtension,'(A)') '.dm3'
      END SELECT
      WRITE(SFilePath ,'(A,A,A,A)') TRIM(SPath),SChemicalFormula(1:ILN),'_+0+0+0',TRIM(SImageExtension)

      ! check if corresponding _+0+0+0.img or _+0+0+0.dm3 image exists
      INQUIRE(FILE=TRIM(SFilePath) ,EXIST=LFileExist)
      IF(LFileExist) THEN
        IFound=1
        !IF (my_rank.EQ.0) PRINT*,"Found experimental image with filepath ",TRIM(SFilePath)
        CALL message(LM, "Found initial experimental image with filepath ",TRIM(SFilePath) )
        EXIT
      ELSE
        !IF (my_rank.EQ.0) PRINT*,"Did not find experimental image with filepath ",TRIM(SFilePath)
        CALL message(LM, "Did not find initial experimental image with filepath ",TRIM(SFilePath) )
      END IF    
    END DO
    ! NB once a file is found, the above do-loop is exited and the variables IFileTypeID, SFilePath and
    ! SPath will have the correct values to continue working with them.

    IF (IFound.EQ.0) THEN
      CALL message(LS, "Did not find experimental image at ",TRIM(SFilePath) )
      IErr = 1
      RETURN
    END IF

    ! if .dm3 allocate raw 4-byte float image matrix
    IF(IFileTypeID.EQ.2.OR.IFileTypeID.EQ.4) ALLOCATE(RImage4ByteFloatDM3(2*IPixelCount,2*IPixelCount),STAT=IErr)
    IF(l_alert(IErr,"ReadExperimentalImages","allocate RImage4ByteFloatDM3")) RETURN

    ! Read in expected image for each LacbedPattern
    DO ind = 1,INoOfLacbedPatterns
      WRITE(SFilename,'(A,A,SP,3(I0),A)') SChemicalFormula(1:ILN),"_",&
            NINT(RInputHKLs(ind,1:3)), TRIM(SImageExtension)
      SFilePath  = TRIM(SPath)//SFilename
      CALL message(LL, dbg7, "SFilename = ", SFilePath )

      ! do corresponding read-in process for .img or .dm3
      SELECT CASE(IFileTypeID)
        CASE(1,3) ! .img
          OPEN(UNIT=IChInImage, STATUS= 'UNKNOWN', FILE=TRIM(SFilePath), &
                FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=2*IPixelCount*IByteSize)
          IF(l_alert(IErr,"ReadExperimentalImages",&
                "OPEN() an experimental image, SFilePath="//TRIM(SFilePath)//&
                ', check input HKLs in felix.hkl')) RETURN
          DO jnd=1,2*IPixelCount
            READ(IChInImage,rec=jnd,IOSTAT=IErr) RImageExpi(jnd,:,ind)
            IF(l_alert(IErr,"ReadExperimentalImages",&
                  "READ() an experimental image, SFilePath="//TRIM(SFilePath)//&
                  ', check input HKLs in felix.hkl')) RETURN
          END DO
          CLOSE(IChInImage,IOSTAT=IErr)
          IF(l_alert(IErr,"ReadExperimentalImages",&
                "CLOSE() an experimental image, SFilePath="//TRIM(SFilePath))) RETURN

        CASE(2,4) ! .dm3
          ! read in .dm3
          !IF (my_rank.EQ.0) PRINT*,"Reading ",TRIM(SFilePath)
          CALL ReadDM3TagsAndImage( SFilePath, 2*IPixelCount, 2*IPixelCount, IErr, RImage4ByteFloatDM3 )
          IF(l_alert(IErr,"ReadExperimentalImages",&
                "ReadDM3TagsAndImage() where SFilePath="//TRIM(SFilePath))) RETURN
          RImageExpi(:,:,ind) = REAL( RImage4ByteFloatDM3, RKIND )
      END SELECT
    END DO

    ! if .dm3 deallocate raw 4-byte float image matrix
    IF(IFileTypeID.EQ.2.OR.IFileTypeID.EQ.4) THEN
      DEALLOCATE(RImage4ByteFloatDM3,STAT=IErr)
      IF(l_alert(IErr,"ReadExperimentalImages","deallocate RImage4ByteFloatDM3")) RETURN
    END IF

    ! If this is reached, subroutine has not exited early with error and hence
    ! images have been read-in correctly
    WRITE(SPrintString,'(I0,A)') INoOfLacbedPatterns,' experimental images successfully loaded'
    CALL message(LS,TRIM(ADJUSTL(SPrintString)))

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

    ! e.g. CALL ThreeDimVectorReadIn(SIncidentBeamDirection,'[',']',RZDirC)

    USE MyNumbers
   
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: SUnformattedVector,SOpenBracketDummy,SCloseBracketDummy
    REAL(RKIND),INTENT(OUT),DIMENSION(3) :: RFormattedVector
    CHARACTER(1) :: SComma=',',SOpenBracket,SCloseBracket
    CHARACTER(100) :: SFormattedVectorX,SFormattedVectorY,SFormattedVectorZ   
    LOGICAL :: LBACK=.TRUE.   
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
