

!>
!! Procedure-description: Write out the sample input file, when none provided
!!
!! Major-Authors: Keith Evans (2014)
!! 
SUBROUTINE WriteOutInputFile (IErr)

  !?? low priority
  !?? need to test this
  
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








  !Weighting Coefficients for figure of merit combination
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RWeightingCoefficients

!--------------------------------------------------------------------
! set up weighting coefficients
!--------------------------------------------------------------------

! Weighting parameter
ALLOCATE(RWeightingCoefficients(INoOfLacbedPatterns),STAT=IErr) 
SELECT CASE (IWeightingFLAG)
CASE(0) ! uniform weighting
   RWeightingCoefficients = ONE
CASE(1) ! smaller g's more important
   DO ind = 1,INoOfLacbedPatterns
      !?? NB untested, does RgPoolMag(ind)match output reflection (ind)?
      RWeightingCoefficients(ind) = RgPoolMag(ind)/MAXVAL(RgPoolMag)
   END DO
CASE(2) ! larger g's more important
   DO ind = 1,INoOfLacbedPatterns
      RWeightingCoefficients(ind) = MAXVAL(RgPoolMag)/RgPoolMag(ind)
   END DO
END SELECT









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
  
  IMPLICIT NONE

  INTEGER(IKIND):: IErr

  !Calling all functions which read felix.inp, felix.sca and felix.cif
  !(felix.hkl too depending on user preference)
  !ensure all input files are in working directory
  
  !felix.inp
  CALL ReadInpFile(IErr)
  IF(l_alert(IErr,"ReadInput","ReadInpFile")) RETURN  

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
!! Procedure-description: Convert vector movements into atomic coordinates
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!! 
SUBROUTINE ConvertVectorMovementsIntoAtomicCoordinates(IVariableID,RIndependentVariable,IErr)
  !RB this is now redundant, moved up to Update Variables
  USE MyNumbers

  USE SConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,jnd,IVariableID,IVectorID,IAtomID
  REAL(RKIND),DIMENSION(INoOfVariables),INTENT(IN) :: RIndependentVariable

!!$  Use IVariableID to determine which vector is being applied (IVectorID)
  IVectorID = IIterativeVariableUniqueIDs(IVariableID,3)
!!$  Use IVectorID to determine which atomic coordinate the vector is to be applied to (IAtomID)
  IAtomID = IAllowedVectorIDs(IVectorID)
!!$  Use IAtomID to applied the IVectodID Vector to the IAtomID atomic coordinate
  RBasisAtomPosition(IAtomID,:) = RBasisAtomPosition(IAtomID,:) + &
       RIndependentVariable(IVariableID)*RAllowedVectors(IVectorID,:)

END SUBROUTINE ConvertVectorMovementsIntoAtomicCoordinates

!------------------------------------------------------------------------







!>
!! Procedure-description: Update structure factors
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!! 
SUBROUTINE UpdateStructureFactors(RIndependentVariable,IErr)

  USE MyNumbers

  USE SConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels
  USE terminal_output
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,jnd
  REAL(RKIND),DIMENSION(INoOfVariables),INTENT(IN) :: RIndependentVariable
  CHARACTER*200 :: SPrintString

  !NB these are Ug's without absorption
  jnd=1
  DO ind = 1+IUgOffset,INoofUgs+IUgOffset!=== temp changes so real part only***
    IF ( (ABS(REAL(CUniqueUg(ind),RKIND)).GE.RTolerance).AND.&!===
        (ABS(AIMAG(CUniqueUg(ind))).GE.RTolerance)) THEN!use both real and imag parts!===
      CUniqueUg(ind)=CMPLX(RIndependentVariable(jnd),RIndependentVariable(jnd+1))!===
      jnd=jnd+2!===
    ELSEIF ( ABS(AIMAG(CUniqueUg(ind))).LT.RTolerance ) THEN!use only real part!===
      CUniqueUg(ind)=CMPLX(RIndependentVariable(jnd),ZERO)!===
      !===      CUniqueUg(ind)=CMPLX(RIndependentVariable(jnd),AIMAG(CUniqueUg(ind)))!===replacement line, delete to revert
      jnd=jnd+1
    ELSEIF ( ABS(REAL(CUniqueUg(ind),RKIND)).LT.RTolerance ) THEN!===use only imag part
      CUniqueUg(ind)=CMPLX(ZERO,RIndependentVariable(jnd))!===
      jnd=jnd+1!===
    ELSE!should never happen!===
      !todo - warning grouped with errors?
      CALL message(LS, "Warning - zero structure factor! At element = ",ind)
      CALL message(LS, "     CUniqueUg vector element value = ", CUniqueUg(IEquivalentUgKey(ind)))!===
    END IF!===
  END DO
  RAbsorptionPercentage = RIndependentVariable(jnd)!===

END SUBROUTINE UpdateStructureFactors

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













!>
!! Procedure-description: Sets up a pseudo random sequence and selects a number
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
REAL(RKIND) FUNCTION RANDOMNUMBER(IRequestedNumber,IErr)

  USE MyNumbers
  
  USE SConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,values(1:8), k,IRequestedNumber
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: seed
  REAL(RKIND),DIMENSION(IRequestedNumber) :: RRandomNumberSequence
  
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

!-------------------------------------------------------------------












!>
!! Procedure-description: Checks that vector movement applied by the simplex
!! initialisation does not move an atom out fo the unit cell, and if it does
!! the atom is moved back into the unit cell on the opposite side as if the
!! atom had moved from one unit cell into the neighbouring one
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE OutofUnitCellCheck(IVariableID,RProposedMovement,RCorrectedMovement,IErr)

  !?? Can't this just be done in one line with MODULO?
  USE MyNumbers
  
  USE SConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,IErr,IVariableID,IAtomID,IVectorID
  REAL(RKIND),DIMENSION(ITHREE) :: RProposedAtomicCoordinate,RDummyMovement
  REAL(RKIND),INTENT(IN) :: RProposedMovement
  REAL(RKIND),INTENT(OUT) :: RCorrectedMovement
  
  IVectorID = IIterativeVariableUniqueIDs(IVariableID,3)
  
  IAtomID = IAllowedVectorIDs(IVectorID)
 
  RProposedAtomicCoordinate(:) = RBasisAtomPosition(IAtomID,:) + &
       RProposedMovement*RAllowedVectors(IVectorID,:)

  RDummyMovement = RProposedMovement

  IF(ANY(RProposedAtomicCoordinate.GT.ONE).OR.ANY(RProposedAtomicCoordinate.LT.ZERO)) THEN
     DO ind = 1,ITHREE
        IF (RProposedAtomicCoordinate(ind).GT.ONE) THEN
           RDummyMovement(ind) = (ONE-RBasisAtomPosition(IAtomID,ind))/RAllowedVectors(IVectorID,ind)
        ELSEIF(RProposedAtomicCoordinate(ind).LT.ZERO) THEN
           RDummyMovement(ind) = (-RBasisAtomPosition(IAtomID,ind))/RAllowedVectors(IVectorID,ind)
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











SUBROUTINE OpenData(IChOutWrite, prefix, surname, IErr)

  USE MyNumbers
  USE WriteToScreen

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  CHARACTER*27 :: surname, surnamelength
  CHARACTER*2 :: prefix,postfix
  INTEGER(IKIND) :: IChOutWrite, IErr
  CHARACTER*34 :: filename
  INTEGER(IKIND) :: index

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
  IErr= 1
  RETURN
  
END SUBROUTINE OpenData













SUBROUTINE OpenImageForReadIn(IErr,filename)
!this is redundant
  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(KIND=IKIND) IChOutWrite, IErr

  CHARACTER*34 filename

  IF((IWriteFLAG.GE.10.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"OpenImageForReadIn()"
  END IF
  IF((IWriteFLAG.GE.10.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
    PRINT*,filename
  END IF

  OPEN(UNIT= IChInImage, ERR= 10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),FORM='UNFORMATTED',&
       ACCESS='DIRECT',IOSTAT=Ierr,RECL=2*IPixelCount*8)
  RETURN

  ! error in OPEN detected
10 PRINT*,"OpenImageForReadIn(): ERR in OPEN()"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenImageForReadIn









SUBROUTINE ReadImageForRefinement(IErr)
!this is redundant
  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind

  IF((IWriteFLAG.GE.10.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ReadImageForRefinement()"
  END IF

  DO ind=1,2*IPixelCount
     READ(IChInImage,rec=ind,ERR=10) RImageIn(ind,:)
  END DO
  RETURN

10 IErr=1
  PRINT*,"ReadImageForRefinement (", my_rank, ") error in READ()",IErr
  RETURN
  
END SUBROUTINE ReadImageForRefinement














SUBROUTINE WriteCif(IErr)

  USE MyNumbers
  USE WriteToScreen

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INCLUDE       'ciftbx-f90.cmn'

  INTEGER(IKIND) :: IErr
  REAL(RKIND),DIMENSION(SIZE(RBasisAtomPosition,DIM=1),SIZE(RBasisAtomPosition,DIM=2)) :: &
       ROutputData
  REAL(RKIND),DIMENSION(2,3) :: RUnitCellParameters
  LOGICAL :: f1

  IF(.NOT.dict_('cif_core.dic','valid')) THEN
     PRINT*,"Requested Core Dictionary not Present"
  END IF

  IF(.NOT.pfile_('felixoutput.cif')) THEN
     PRINT*,"Cif file already exists"
  END IF

  f1 = pdata_('DataBlock') !Open a Data Block

  call close_

END SUBROUTINE WriteCif










SUBROUTINE MontageSetup(RMontageImages,RIndividualReflectionImages,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  %
!!$  %    Creates Montage Images
!!$  %
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  USE WriteToScreen
  USE MyNumbers
  USE IConst
  
  USE MPI
  USE MyMPI
  
  USE IPara; USE RPara
  
  IMPLICIT NONE

  INTEGER(IKIND):: IErr,IThicknessIndex,knd,ind,jnd
  REAL(RKIND),DIMENSION(MAXVAL(IImageSizeXY),&
       MAXVAL(IImageSizeXY),IThicknessCount):: RMontageImages
  REAL(RKIND),DIMENSION(INoOfLacbedPatterns,IThicknessCount,IPixelTotal):: &
       RIndividualReflectionImages

     DO IThicknessIndex =1,IThicknessCount
        DO knd = 1,IPixelTotal
           jnd = IPixelLocations(knd,1)
           ind = IPixelLocations(knd,2)
           CALL MontageInitialisation(ind,jnd,IThicknessIndex,RMontageImages,&
                RIndividualReflectionImages(:,IThicknessIndex,knd),IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"MontageSetup(", my_rank, ") error ", IErr, &
                   " in MakeMontagePixel"
              RETURN
           ENDIF
        END DO
     END DO

END SUBROUTINE MontageSetup







SUBROUTINE MicroscopySettings( IErr )
!This is now redundant
  USE MyNumbers
  USE WriteToScreen
  
  USE SConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  REAL(RKIND)::norm,dummy,ROneThousand
  INTEGER :: IErr

  ROneThousand = 1000.0_RKIND

  !Electron velocity in metres per second
  RElectronVelocity= RSpeedOfLight*SQRT(ONE-((RElectronMass*RSpeedOfLight**2)/ &
       (RElectronCharge*RAcceleratingVoltage*ROneThousand + &
        RElectronMass*RSpeedOfLight**2) )**2 )
  
  !Electron wavelength in Angstroms
  RElectronWaveLength= RPlanckConstant / &
       ( SQRT(TWO*RElectronMass*RElectronCharge*RAcceleratingVoltage*ROneThousand) * &
         SQRT(ONE + (RElectronCharge*RAcceleratingVoltage*ROneThousand) / &
          (TWO*RElectronMass*RSpeedOfLight**2) ))* RAngstromConversion
  
  !(NB k=2pi/lambda and exp(i*k.r), physics convention)
  RElectronWaveVectorMagnitude=TWOPI/RElectronWaveLength
  RRelativisticCorrection= ONE/SQRT( ONE - (RElectronVelocity/RSpeedOfLight)**2 )
  RRelativisticMass= RRelativisticCorrection*RElectronMass
     
  RETURN

END SUBROUTINE MicroscopySettings









SUBROUTINE ReadScaFile( IErr )
  !This is now redundant
  !completely pointless subroutine that just calls another subroutine
  USE MyNumbers
  USE WriteToScreen

  USE SConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE


  INTEGER IErr, ILine

  CALL ScatteringFactors(IScatterFactorMethodFLAG,IErr)

END SUBROUTINE ReadScaFile







SUBROUTINE RecoverSavedSimplex(RSimplexVariable,RSimplexFoM,RStandardDeviation,RMean,Iter,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % This Subroutine reads the fr-simplex.txt file from a previous
!!$  % refinement run, and recreates the simplex volume and tolerances
!!$  % allowing for the continuation of a previous refinement
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE SConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind,Iter
  REAL(RKIND),DIMENSION(INoOfVariables+1,INoOfVariables) :: RSimplexVariable
  REAL(RKIND),DIMENSION(INoOfVariables+1) :: RSimplexFoM
  REAL(RKIND) :: RStandardDeviation,RMean
  CHARACTER*200 :: CSizeofData,SFormatString,filename

  WRITE(filename,*) "fr-Simplex.txt"

  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
        FILE=TRIM(ADJUSTL(filename)))
  
  WRITE(CSizeofData,*) INoOfVariables+1
  WRITE(SFormatString,*) "("//TRIM(ADJUSTL(CSizeofData))//"(1F6.3,1X),A1)"

  DO ind = 1,(INoOfVariables+1)
     READ(IChOutSimplex,FMT=SFormatString) RSimplexVariable(ind,:),RSimplexFoM(ind)
  END DO
    
  READ(IChOutSimplex,FMT="(2(1F6.3,1X),I5.1,I5.1,A1)") RStandardDeviation,RMean,IStandardDeviationCalls,Iter

  CLOSE(IChOutSimplex)

END SUBROUTINE RecoverSavedSimplex











SUBROUTINE RandomSequence(RRandomSequence,IRandomSequenceLength,ISeedModifier,IErr)
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Sets up a pseudo random sequence and selects a number

  USE MyNumbers
  
  USE SConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,Ivalues(1:8), k,IRandomSequenceLength,ISeedModifier
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: seed
  REAL(RKIND),DIMENSION(IRandomSequenceLength) :: RRandomSequence
  
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










SUBROUTINE InitialiseAtomicVectorMagnitudes(IVariableID,RCorrectedMovement,IErr)
  
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Creates pseudo random movements of atoms using allowed vectors
!!$  % to initialise the simplex, proposed movements which exit the unit
!!$  $ cell are corrected to bring the atom back in on the opposite side
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE SConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,IVariableID
  REAL(RKIND) :: RNegativeMovement,RPositiveMovement,RCorrectedMovement,RANDOMNUMBER

  RNegativeMovement = RSimplexLengthScale*(-1.0_RKIND)
  RPositiveMovement = RSimplexLengthScale
!RB this check can be done in less lines than it takes to call the subroutine
  IF(RANDOMNUMBER(IVariableID,IErr).LT.0.5_RKIND) THEN
     CALL OutofUnitCellCheck(IVariableID,RNegativeMovement,RCorrectedMovement,IErr)
  ELSE
     CALL OutofUnitCellCheck(IVariableID,RPositiveMovement,RCorrectedMovement,IErr)
  END IF

END SUBROUTINE InitialiseAtomicVectorMagnitudes







SUBROUTINE RefinementVariableSetup(RIndependentVariable,IErr)
!This is redundant  
  USE MyNumbers
  
  USE SConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,ind,IVariableType
  REAL(RKIND),DIMENSION(INoOfVariables),INTENT(OUT) :: RIndependentVariable
  
  IF((IWriteFLAG.GE.10.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RefinementVariableSetup(",my_rank,")"
  END IF
  
!!$  Fill the Independent Value array with values

  DO ind = 1,INoOfVariables
     IVariableType = IIterativeVariableUniqueIDs(ind,2)
     SELECT CASE (IVariableType)
     CASE(1)
	    !Structure factor refinement, define in SymmetryRelatedStructureFactorDetermination
    CASE(2)
        RIndependentVariable(ind) = RAllowedVectorMagnitudes(IIterativeVariableUniqueIDs(ind,3))
     CASE(3)
        RIndependentVariable(ind) = RBasisOccupancy(IIterativeVariableUniqueIDs(ind,3))
     CASE(4)
        RIndependentVariable(ind) = RBasisIsoDW(IIterativeVariableUniqueIDs(ind,3))
     CASE(5)
        RIndependentVariable(ind) = RAnisotropicDebyeWallerFactorTensor(&
             IIterativeVariableUniqueIDs(ind,3),&
             IIterativeVariableUniqueIDs(ind,4),&
             IIterativeVariableUniqueIDs(ind,5))
     CASE(6)
        SELECT CASE(IIterativeVariableUniqueIDs(ind,3))
        CASE(1)
           RIndependentVariable(ind) = RLengthX
        CASE(2)
           RIndependentVariable(ind) = RLengthY
        CASE(3)
           RIndependentVariable(ind) = RLengthZ
        END SELECT
     CASE(7)
        SELECT CASE(IIterativeVariableUniqueIDs(ind,3))
        CASE(1)
           RIndependentVariable(ind) = RAlpha
        CASE(2)
           RIndependentVariable(ind) = RBeta
        CASE(3)
           RIndependentVariable(ind) = RGamma
        END SELECT
     CASE(8)
        RIndependentVariable(ind) = RConvergenceAngle
     CASE(9)
        RIndependentVariable(ind) = RAbsorptionPercentage
     CASE(10)
        RIndependentVariable(ind) = RAcceleratingVoltage
     CASE(11)
        RIndependentVariable(ind) = RRSoSScalingFactor
     CASE(12)
	    !Structure factor refinement, define in SymmetryRelatedStructureFactorDetermination
     END SELECT
  END DO

END SUBROUTINE RefinementVariableSetup






REAL(RKIND) FUNCTION RStandardError(RStandardDeviation,RMean,IErr)

  USE MyNumbers

  USE SConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr
  REAL(RKIND),INTENT(INOUT) :: RStandardDeviation,RMean

  IF (IStandardDeviationCalls.GT.1) THEN
     RMean = (RMean*REAL(IStandardDeviationCalls,RKIND) + &
          RFigureofMerit)/REAL(IStandardDeviationCalls+1,RKIND)
     RStandardDeviation = SQRT(&
          ((REAL(IStandardDeviationCalls,RKIND)*RStandardDeviation**2)+&
          (RFigureofMerit-RMean)**2)/ &
          REAL(IStandardDeviationCalls+1,RKIND))   
  ELSE
     RMean = RFigureofMerit
     RStandardDeviation = ZERO
  END IF
  RStandardError = RStandardDeviation/SQRT(REAL(IStandardDeviationCalls+1,RKIND))
  IStandardDeviationCalls = IStandardDeviationCalls + 1

END FUNCTION  RStandardError










SUBROUTINE ExperimentalSetup (IErr)

  USE WriteToScreen
  USE MyNumbers
  USE IConst
  
  USE IPara; USE RPara; USE SPara; USE CPara
 
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr

  !--------------------------------------------------------------------
  ! crystallography settings
  CALL CrystallographyInitialisation( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"ExperimentalSetup(",my_rank,")error",IErr,"in CrystallographyInitialisation()"
     RETURN
  ENDIF

  !--------------------------------------------------------------------
  ! diffraction initialization
  !--------------------------------------------------------------------
  CALL DiffractionPatternInitialisation( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"ExperimentalSetup(",my_rank,")error",IErr,"in DiffractionPatternInitialisation()"
     RETURN
  ENDIF

END SUBROUTINE ExperimentalSetup












SUBROUTINE DiffractionPatternInitialisation

  USE WriteToScreen
  USE MyNumbers
  USE IConst

  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr

  CALL ReflectionDetermination (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternInitialisation(",my_rank,")error", IErr, &
          "in ReflectionDetermination()"
     RETURN
  ENDIF

  CALL SpecificReflectionDetermination (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternInitialisation(", my_rank, ") error", IErr, &
          "in SpecificReflectionDetermination()"
     RETURN
  ENDIF

  CALL DiffractionPatternCalculation (IErr)
   IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternInitialisation(", my_rank, ") error", IErr, &
          "in DiffractionPatternCalculation()"
     RETURN
  ENDIF

END SUBROUTINE DiffractionPatternInitialisation












SUBROUTINE NewHKLMake(Ihklmax,Rhkl0Vec,RHOLZAcceptanceAngle,IErr)
  
  USE MyNumbers
  USE WriteToScreen
  
  USE SConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  USE MPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr, Ihklmax,ind,jnd,knd,INhkl
  REAL(RKIND) :: RHOLZAcceptanceAngle
  REAL(RKIND), DIMENSION(ITHREE) :: Rhkl0Vec,RhklDummyUnitVec,RhklDummyVec,Rhkl0UnitVec

  INhkl = 0
  
  Rhkl0UnitVec= Rhkl0Vec/SQRT(DOT_PRODUCT(REAL(Rhkl0Vec,RKIND),REAL(Rhkl0Vec,RKIND)))

!First count the number of reflections in the acceptance angle
!??? doesn't this only work for cubic systems?  Where are the magnitudes of the reciprocal lattice vectors?
  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1          
           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)          
           IF( (ind.NE.0).OR.(jnd.NE.0).OR.(knd.NE.0) ) THEN
              RhklDummyUnitVec= RhklDummyVec / &
                   SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))
           ELSE
              RhklDummyUnitVec = RhklDummyVec
           END IF
           
           SELECT CASE(SSpaceGroupName)
              
           CASE("F") !Face Centred
              IF(((ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),TWO))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),TWO))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),TWO))).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),TWO))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),TWO))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),TWO))).GT.TINY))) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1
                 END IF
              END IF
			  
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("P")! Primitive
              IF(IHolzFLAG.EQ.0) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                   .LE.SIN(RHOLZAcceptanceAngle)) THEN
                 INhkl = INhkl +1       
              END IF
           CASE DEFAULT
              PRINT*,"HKLMake: unknown space group", SSpaceGroupName, ", aborting"
              IErr=1
              RETURN
           END SELECT

        END DO
     END DO
  END DO

!RB now allocate the hkl list...  
  ALLOCATE(Rhkl(INhkl,ITHREE),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"hklMake(",my_rank,")error allocating Rhkl"
     RETURN
  END IF

!RB ...and calculate it all again, filling Rhkl  
  INhkl = 0

  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1

           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)

           IF(ind.NE.0.AND.jnd.NE.0.AND.knd.NE.0) THEN
              RhklDummyUnitVec= RhklDummyVec / &
                   SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))
           ELSE
              RhklDummyUnitVec = RhklDummyVec
           END IF
           
           SELECT CASE(SSpaceGroupName)
		   
           CASE("F") !Face Centred
              IF(((ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),TWO))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),TWO))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),TWO))).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),TWO))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),TWO))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),TWO))).GT.TINY))) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("P")! Primitive

		   IF(IHolzFLAG.EQ.0) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                    Rhkl(INhkl,:)= RhklDummyVec
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                 INhkl =  INhkl + 1
                 Rhkl(INhkl,:) = RhklDummyVec                 
              END IF
           CASE DEFAULT
              PRINT*,"HKLMake(): unknown space group", SSpaceGroupName, "--- aborting"
              IErr=1
              RETURN
           END SELECT

        END DO
     END DO
  END DO

END SUBROUTINE NewHKLmake


















SUBROUTINE DiffractionPatternCalculation (IErr)
  !RB This is now redundant
  USE MyNumbers
  USE WriteToScreen
  USE IConst
  
  USE IPara; USE RPara;
  
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,IErr
  CHARACTER*20 :: Sind
  
  DO ind =1,SIZE(Rhkl,DIM=1)
     RgDotNorm(ind) = DOT_PRODUCT(RgPool(ind,:),RNormDirM)
  END DO   
  
  ! smallest g is gmag(2) IF 000 beam is included !!!add error catch here
  RMinimumGMag = RgPoolMag(2)
  
  IF (nReflections.LT.INoOfLacbedPatterns) THEN
     nReflections = INoOfLacbedPatterns
  END IF
  
  ! resolution in k space
  RDeltaK = RMinimumGMag*RConvergenceAngle/REAL(IPixelCount,RKIND)

  RETURN

END SUBROUTINE DiffractionPatternCalculation















SUBROUTINE CrystallographyInitialisation( IErr )

  USE MyNumbers
  USE WriteToScreen
  
  USE SConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr, ind

  !Removed PseudoCubic translation

  CALL ReciprocalLattice(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
             " in ReciprocalLattice"
     RETURN
  ENDIF

  CALL UniqueAtomPositions(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
             " in UniqueAtomPositions "
     RETURN
  ENDIF

  !zz CALL CrystalUniqueFractionalAtomicPostitionsCalculation(IErr)
!zz    IF( IErr.NE.0 ) THEN
!zz      PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
!zz              " in CrystalUniqueFractionalAtomicPostitionsCalculation "
!zz      RETURN
!zz   ENDIF

END SUBROUTINE CrystallographyInitialisation








