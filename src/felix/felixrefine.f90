!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: 2.0
! Date: 31-08-2022
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build: cRED 
! Author:  r.beanland@warwick.ac.uk
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

!>
!! felixrefine
!!
PROGRAM Felixrefine
  
  USE MyNumbers
  USE message_mod
  USE MPI
  USE MyMPI
  USE read_files_mod
  USE read_cif_mod
  USE set_scatter_factors_mod
  USE setup_reflections_mod
  USE setup_space_group_mod
  USE crystallography_mod
  USE ug_matrix_mod
  USE bloch_mod
  USE write_output_mod

  USE IConst; USE RConst; USE SConst
  USE IPara;  USE RPara;  USE CPara; USE SPara;
!  USE BlochPara 
  USE IChannels

  ! local variable definitions
  IMPLICIT NONE
 
  INTEGER(IKIND) :: IErr,ind,jnd,knd,lnd,mnd,nnd,Iter,ICutOff,IHOLZgPoolMag,&
        IBSMaxLocGVecAmp,ILaueLevel,INumTotalReflections,ITotalLaueZoneLevel,&
        IPrintFLAG,ICycle,INumInitReflections,IZerothLaueZoneLevel,xnd,&
        INumFinalReflections,IThicknessIndex,IVariableType,IArrayIndex,&
        IAnisotropicDebyeWallerFactorElementNo,IStartTime,IStartTime2
  INTEGER(4) :: IErr4
  REAL(RKIND) :: RGlimit,RLaueZoneGz,RMaxGMag,RKn,RThickness,&
        RScale,RMaxUgStep,Rdx,RStandardDeviation,RMean,RGzUnitVec,&
        RMinLaueZoneValue,Rdf,RLastFit,RBestFit,RMaxLaueZoneValue,&
        RMaxAcceptanceGVecMag,RandomSign,RLaueZoneElectronWaveVectorMag,&
        RvarMin,RfitMin,RFit0,Rconvex,Rtest,Rplus,Rminus,RdeltaX,RdeltaY
  REAL(RKIND),DIMENSION(ITHREE) :: RXDirOn,RZDirOn

  CHARACTER(40) :: my_rank_string
  CHARACTER(20) :: h,k,l
  CHARACTER(14) :: Sest

  !--------------------------------------------------------------------
  ! startup
  !--------------------------------------------------------------------

  ! initialise constants
  CALL Init_Numbers ! constants for calculations
  CALL InitialiseMessage ! constants required for formatted terminal output
  IErr=0

  ! MPI initialization
  CALL MPI_Init(IErr4) 
  IF(l_alert(INT(REAL(IErr4)),"felixrefine","MPI_Init")) CALL abort
  CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,IErr4) ! get rank of the current process
  IF(l_alert(INT(REAL(IErr4)),"felixrefine","MPI_Comm_rank")) CALL abort
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,IErr4) ! get size of the current communicator
  IF(l_alert(INT(REAL(IErr4)),"felixrefine","MPI_Comm_size")) CALL abort

  ! startup terminal output
  CALL message(LS,"-----------------------------------------------------------------")
  CALL message(LS,"felixrefine: ", RStr)
  CALL message(LS,"             ", DStr)
  CALL message(LS,"             ", AStr)
  CALL message(LS,"-----------------------------------------------------------------")
  CALL message(LS,"total number of MPI ranks ", p, ", screen messages via rank", my_rank)
  CALL message(LS,"-----------------------------------------------------------------")

  ! timing setup
  CALL SYSTEM_CLOCK( IStartTime,IClockRate )

  !--------------------------------------------------------------------
  ! input section 
  !--------------------------------------------------------------------

  CALL read_cif(IErr) ! felix.cif ! some allocations are here
  IF(l_alert(IErr,"felixrefine","ReadCif")) CALL abort

  CALL ReadHklFile(IErr) ! the list of hkl's to input/output
  IF(l_alert(IErr,"felixrefine","ReadHklFile")) CALL abort
  !--------------------------------------------------------------------
  ! allocations for arrays to track frame simulations
  ALLOCATE(IhklsFrame(INoOfHKLsAll),STAT=IErr) ! Legacy list
  IF(l_alert(IErr,"felixrefine","allocate IhklsFrame")) CALL abort
  ALLOCATE(IhklsAll(INoOfHKLsAll),STAT=IErr)! List for full sim
  IF(l_alert(IErr,"felixrefine","allocate IhklsAll")) CALL abort
  ALLOCATE(ILiveList(INoOfHKLsAll),STAT=IErr)! List of current outputs
  IF(l_alert(IErr,"felixrefine","allocate IhklsAll")) CALL abort
  ILiveList = 0 !see write_outputs for all meanings of this flag

  CALL ReadInpFile(IErr) ! felix.inp
  IF(l_alert(IErr,"felixrefine","ReadInpFile")) CALL abort
  CALL SetMessageMode( IWriteFLAG, IErr )
  IF(l_alert(IErr,"felixrefine","set_message_mod_mode")) CALL abort

  !--------------------------------------------------------------------
  ! allocations, now we know the size of the beam pool
  ! RgPool is a list of g-vectors in the microscope ref frame,
  ! units of 1/A (NB exp(-i*q.r),  physics negative convention)
  ALLOCATE(RgPool(INhkl,ITHREE),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RgPool")) CALL abort
  ! g-vector magnitudes
  ! in reciprocal Angstrom units, in the Microscope reference frame
  ALLOCATE(RgPoolMag(INhkl),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RgPoolMag")) CALL abort
  ! g-vector components parallel to the surface unit normal
  ALLOCATE(RgDotNorm(INhkl),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RgDotNorm")) CALL abort
  ! Matrix of 2pi*g-vectors that corresponds to the Ug matrix
  ALLOCATE(RgMatrix(INhkl,INhkl,ITHREE),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RgMatrix")) CALL abort
  ! NB Rhkl are in INTEGER form [h,k,l] but are REAL to allow dot products etc.
  ALLOCATE(Rhkl(INhkl,ITHREE),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate Rhkl")) CALL abort
  ! Deviation parameter
  ALLOCATE(RDevPara(INhkl),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate Rhkl")) CALL abort
  ! allocate Ug arrays
  ALLOCATE(CUgMatNoAbs(INhkl,INhkl),STAT=IErr)! Ug Matrix without absorption
  IF(l_alert(IErr,"felixrefine","allocate CUgMatNoAbs")) CALL abort
  ALLOCATE(CUgMatPrime(INhkl,INhkl),STAT=IErr)! U'g Matrix of just absorption
  IF(l_alert(IErr,"felixrefine","allocate CUgMatPrime")) CALL abort
  ALLOCATE(CUgMat(INhkl,INhkl),STAT=IErr)! Ug+U'g Matrix, including absorption
  IF(l_alert(IErr,"felixrefine","allocate CUgMat")) CALL abort
  ! Matrix with numbers marking equivalent Ug's
  ALLOCATE(ISymmetryRelations(INhkl,INhkl),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate ISymmetryRelations")) CALL abort

  IThicknessCount= NINT((RFinalThickness-RInitialThickness)/RDeltaThickness) + 1

  !--------------------------------------------------------------------
  ! set up scattering factors, k-space resolution
  !--------------------------------------------------------------------
  CALL SetScatteringFactors(IScatterFactorMethodFLAG,IErr)
  IF(l_alert(IErr,"felixrefine","SetScatteringFactors")) CALL abort
  ! returns global RScattFactors depeding upon scattering method: Kirkland, Peng, etc.

  ! Calculate wavevector magnitude k and relativistic mass
  ! Electron Velocity in metres per second
  RElectronVelocity = &
        RSpeedOfLight*SQRT( ONE - ((RElectronMass*RSpeedOfLight**2) / &
        (RElectronCharge*RAcceleratingVoltage*THOUSAND+RElectronMass*RSpeedOfLight**2))**2 )
  ! Electron WaveLength in Angstroms
  RElectronWaveLength = RPlanckConstant / &
        (  SQRT(TWO*RElectronMass*RElectronCharge*RAcceleratingVoltage*THOUSAND) * &
        SQRT( ONE + (RElectronCharge*RAcceleratingVoltage*THOUSAND) / &
        (TWO*RElectronMass*RSpeedOfLight**2) )  ) * RAngstromConversion
  ! NB --- k=2pi/lambda and exp(i*k.r), physics convention, in reciprocal Angstroms
  RElectronWaveVectorMagnitude = TWOPI/RElectronWaveLength
  !resolution in k-space N.B. in cRED we define convergence angle as half the y-size
  RDeltaK = TWOPI*DEG2RADIAN*RFrameAngle/(RElectronWaveLength*REAL(ISizeX,RKIND))
  ! y-dimension of simulation, taking the input RConvergenceAngle as half-convergence angle
  ISizeY = NINT(TWOPI*TWO*RConvergenceAngle/RDeltaK)
  WRITE(SPrintString, FMT='(A11,I3,1x,A2,I3,A7)') "Simulation ",ISizeX,"x ",ISizeY," pixels"
  CALL message(LS,SPrintString)
  RRelativisticCorrection = ONE/SQRT( ONE - (RElectronVelocity/RSpeedOfLight)**2 )
  RRelativisticMass = RRelativisticCorrection*RElectronMass
  !conversion from Vg to Ug, h^2/(2pi*m0*e), see e.g. Kirkland eqn. C.5
  RScattFacToVolts = (RPlanckConstant**2)*(RAngstromConversion**2)/&
  (TWOPI*RElectronMass*RElectronCharge*RVolume)

  !--------------------------------------------------------------------
  ! ImageInitialisation
  !--------------------------------------------------------------------
  IPixelTotal = ISizeX*ISizeY
  ALLOCATE(IPixelLocation(IPixelTotal,2),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate IPixelLocation")) CALL abort
  ! we keep track of where a calculation goes in the image using the
  ! IPixelLocation array.  Remember fortran indexing is [row,col]=[y,x]
  lnd = 0
  DO ind = 1,ISizeY
    DO jnd = 1,ISizeX
      lnd = lnd + 1
      IPixelLocation(lnd,1) = ind
      IPixelLocation(lnd,2) = jnd
    END DO
  END DO

  !--------------------------------------------------------------------
  ! set up unit cell 
  ! total possible number of atoms/unit cell
  IMaxPossibleNAtomsUnitCell=SIZE(RBasisAtomPosition,1)*SIZE(RSymVec,1)
  ! over-allocate since actual size not known before calculation 
  ! (atoms in special positions will be duplicated)
  ! atoms, fractional unit cell
  ALLOCATE(RAtomPosition(IMaxPossibleNAtomsUnitCell,ITHREE),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RAtomPosition")) CALL abort
  ! atoms,in microscope reference frame, in Angstrom units
  ALLOCATE(RAtomCoordinate(IMaxPossibleNAtomsUnitCell,ITHREE),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RAtomCoordinate")) CALL abort
  ALLOCATE(SAtomLabel(IMaxPossibleNAtomsUnitCell),STAT=IErr) ! atom label
  IF(l_alert(IErr,"felixrefine","allocate SAtomLabel")) CALL abort
  ALLOCATE(SAtomName(IMaxPossibleNAtomsUnitCell),STAT=IErr) ! atom name
  IF(l_alert(IErr,"felixrefine","allocate SAtomName")) CALL abort
  ALLOCATE(RIsoDW(IMaxPossibleNAtomsUnitCell),STAT=IErr) ! Isotropic Debye-Waller factor
  IF(l_alert(IErr,"felixrefine","allocate RIsoDW")) CALL abort
  ALLOCATE(ROccupancy(IMaxPossibleNAtomsUnitCell),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate ROccupancy")) CALL abort
  ALLOCATE(IAtomicNumber(IMaxPossibleNAtomsUnitCell),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate IAtomicNumber")) CALL abort
  ! Anisotropic Debye-Waller factor, not yet functioning
  ALLOCATE(IAnisoDW(IMaxPossibleNAtomsUnitCell),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate IAnisoDW")) CALL abort
  !--------------------------------------------------------------------
  ! fill unit cell from basis and symmetry, then remove duplicates at special positions
  CALL UniqueAtomPositions(IErr)
  IF(l_alert(IErr,"felixrefine","UniqueAtomPositions")) CALL abort
  !?? RB could re-allocate RAtomCoordinate,SAtomName,RIsoDW,ROccupancy,
  !?? IAtomicNumber,IAnisoDW to match INAtomsUnitCell?

  ! From the unit cell we produce RaVecO, RbVecO, RcVecO in an orthogonal reference frame O
  ! with Xo // a and Zo perpendicular to the ab plane, in Angstrom units
  ! and reciprocal lattice vectors RarVecO, RbrVecO, RcrVecO in the same reference frame
  CALL ReciprocalLattice(IErr)
  IF(l_alert(IErr,"felixrefine","ReciprocalLattice")) CALL abort

  !--------------------------------------------------------------------
  ! Start of simulation-specific calculations, depending on crystal orientation
  !--------------------------------------------------------------------
  ! X, Y and Z are orthogonal vectors that defines the simulation
  ! Also referred to as the microscope reference frame M.
  ! The electron beam propagates along +Zm.
  ! The alpha rotation axis is along Ym.  Positive alpha rotation moves the field
  ! of view of the simulation along +Xm.
  ! In the crystal reference frame we read in reciprocal vectors RXDirC and RZDirC
  ! NB No check has been made to ensure that they are perpendicular
  ! RXDirO,RYDirO,RZDirO are UNIT reciprocal lattice vectors parallel to X,Y,Z
  RXDirO = RXDirC_0(1)*RarVecO + RXDirC_0(2)*RbrVecO + RXDirC_0(3)*RcrVecO
  RXDirO = RXDirO/SQRT(DOT_PRODUCT(RXDirO,RXDirO))
  RZDirO = RZDirC_0(1)*RaVecO + RZDirC_0(2)*RbVecO + RZDirC_0(3)*RcVecO
  RZDirO = RZDirO/SQRT(DOT_PRODUCT(RZDirO,RZDirO))
  RYDirO = CROSS(RZDirO,RXDirO)

  ! frame counter
  IFrame = 1
  DO WHILE(IFrame.LE.INFrames)
    WRITE(SPrintString, FMT='(A6,I3,A3)') "Frame ",IFrame,"..."
    CALL message(LS,dbg3,SPrintString)
    ! Increment frame angle, if it's not the first 
    IF(IFrame.GT.1) THEN
      RXDirOn = RXDirO-RZDirO*TAN(DEG2RADIAN*RFrameAngle)
      RZDirOn = RZDirO+RXDirO*TAN(DEG2RADIAN*RFrameAngle)
      RXDirO = RXDirOn/SQRT(DOT_PRODUCT(RXDirOn,RXDirOn))
      RZDirO = RZDirOn/SQRT(DOT_PRODUCT(RZDirOn,RZDirOn))
    END IF

    ! Create reciprocal lattice vectors in Microscope reference frame
    CALL CrystalOrientation(IErr)
    IF(l_alert(IErr,"felixrefine","CrystalOrientation")) CALL abort
    !--------------------------------------------------------------------
    ! Fill the list of reflections Rhkl
    Rhkl = ZERO
    RGlimit = 10.0*TWOPI    
    CALL HKLMake(RGlimit,IErr)
    IF(l_alert(IErr,"felixrefine","HKLMake")) CALL abort
    CALL message(LL,dbg7,"Rhkl matrix: ",NINT(Rhkl(1:INhkl,:)))

    !--------------------------------------------------------------------
    ! sort hkl in descending order of magnitude (not sure this is needed, really)
    CALL HKLSort(Rhkl,INhkl,IErr) 
    IF(l_alert(IErr,"felixrefine","SortHKL")) CALL abort
    ! Assign numbers to different reflections -> IhklsFrame, IhklsAll, INoOfHKLsFrame
    CALL HKLList(IErr)
    IF(l_alert(IErr,"felixrefine","SpecificReflectionDetermination")) CALL abort
    !--------------------------------------------------------------------
    ! Now we know which reflections are in this frame
    ! allocate image arrays for pixel-parallel simulations
    ! All the individual calculations go into RSimulatedPatterns later with MPI_GATHERV
    ! NB RSimulatedPatterns is a vector with respect to pixels, not a 2D image
    ALLOCATE(RSimulatedPatterns(INoOfHKLsFrame,IThicknessCount,IPixelTotal),STAT=IErr)
    IF(l_alert(IErr,"felixrefine","allocate RSimulatedPatterns")) CALL abort
    ! Images, NB Fortan arrays are [row,column]=[y,x]
    ALLOCATE(RImageSimi(ISizeY,ISizeX,INoOfHKLsFrame,IThicknessCount),STAT=IErr)
    IF(l_alert(IErr,"felixrefine","allocate RImageSimi")) CALL abort

    !--------------------------------------------------------------------
    ! calculate g vector magnitudes and components parallel to the surface
    CALL gVectors(IErr)
    IF(l_alert(IErr,"felixrefine","gVectors")) CALL abort
    !--------------------------------------------------------------------
    ! structure factor initialization
    ! Calculate Ug matrix for each entry in CUgMatNoAbs(1:INhkl,1:INhkl)
    CALL StructureFactorInitialisation(IErr)
    IF(l_alert(IErr,"felixrefine","StructureFactorInitialisation")) CALL abort
    ! NB IEquivalentUgKey and CUniqueUg allocated in here
    ! CUniqueUg vector produced here to later fill RIndependentVariable
      
    !--------------------------------------------------------------------
    ! calculate absorptive scattering factors
    !--------------------------------------------------------------------
    CALL SYSTEM_CLOCK( IStartTime2 )
    IF(IFrame.EQ.1) CALL message(LS,dbg3,"Starting absorption calculation... ")
    CALL Absorption (IErr)
    CALL message( LM, "Initial Ug matrix, with absorption (nm^-2)" )
    DO ind = 1,6
      WRITE(SPrintString,FMT='(3(I2,1X),A2,1X,6(F7.4,1X,F7.4,2X))') NINT(Rhkl(ind,:)),": ",100*CUgMat(ind,1:6)
      CALL message( LM,dbg3, SPrintString)
    END DO
    IF(l_alert(IErr,"felixrefine","Absorption")) CALL abort
    IF(IFrame.EQ.1) THEN 
      CALL PrintEndTime(LS,IStartTime2, "Absorption" )
      CALL SYSTEM_CLOCK( IStartTime2 )
    END IF

    !--------------------------------------------------------------------
    ! set up arrays for pixel-parallel simulations
    !--------------------------------------------------------------------
    RSimulatedPatterns = ZERO
    ! The pixels to be calculated by this core  
    ILocalPixelCountMin= (IPixelTotal*(my_rank)/p)+1
    ILocalPixelCountMax= (IPixelTotal*(my_rank+1)/p) 
    ALLOCATE(RIndividualReflections(INoOfHKLsFrame,IThicknessCount,&
           (ILocalPixelCountMax-ILocalPixelCountMin)+1),STAT=IErr)
    IF(l_alert(IErr,"felixrefine","allocate RIndividualReflections")) CALL abort

    ! position of pixels calculated by this core, ILocalPixelOffset & ILocalNPix are global variables
    ALLOCATE(ILocalPixelOffset(p),ILocalNPix(p),STAT=IErr)
    IF(l_alert(IErr,"felixrefine","allocate ILocalPixelOffset")) CALL abort
    DO ind = 1,p
      ILocalPixelOffset(ind) = (IPixelTotal*(ind-1)/p)*INoOfHKLsFrame*IThicknessCount
      ILocalNPix(ind) = (((IPixelTotal*(ind)/p) - (IPixelTotal*(ind-1)/p)))* &
            INoOfHKLsFrame*IThicknessCount    
    END DO

    !--------------------------------------------------------------------
    ! simulation (different local pixels for each core)
    !--------------------------------------------------------------------
    IF(IFrame.EQ.1) CALL message(LS,"Bloch wave calculation...")    
    ! Reset simulation   
    RIndividualReflections = ZERO
    DO knd = ILocalPixelCountMin,ILocalPixelCountMax,1
      ! fills array for each pixel number knd (y & x coordinates are in IPixelLocation)
      CALL BlochCoefficientCalculation(IPixelLocation(knd,1),IPixelLocation(knd,2),knd, &
              ILocalPixelCountMin, nBeams, RThickness,RKn, IErr)
      IF(l_alert(IErr,"Simulate","BlochCoefficientCalculation")) CALL abort
    END DO
    !===================================== ! MPI gatherv into RSimulatedPatterns
    CALL MPI_GATHERV(RIndividualReflections,SIZE(RIndividualReflections),MPI_DOUBLE_PRECISION,&
         RSimulatedPatterns,ILocalNPix,ILocalPixelOffset,MPI_DOUBLE_PRECISION,&
         root,MPI_COMM_WORLD,IErr)
    !This brodcast is not strictly necessary but keeps all cores synchronised
    CALL MPI_BCAST(RIndividualReflections,SIZE(RIndividualReflections),MPI_DOUBLE_PRECISION,&
         root,MPI_COMM_WORLD,IErr)
    !=====================================
    IF(l_alert(IErr,"SimulateAndFit","MPI_GATHERV")) CALL abort
    ! put 1D array RSimulatedPatterns into 2D image RImageSimi (should be done with RESHAPE?)
    ! NB dimensions of RSimulatedPatterns(INoOfHKLsFrame,IThicknessCount,IPixelTotal)
    ! and RImageSimi(height, width, INoOfHKLsFrame,IThicknessCount )
    RImageSimi = ZERO
    lnd = 0
    DO ind = 1,ISizeY
      DO jnd = 1,ISizeX
        lnd = lnd+1
        RImageSimi(ind,jnd,:,:) = RSimulatedPatterns(:,:,lnd)
      END DO
    END DO

    !--------------------------------------------------------------------
    ! make compound 000 image
    IF(IFrame.EQ.1)THEN!first frame, set up the image
      ALLOCATE(RBrightField(ISizeY,ISizeX), STAT=IErr)
      IF(l_alert(IErr,"felixrefine","allocate RBrightField")) CALL abort
      RBrightField = RImageSimi(:,:,1,1)
    ELSE! subsequent frame, append the image
      ALLOCATE(RTempImage (ISizeY,ISizeX+SIZE(RBrightField,DIM=2)), STAT=IErr)
      IF(l_alert(IErr,"felixrefine","allocate RBrightField")) CALL abort
      RTempImage(:,1:SIZE(RBrightField,DIM=2)) = RBrightField
      RTempImage(:,SIZE(RBrightField,DIM=2)+1:) = RImageSimi(:,:,1,1)
      DEALLOCATE(RBrightField, STAT=IErr)
      ALLOCATE(RBrightField (ISizeY,SIZE(RTempImage,DIM=2)), STAT=IErr)
      RBrightField = RTempImage
      DEALLOCATE(RTempImage, STAT=IErr)
    END IF




    ! Gaussian blur to match experiment using global variable RBlurRadius
    IF (RBlurRadius.GT.TINY) THEN
      DO ind=1,INoOfHKLsFrame
        DO jnd=1,IThicknessCount
          CALL BlurG(RImageSimi(:,:,ind,jnd),ISizeX,ISizeY,RBlurRadius,IErr)
        END DO
      END DO
    END IF

    IF(l_alert(IErr,"felixrefine","Simulate")) CALL abort
    IF(IFrame.EQ.1) CALL PrintEndTime(LS,IStartTime2, "Simulation" )

    ! output, multiple thicknesses
    IF(my_rank.EQ.0) THEN
      WRITE(SPrintString,FMT='(A24,I3,A12)') &
        "Writing simulations for ", IThicknessCount," thicknesses"
      CALL message(LS,SPrintString)
      DO ind = 1,IThicknessCount
        CALL WriteIterationOutput(ind,IErr)
        IF(l_alert(IErr,"felixrefine","WriteIterationOutput")) CALL abort 
      END DO  
    END IF 

    !--------------------------------------------------------------------
    ! deallocate memory used in each frame
    !--------------------------------------------------------------------
    DEALLOCATE(RSimulatedPatterns,STAT=IErr)
    DEALLOCATE(RImageSimi,STAT=IErr)
    DEALLOCATE(IEquivalentUgKey,STAT=IErr)
    DEALLOCATE(CUniqueUg,STAT=IErr)
    DEALLOCATE(RIndividualReflections,STAT=IErr)
    DEALLOCATE(ILocalPixelOffset,STAT=IErr)
    DEALLOCATE(ILocalNPix,STAT=IErr)
    !--------------------------------------------------------------------
    ! frame loop
    IFrame = IFrame + 1
  END DO

  !--------------------------------------------------------------------
  ! finish off: deallocations for variables used in all frames
  !--------------------------------------------------------------------
  DEALLOCATE(RgPoolMag,STAT=IErr)
  DEALLOCATE(RgPool,STAT=IErr)
  DEALLOCATE(RgMatrix,STAT=IErr)
  DEALLOCATE(RgDotNorm,STAT=IErr)
  DEALLOCATE(Rhkl,STAT=IErr)
  DEALLOCATE(CUgMat,STAT=IErr)
  DEALLOCATE(CUgMatNoAbs,STAT=IErr)
  DEALLOCATE(CUgMatPrime,STAT=IErr)
  DEALLOCATE(ISymmetryRelations,STAT=IErr)
  DEALLOCATE(RAtomPosition,STAT=IErr)
  DEALLOCATE(SAtomName,STAT=IErr)
  DEALLOCATE(RIsoDW,STAT=IErr)
  DEALLOCATE(ROccupancy,STAT=IErr)
  DEALLOCATE(IAtomicNumber,STAT=IErr)
  DEALLOCATE(IAnisoDW,STAT=IErr)
  DEALLOCATE(RAtomCoordinate,STAT=IErr)
  DEALLOCATE(IPixelLocation,STAT=IErr)

  CALL message( LS, "--------------------------------" )
  CALL PrintEndTime( LS, IStartTime, "Calculation" )
  CALL message( LS, "--------------------------------")
  CALL message( LS, "||||||||||||||||||||||||||||||||")
  
  ! clean shutdown
  CALL MPI_Finalize(IErr)
  STOP
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  CONTAINS

  SUBROUTINE abort
    IErr=1
    IF(l_alert(IErr,"felixrefine","ABORTING")) CONTINUE
    CALL MPI_Abort(MPI_COMM_WORLD,1,IErr)
    STOP
  END SUBROUTINE abort


  SUBROUTINE BlurG(RImageToBlur,IPixX,IPixY,RBlurringRadius,IErr)

    USE MyNumbers
    USE MPI
    USE message_mod

    IMPLICIT NONE

    REAL(RKIND),DIMENSION(IPixX,IPixY),INTENT(INOUT) :: RImageToBlur
    INTEGER(IKIND),INTENT(IN) :: IPixX,IPixY
    REAL(RKIND),INTENT(IN) :: RBlurringRadius
    INTEGER(IKIND),INTENT(OUT) :: IErr

    REAL(RKIND),DIMENSION(IPixX,IPixY) :: RTempImage,RShiftImage
    INTEGER(IKIND) :: ind,jnd,IKernelRadius,IKernelSize
    REAL(RKIND),DIMENSION(:), ALLOCATABLE :: RGauss1D
    REAL(RKIND) :: Rind,Rsum,Rmin,Rmax

    ! get min and max of input image
    Rmin=MINVAL(RImageToBlur)
    Rmax=MAXVAL(RImageToBlur)

    ! set up a 1D kernel of appropriate size  
    IKernelRadius=NINT(3*RBlurringRadius)
    ALLOCATE(RGauss1D(2*IKernelRadius+1),STAT=IErr)!ffs
    Rsum=0
    DO ind=-IKernelRadius,IKernelRadius
      Rind=REAL(ind)
      RGauss1D(ind+IKernelRadius+1)=EXP(-(Rind**2)/(2*(RBlurringRadius**2)))
      Rsum=Rsum+RGauss1D(ind+IKernelRadius+1)
      IF(ind==0) IErr=78 
    END DO
    RGauss1D=RGauss1D/Rsum!normalise
    RTempImage=RImageToBlur*0_RKIND !reset the temp image 

    ! apply the kernel in direction 1
    DO ind = -IKernelRadius,IKernelRadius
       IF (ind.LT.0) THEN
          RShiftImage(1:IPixX+ind,:)=RImageToBlur(1-ind:IPixX,:)
          DO jnd = 1,1-ind!edge fill on right
             RShiftImage(IPixX-jnd+1,:)=RImageToBlur(IPixX,:)
          END DO
       ELSE
          RShiftImage(1+ind:IPixX,:)=RImageToBlur(1:IPixX-ind,:)
          DO jnd = 1,1+ind!edge fill on left
             RShiftImage(jnd,:)=RImageToBlur(1,:)
          END DO
       END IF
       RTempImage=RTempImage+RShiftImage*RGauss1D(ind+IKernelRadius+1)
    END DO

    ! make the 1D blurred image the input for the next direction
    RImageToBlur=RTempImage
    RTempImage=RImageToBlur*0_RKIND ! reset the temp image

    ! apply the kernel in direction 2  
    DO ind = -IKernelRadius,IKernelRadius
       IF (ind.LT.0) THEN
          RShiftImage(:,1:IPixY+ind)=RImageToBlur(:,1-ind:IPixY)
          DO jnd = 1,1-ind!edge fill on bottom
             RShiftImage(:,IPixY-jnd+1)=RImageToBlur(:,IPixY)
          END DO
       ELSE
          RShiftImage(:,1+ind:IPixY)=RImageToBlur(:,1:IPixY-ind)
          DO jnd = 1,1+ind!edge fill on top
             RShiftImage(:,jnd)=RImageToBlur(:,1)
          END DO
       END IF
       RTempImage=RTempImage+RShiftImage*RGauss1D(ind+IKernelRadius+1)
    END DO
    DEALLOCATE(RGauss1D,STAT=IErr)

    ! set intensity range of output image to match that of the input image
    RTempImage=RTempImage-MINVAL(RTempImage)
    RTempImage=RTempImage*(Rmax-Rmin)/MAXVAL(RTempImage)+Rmin
    ! return the blurred image
    RImageToBlur=RTempImage;

  END SUBROUTINE BlurG    

END PROGRAM Felixrefine

