!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: 1.3
! Date: 13-05-2024
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build: g-vector limit 
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
  USE image_initialisation_mod
  USE setup_space_group_mod
  USE crystallography_mod
  USE ug_matrix_mod
  USE refinementcontrol_mod       ! CONTAINS Simulate and SimulateAndFit
  USE write_output_mod
  USE simplex_mod

  USE IConst; USE RConst; USE SConst
  USE IPara;  USE RPara;  USE CPara; USE SPara;
  USE BlochPara 
  USE IChannels

  ! local variable definitions
  IMPLICIT NONE
 
  INTEGER(IKIND) :: IErr,ind,jnd,knd,lnd,mnd,nnd,Iter,ICutOff,IHOLZgPoolMag,&
        IBSMaxLocGVecAmp,ILaueLevel,INumTotalReflections,ITotalLaueZoneLevel,&
        IPrintFLAG,ICycle,INumInitReflections,IZerothLaueZoneLevel,xnd,&
        INumFinalReflections,IThicknessIndex,IVariableType,IArrayIndex,&
        IAnisotropicDebyeWallerFactorElementNo,IStartTime,IStartTime2
  INTEGER(4) :: IErr4
  REAL(RKIND) :: REmphasis,RLaueZoneGz,RMaxGMag,&
        RPvecMag,RScale,RMaxUgStep,Rdx,RStandardDeviation,RMean,RGzUnitVec,&
        RMinLaueZoneValue,Rdf,RLastFit,RBestFit,RMaxLaueZoneValue,&
        RMaxAcceptanceGVecMag,RandomSign,RLaueZoneElectronWaveVectorMag,&
        RvarMin,RfitMin,RFit0,Rconvex,Rtest,Rplus,Rminus,RdeltaX,RdeltaY
  REAL(RKIND),DIMENSION(100) :: RTemp!temporary holder for refinement variables
  INTEGER(IKIND),DIMENSION(100) :: ITemp,Itemp2!temporary holder for refinement type and atom
  REAL(RKIND),DIMENSION(ITHREE) :: R3var,R3fit
  INTEGER(IKIND),DIMENSION(10) :: INoOfVariablesForRefinementType

  ! allocatable arrays
  INTEGER(IKIND),DIMENSION(:),ALLOCATABLE :: IOriginGVecIdentifier
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RSimplexFoM,RIndependentVariable,&
        RCurrentVar,RVar0,RLastVar,RPvec,RFitVec,RLastVec
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: RSimplexVariable,RgDummyVecMat,&
        RgPoolMagLaue,RTestImage,ROnes,RVarMatrix,RSimp
  
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
  IInitialSimulationFLAG = 1

  ! MPI initialization
  CALL MPI_Init(IErr4) 
  IF(l_alert(INT(REAL(IErr4)),"felixrefine","MPI_Init")) CALL abort
  CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,IErr4) ! get rank of the current process
  IF(l_alert(INT(REAL(IErr4)),"felixrefine","MPI_Comm_rank")) CALL abort
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,IErr4) ! get size of the current communicator
  IF(l_alert(INT(REAL(IErr4)),"felixrefine","MPI_Comm_size")) CALL abort

  ! startup terminal output
  CALL message(LS,"-----------------------------------------------------------------")
  INCLUDE "version.txt"
!!$#ifdef git
!!$  CALL message(LS,"felixrefine: ", TRIM("GITVERSION"))
!!$  CALL message(LS,"             ", TRIM("GITBRANCH"))
!!$  CALL message(LS,"             ", TRIM("COMPILED"))
!!$#else
!!$  CALL message(LS,"felixrefine: ", "see https://github.com/WarwickMicroscopy/Felix for version")
!!$#endif
  CALL message(LS,"-----------------------------------------------------------------")
  CALL message(LS,"total number of MPI ranks ", p, ", screen messages via rank", my_rank)
  CALL message(LS,"-----------------------------------------------------------------")

  ! timing setup
  CALL SYSTEM_CLOCK( IStartTime,IClockRate )

  !--------------------------------------------------------------------
  ! input section 
  !--------------------------------------------------------------------
  
  CALL ReadInpFile(IErr) ! felix.inp
  IF(l_alert(IErr,"felixrefine","ReadInpFile")) CALL abort
  CALL SetMessageMode( IWriteFLAG, IErr )
  IF(l_alert(IErr,"felixrefine","set_message_mod_mode")) CALL abort

  CALL read_cif(IErr) ! felix.cif ! some allocations are here
  IF(l_alert(IErr,"felixrefine","ReadCif")) CALL abort

  CALL ReadHklFile(IErr) ! the list of hkl's to input/output
  IF(l_alert(IErr,"felixrefine","ReadHklFile")) CALL abort

  ! read experimental images (if in refinement mode)
  IF (ISimFLAG.EQ.0) THEN ! it's a refinement, so read experimental images
    ALLOCATE(RImageExpi(2*IPixelCount,2*IPixelCount,INoOfLacbedPatterns),STAT=IErr)  
    IF(l_alert(IErr,"felixrefine","allocate RImageExpi")) CALL abort
    CALL ReadExperimentalImages(IErr)
    IF(l_alert(IErr,"felixrefine","ReadExperimentalImages")) CALL abort
  ELSE ! not refinement, simulation only
    CALL message(LS,"Simulation only")
  END IF
  !--------------------------------------------------------------------
  ! set up scattering factors, relativistic electrons, reciprocal lattice
  !--------------------------------------------------------------------
  CALL SetScatteringFactors(IScatterFactorMethodFLAG,IErr)
  IF(l_alert(IErr,"felixrefine","SetScatteringFactors")) CALL abort
  ! returns global RScattFactors depeding upon scattering method: Kirkland, Peng, etc.

  ! Calculate wavevector magnitude k and relativistic mass
  ! Electron Velocity in metres per second
  RElectronVelocity = &
        RSpeedOfLight*SQRT( ONE - ((RElectronMass*RSpeedOfLight**2) / &
        (RElectronCharge*RAcceleratingVoltage*THOUSAND+RElectronMass*RSpeedOfLight**2))**2 )
  ! Electron WaveLength in metres
  RElectronWaveLength = RPlanckConstant / &
        (  SQRT(TWO*RElectronMass*RElectronCharge*RAcceleratingVoltage*THOUSAND) * &
        SQRT( ONE + (RElectronCharge*RAcceleratingVoltage*THOUSAND) / &
        (TWO*RElectronMass*RSpeedOfLight**2) )  ) * RAngstromConversion
  ! NB --- k=2pi/lambda and exp(i*k.r), physics convention
  RElectronWaveVectorMagnitude=TWOPI/RElectronWaveLength
  RRelativisticCorrection = ONE/SQRT( ONE - (RElectronVelocity/RSpeedOfLight)**2 )
  RRelativisticMass = RRelativisticCorrection*RElectronMass
  !conversion from Vg to Ug, h^2/(2pi*m0*e), see e.g. Kirkland eqn. C.5
  RScattFacToVolts=(RPlanckConstant**2)*(RAngstromConversion**2)/&
  (TWOPI*RElectronMass*RElectronCharge*RVolume)
  ! Creates reciprocal lattice vectors in Microscope reference frame
  CALL ReciprocalLattice(IErr)
  IF(l_alert(IErr,"felixrefine","ReciprocalLattice")) CALL abort
  !resolution in k-space
  RDeltaK = TWOPI*RConvergenceAngle/REAL(IPixelCount,RKIND)
!IF(my_rank.EQ.0)PRINT*, "Delta K", RDeltaK  

  !--------------------------------------------------------------------
  ! allocate atom and Debye-Waller factor arrays
  !--------------------------------------------------------------------
  ! Reset the basis so that atomic coordinate refinement is possible 
  IF(IRefineMode(2).EQ.1) CALL PreferredBasis(IErr)!A crystallography subroutine
  ! total possible atoms/unit cell
  IMaxPossibleNAtomsUnitCell=SIZE(RBasisAtomPosition,1)*SIZE(RSymVec,1)
  ! over-allocate since actual size not known before calculation of unique positions 
  ! (atoms in special positions will be duplicated)

  ! allocations using that RBasisAtomPosition, RSymVec have now been setup
  ! fractional unit cell coordinates are used for RAtomPosition, like BasisAtomPosition
  ALLOCATE(RAtomPosition(IMaxPossibleNAtomsUnitCell,ITHREE),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RAtomPosition")) CALL abort
  ! RAtomCoordinate is in the microscope reference frame in Angstrom units
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
  ! Anisotropic Debye-Waller factor
  ALLOCATE(IAnisoDW(IMaxPossibleNAtomsUnitCell),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate IAnisoDW")) CALL abort

  !--------------------------------------------------------------------
  ! set up unique atom positions, reflection pool
  ! fills unit cell from basis and symmetry, removes duplicate atoms at special positions
  CALL UniqueAtomPositions(IErr)
  IF(l_alert(IErr,"felixrefine","UniqueAtomPositions")) CALL abort
  !?? RB could re-allocate RAtomCoordinate,SAtomName,RIsoDW,ROccupancy,
  !?? IAtomicNumber,IAnisoDW to match INAtomsUnitCell?

  ! Fill the list of reflections Rhkl (global variable)
  ! NB Rhkl are in INTEGER form [h,k,l] but are REAL to allow dot products etc.
  ! ** now done in HKLMake ** ALLOCATE(Rhkl(INhkl,ITHREE),STAT=IErr)
  !IF(l_alert(IErr,"felixrefine","allocate Rhkl")) CALL abort
  ! make the beam pool, uses RgLimit    
  CALL HKLMake(IErr)
  IF(l_alert(IErr,"felixrefine","HKLMake")) CALL abort
  CALL message(LL,dbg7,"Rhkl matrix: ",NINT(Rhkl(1:INhkl,:)))

  !--------------------------------------------------------------------
  ! sort hkl in descending order of magnitude
  CALL HKLSort(Rhkl,INhkl,IErr) 
  IF(l_alert(IErr,"felixrefine","SortHKL")) CALL abort
  !?? RB may result in an error when the reflection pool does not reach
  !?? the highest hkl of the experimental data? 

  ! Assign numbers to different reflections -> IOutputReflections, INoOfLacbedPatterns
  CALL HKLList(IErr)
  IF(l_alert(IErr,"felixrefine","SpecificReflectionDetermination")) CALL abort

  !--------------------------------------------------------------------
  ! allocations
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
  ! Matrix of their magnitudes 
!  ALLOCATE(RgMatrixMagnitude(INhkl,INhkl),STAT=IErr)
!  IF(l_alert(IErr,"felixrefine","allocate RgMatrixMagnitude")) CALL abort
 
  !--------------------------------------------------------------------
  ! calculate g vector list, magnitudes and components parallel to the surface
  CALL gVectors(IErr)
  IF(l_alert(IErr,"felixrefine","gVectors")) CALL abort

  !--------------------------------------------------------------------
  ! sort into Laue Zones
  !Dummy matrix, used in HOLZ calculation (not working)
  ALLOCATE(RgDummyVecMat(INhkl,ITHREE),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RgDummyVecMat")) CALL abort
  ICutOff = 1
  DO ind=1,INhkl
    ! If a g-vector has a non-zero z-component it is not in the ZOLZ
    IF(ABS(RgPool(ind,3)).GT.TINY.AND.ICutOff.NE.0) THEN
      RGzUnitVec=ABS(RgPool(ind,3))
      ICutOff=0
    END IF
  ENDDO
  
  ! This should occur when no higher Laue zones (IHolzFLAG off)
  IF(ICutOff.EQ.1) THEN ! all g-vectors have z-component equal to zero
    RGzUnitVec=ZERO
  END IF

  IF( ICutOff.EQ.1 .AND. IHolzFLAG.EQ.1 ) IErr = 1
  IF( l_alert(IErr, "felixrefine", &
        "fill Laue Zones. IHolzFLAG = 1 in felix.inp, however no higher order g-vectors " &
        //"were found. Continuing with zeroth-order Laue zone only.") ) IErr = 0
  
  RgDummyVecMat=RgPool
  WHERE(ABS(RgPool(:,3)).GT.TINY) ! higher order Laue zones cases
    RgDummyVecMat(:,3)=RgDummyVecMat(:,3)/RGzUnitVec 
  END WHERE ! divide zero is not a concern as condition matches above

  ! min & max Laue Zones 
  RMaxLaueZoneValue=MAXVAL(RgDummyVecMat(:,3),DIM=1)
  RMinLaueZoneValue=MINVAL(RgDummyVecMat(:,3),DIM=1)
  ITotalLaueZoneLevel=NINT(RMaxLaueZoneValue+ABS(RMinLaueZoneValue)+1,IKIND)
  
  DEALLOCATE(RgDummyVecMat, STAT=IErr)
  IF(l_alert(IErr,"felixrefine","deallocate RgDummyVecMat")) CALL abort

  !--------------------------------------------------------------------
  ! (optional) calculate higher order Laue zone 
  ALLOCATE(RgPoolMagLaue(INhkl,ITotalLaueZoneLevel),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RgPoolMagLaue")) CALL abort

  ! calculate higher order Laue zones
  IF(RAcceptanceAngle.NE.ZERO.AND.IHOLZFLAG.EQ.1) THEN
    !?? currently not working, unused variables here, lots needs to be checked
    INumtotalReflections=0
    DO ind=1,ITotalLaueZoneLevel
      ILaueLevel=ind-IZerothLaueZoneLevel
      RLaueZoneGz=RGzUnitVec*ILaueLevel
      DO jnd=1,INhkl
        IF(RgPool(jnd,3).GE.(RLaueZoneGz-TINY).AND. &
             RgPool(jnd,3).LE.(RLaueZoneGz+TINY)) THEN
          RgPoolMagLaue(jnd,ind)=SQRT((RgPool(jnd,1)**2)+(RgPool(jnd,2)**2))              
        ELSE
          RgPoolMagLaue(jnd,ind)=NEGHUGE
        END IF
      END DO
      INumInitReflections = COUNT(RgPoolMagLaue(:,ind).NE.NEGHUGE)
      RLaueZoneElectronWaveVectorMag = RElectronWaveVectorMagnitude-ABS(RLaueZoneGz)           
      RMaxAcceptanceGVecMag = (RLaueZoneElectronWaveVectorMag * &
            TAN(RAcceptanceAngle*DEG2RADIAN))
      WHERE(ABS(RgPoolMagLaue(:,ind)).GT.RMaxAcceptanceGVecMag)
        RgPoolMagLaue(:,ind)=NEGHUGE
      END WHERE
      INumFinalReflections=COUNT(RgPoolMagLaue(:,ind).NE.NEGHUGE)
      INumTotalReflections=INumTotalReflections+INumInitReflections
    END DO

    knd=0
    DO ind=1,INhkl
      IF(SUM(RgPoolMagLaue(ind,:))/REAL(ITotalLaueZoneLevel,RKIND).GT.NEGHUGE) THEN
        knd=knd+1
      END IF
    END DO
    IHOLZgPoolMag=knd
    ALLOCATE(IOriginGVecIdentifier(IHOLZgPoolMag),STAT=IErr)
    IF(l_alert(IErr,"felixrefine","allocate IOriginGVecIdentifier")) CALL abort
    IOriginGVecIdentifier=0
    knd=1
    DO ind=1,INhkl
      IF((SUM(RgPoolMagLaue(ind,:))/REAL(ITotalLaueZoneLevel,RKIND)).GT.NEGHUGE) THEN
        IOriginGVecIdentifier(knd)=ind
        knd=knd+1
      END IF
    END DO

  END IF

  !--------------------------------------------------------------------
  ! reduce number of reflections for finite acceptance angle 
  !--------------------------------------------------------------------
  ! acceptance angle (NB input of ZERO actually means unlimited)
  IF(RAcceptanceAngle.NE.ZERO) THEN!if acceptance angle is small, g-vectors will be limited
    IF (IHOLZFLAG.EQ.0) THEN!ZOLZ only
      RMaxAcceptanceGVecMag=(RElectronWaveVectorMagnitude*TAN(RAcceptanceAngle*DEG2RADIAN))
      IF(RgPoolMag(INhkl).GT.RMaxAcceptanceGVecMag) RMaxGMag = RMaxAcceptanceGVecMag 
    ELSE!HOLZ too(not working?)
      IBSMaxLocGVecAmp=MAXVAL(IOriginGVecIdentifier)
      RMaxGMag=RgPoolMag(IBSMaxLocGVecAmp)
      IF(RgPoolMag(IBSMaxLocGVecAmp).GE.RgPoolMag(INhkl)) &
        RMaxGMag = RgPoolMag(INhkl)
    END IF
    ! count reflections up to cutoff magnitude
    jnd=INhkl
    INhkl=0
    DO ind=1,jnd
      IF (ABS(RgPoolMag(ind)).LE.RMaxGMag) INhkl=INhkl+1
    END DO
    !check the acceptance angle is big enough to produce the LACBED patterns
    IF (INhkl.LT.INoOfLacbedPatterns) THEN
      IErr=1
      IF(l_alert(IErr,"felixrefine","Acceptance angle is too small! Please increase it or set to 0.0")) CALL abort
    END IF
    
  END IF
  
  
  ! deallocation
  DEALLOCATE(RgPoolMagLaue)!
  IF (RAcceptanceAngle.NE.ZERO.AND.IHOLZFLAG.EQ.1) THEN
    DEALLOCATE(IOriginGVecIdentifier)
  END IF

  !--------------------------------------------------------------------
  ! allocate Ug arrays
  !--------------------------------------------------------------------
  ! Ug matrix etc.
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
  CALL message(LS,dbg3,"Starting absorption calculation... ")
  CALL Absorption (IErr)
!  CALL message( LM, "Initial Ug matrix, with absorption (nm^-2)" )
!  DO ind = 1,6
!      WRITE(SPrintString,FMT='(3(I2,1X),A2,1X,6(F7.4,1X,F7.4,2X))') NINT(Rhkl(ind,:)),": ",100*CUgMat(ind,1:6)
!    CALL message( LM,dbg3, SPrintString)
!  END DO
  IF(l_alert(IErr,"felixrefine","Absorption")) CALL abort
  CALL PrintEndTime(LS,IStartTime2, "Absorption" )
  CALL SYSTEM_CLOCK( IStartTime2 )

  !--------------------------------------------------------------------
  ! INoOfVariables calculated depending upon Ug and non-Ug refinement
  !--------------------------------------------------------------------
  ! Ug refinement is a special case, cannot do any other refinement alongside
  ! To count the number of Independent Variables:
  ! Put the variable into Rtemp, the variable type into Itemp, the atom into
  ! Ttemp2.  Increment jnd.  Max number of refinement variables is 100!
  ! Once counting is done, jnd is used to allocate the size of 'proper' refinement variables
  ! RIndependentVariable and IIndependentVariableType.  Itemp,Rtemp are then copied over
  !Itemp2=0 means the atom is not being refined (or the refinement isn't for an
  !atom), initialise at zero
  ITemp2=0
  IF(ISimFLAG.EQ.0) THEN
    jnd=1
    IF(IRefineMode(1).EQ.1) THEN ! It's a Ug refinement, A
      IUgOffset=1  ! can skip Ug's in refinement using offset, 1 is inner potential
      DO ind = 1+IUgOffset,INoofUgs+IUgOffset
        IF ( ABS(REAL(CUniqueUg(ind),RKIND)).GE.RTolerance ) THEN
          Rtemp(jnd) = REAL(CUniqueUg(ind),RKIND)
          Itemp(jnd)=1
          jnd=jnd+1
        END IF
        IF ( ABS(AIMAG(CUniqueUg(ind))).GE.RTolerance ) THEN 
          Rtemp(jnd) = AIMAG(CUniqueUg(ind))
          Itemp(jnd)=1
          jnd=jnd+1
        END IF
      END DO
      IF (IAbsorbFLAG.EQ.1) THEN ! proportional absorption
        INoOfVariables = jnd ! the last variable is for absorption
        Rtemp(jnd) = RAbsorptionPercentage
        Itemp(jnd)=1
        jnd=jnd+1
      END IF

    ELSE ! It's not a Ug refinement, so count refinement variables
      ! Variables can be refined together
      
      IF (IRefineMode(2).EQ.1) THEN! Atom coordinate refinement, B
        CALL SetupAtomMovements(IErr)!returns IAtomMoveList and RVector
        IF(l_alert(IErr,"felixrefine","SetupAtomMovements")) CALL abort
        DO ind=1,SIZE(IAtomMoveList)
            Rtemp(jnd)=DOT_PRODUCT(RBasisAtomPosition(IAtomMoveList(ind),:),RVector(ind,:))
            Itemp(jnd)=2
            Itemp2(jnd)=IAtomMoveList(ind)!don't actually need this, just for completeness at present
            jnd=jnd+1
        END DO
      END IF

      IF (IRefineMode(3).EQ.1) THEN ! Occupancy, C
        DO ind=1,SIZE(IAtomsToRefine)
            Rtemp(jnd)=RBasisOccupancy(IAtomsToRefine(ind))
            Itemp(jnd)=3
            Itemp2(jnd)=IAtomsToRefine(ind)
            jnd=jnd+1
        END DO
      END IF

      IF (IRefineMode(4).EQ.1) THEN ! Isotropic DW, D
        DO ind=1,SIZE(IAtomsToRefine)
            Rtemp(jnd)=RBasisIsoDW(IAtomsToRefine(ind))
            Itemp(jnd)=4
            Itemp2(jnd)=IAtomsToRefine(ind)
            jnd=jnd+1
        END DO
      END IF
      
      IF (IRefineMode(5).EQ.1) THEN  ! Anisotropic DW, E
        IErr=1!not yet implemented!!!
        IF(l_alert(IErr,"felixrefine",&
            "Anisotropic Debye-Waller factor refinement not yet implemented, sorry")) CALL abort
      END IF
      
      IF (IRefineMode(6).EQ.1) THEN ! Lattice parameters, F
        IF(l_alert(IErr,"felixrefine","ConvertSpaceGroupToNumber")) CALL abort
        !This section needs work to include rhombohedral cells and non-standard
        !settings!!!
        Rtemp(jnd)=RLengthX!This is a free parameter for all lattice types
        Itemp(jnd)=6
        jnd=jnd+1
        IF (ISpaceGrp.LT.75) THEN!triclinic,monoclinic,orthorhombic
          Rtemp(jnd)=RLengthY
          Itemp(jnd)=16
          jnd=jnd+1
          Rtemp(jnd)=RLengthZ
          Itemp(jnd)=26
          jnd=jnd+1
        ELSE IF (ISpaceGrp.GT.142.AND.ISpaceGrp.LT.168) THEN!rhombohedral
          IErr=1!need to work out R- vs H- settings!!!
          PRINT*,"Rhombohedral R- and H- cells not yet implemented for unit cell refinement"
        ELSE IF ((ISpaceGrp.GT.167.AND.ISpaceGrp.LT.195).OR.&!Hexagonal
                   (ISpaceGrp.GT.74.AND.ISpaceGrp.LT.143)) THEN!Tetragonal
          Rtemp(jnd)=RLengthZ
          Itemp(jnd)=26
          jnd=jnd+1
        END IF
      END IF

      IF (IRefineMode(7).EQ.1) THEN ! Unit cell angles, G
        IErr=1!not yet implemented!!!
        IF(l_alert(IErr,"felixrefine",&
            "Unit cell angle refinement not yet implemented, sorry")) CALL abort
        jnd=jnd+3
      END IF
      
      IF (IRefineMode(8).EQ.1) THEN ! Convergence angle, H
        Rtemp(jnd)=RConvergenceAngle
        Itemp(jnd)=8
        jnd=jnd+1
      END IF

      IF (IRefineMode(9).EQ.1) THEN ! kV, I
        jnd=jnd+1
        IErr=1!not yet implemented!!!
        IF(l_alert(IErr,"felixrefine",&
            "kV refinement not yet implemented, sorry")) CALL abort
      END IF

    END IF
    ! Total number of independent variables
    INoOfVariables = jnd-1
    IF (INoOfVariables.EQ.0) THEN 
      ! there's no refinement requested, say so and quit
      IErr = 1
      IF(l_alert(IErr,"felixrefine",&
            "No refinement variables! Check IRefineModeFLAG in felix.inp. "// &
            "Valid refine modes are A,B,C,D,E,F,G,H,I,S")) CALL abort
    END IF
    IF (INoOfVariables.EQ.1 ) THEN 
      CALL message(LS,"Only one independent variable")
    ELSE
      CALL message(LS,"Number of independent variables = ",INoOfVariables)
    END IF
    ALLOCATE(RIndependentVariable(INoOfVariables),STAT=IErr) 
    ALLOCATE(RIndependentDelta(INoOfVariables),STAT=IErr)
    ALLOCATE(IIndependentVariableType(INoOfVariables),STAT=IErr)
    ALLOCATE(IIndependentVariableAtom(INoOfVariables),STAT=IErr)
    IF(l_alert(IErr,"felixrefine","allocate IndependentVariable")) CALL abort
    RIndependentVariable=Rtemp(1:INoOfVariables)
    IIndependentVariableType=Itemp(1:INoOfVariables)
    IIndependentVariableAtom=Itemp2(1:INoOfVariables)
  END IF

  !--------------------------------------------------------------------
  ! allocate, ImageInitialisation, ImageMaskInitialisation
  !--------------------------------------------------------------------
  ALLOCATE(RhklPositions(INhkl,2),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RhklPositions")) CALL abort

  ! creates circular or square image mask depending upon IMaskFLAG and assign 
  ! IPixelLocations ALLOCATED here
  ! Have removed IMaskFLAG
  CALL ImageMaskInitialisation(IErr)
  IF(l_alert(IErr,"felixrefine","ImageMaskInitialisation")) CALL abort

  !--------------------------------------------------------------------
  ! allocate & setup image arrays for pixel-parallel simulations
  !--------------------------------------------------------------------
  ! All the individual calculations go into RSimulatedPatterns later with MPI_GATHERV
  ! NB RSimulatedPatterns is a vector with respect to pixels, not a 2D image
  ALLOCATE(RSimulatedPatterns(INoOfLacbedPatterns,IThicknessCount,IPixelTotal),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RSimulatedPatterns")) CALL abort
  ! Images to match RImageExpi
  ALLOCATE(RImageSimi(2*IPixelCount,2*IPixelCount,INoOfLacbedPatterns,IThicknessCount),&
        STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RImageSimi")) CALL abort

  RSimulatedPatterns = ZERO
  ! The pixels to be calculated by this core  
  ILocalPixelCountMin= (IPixelTotal*(my_rank)/p)+1
  ILocalPixelCountMax= (IPixelTotal*(my_rank+1)/p) 
  ALLOCATE(RIndividualReflections(INoOfLacbedPatterns,IThicknessCount,&
         (ILocalPixelCountMax-ILocalPixelCountMin)+1),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RIndividualReflections")) CALL abort

  ! position of pixels calculated by this core, IDisplacements & ICount are global variables
  ALLOCATE(IDisplacements(p),ICount(p),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate IDisplacements")) CALL abort
  DO ind = 1,p
    IDisplacements(ind) = (IPixelTotal*(ind-1)/p)*INoOfLacbedPatterns*IThicknessCount
    ICount(ind) = (((IPixelTotal*(ind)/p) - (IPixelTotal*(ind-1)/p)))* &
          INoOfLacbedPatterns*IThicknessCount    
  END DO

  !--------------------------------------------------------------------
  ! baseline simulation
  !--------------------------------------------------------------------
  RFigureofMerit=666.666 ! Initial large value, diabolically
  Iter = 0
  ! baseline simulation with timer
  CALL Simulate(IErr)
  IF(l_alert(IErr,"felixrefine","Simulate")) CALL abort
  CALL PrintEndTime(LS,IStartTime2, "Simulation" )

  IF (ISimFLAG.EQ.1) THEN ! Simulation only mode   
    ! simulate multiple thicknesses
    IF(my_rank.EQ.0) THEN
      WRITE(SPrintString,FMT='(A24,I3,A12)') &
        "Writing simulations for ", IThicknessCount," thicknesses"
      CALL message(LS,SPrintString)
      DO ind = 1,IThicknessCount
        CALL WriteIterationOutput(Iter,ind,IErr)
        IF(l_alert(IErr,"felixrefine","WriteIterationOutput")) CALL abort 
      END DO  
    END IF 
   ELSE ! Refinement Mode
    IF(my_rank.EQ.0) THEN!outputs come from core 0 only
      ! Figure of merit is passed back as a global variable
      CALL FigureOfMeritAndThickness(Iter,IThicknessIndex,IErr)
      IF(l_alert(IErr,"felixrefine",&
            "FigureOfMeritAndThickness")) CALL abort 
      CALL message ( LS, "Writing output; baseline simulation" )
      CALL WriteIterationOutput(Iter,IThicknessIndex,IErr)
      IF(l_alert(IErr,"felixrefine","WriteIterationOutput")) CALL abort
    END IF
    
    !===================================== ! Send the fit index to all cores
    CALL MPI_BCAST(RFigureofMerit,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
    !=====================================
  
    !--------------------------------------------------------------------
    ! Iterations and refinement begin
    !--------------------------------------------------------------------

    !--------------------------------------------------------------------
    ! Branch depending upon refinement method
    !--------------------------------------------------------------------
    ! We have INoOfVariables to refine, held in RIndependentVariable(1:INoOfVariables)
    ! Their type is held in IIndependentVariableType(1:INoOfVariables)

    SELECT CASE(IRefineMethodFLAG)

    CASE(1)
      CALL SimplexRefinement
      IF(l_alert(IErr,"felixrefine","SimplexRefinement")) CALL abort 
    
    CASE(2)
      CALL DownhillRefinement
      IF(l_alert(IErr,"felixrefine","DownhillRefinement")) CALL abort 
      
    CASE(3)
      CALL MaxGradientRefinement
      IF(l_alert(IErr,"felixrefine","MaxGradientRefinement")) CALL abort 
      
    CASE(4)
      CALL PairwiseRefinement
      IF(l_alert(IErr,"felixrefine","PairwiseRefinement")) CALL abort 
     
    CASE DEFAULT ! Simulation only, should never happen
      CALL message( LS, "No refinement, simulation only")
       
    END SELECT 

  END IF

  !--------------------------------------------------------------------
  ! deallocate Memory
  !--------------------------------------------------------------------

  DEALLOCATE(CUgMat,STAT=IErr) 
  DEALLOCATE(CUgMatNoAbs,STAT=IErr)
  DEALLOCATE(CUgMatPrime,STAT=IErr)
  DEALLOCATE(ISymmetryRelations,STAT=IErr)
  DEALLOCATE(IEquivalentUgKey,STAT=IErr)
  DEALLOCATE(CUniqueUg,STAT=IErr)
  DEALLOCATE(RIndividualReflections,STAT=IErr)
  DEALLOCATE(IDisplacements,STAT=IErr)
  DEALLOCATE(ICount,STAT=IErr)
  DEALLOCATE(Rhkl,STAT=IErr)
  DEALLOCATE(RgPoolMag,STAT=IErr)
  DEALLOCATE(RgPool,STAT=IErr)
  DEALLOCATE(RgMatrix,STAT=IErr)
  DEALLOCATE(RSimulatedPatterns,STAT=IErr)
  DEALLOCATE(RAtomPosition,STAT=IErr)
  DEALLOCATE(SAtomName,STAT=IErr)
  DEALLOCATE(RIsoDW,STAT=IErr)
  DEALLOCATE(ROccupancy,STAT=IErr)
  DEALLOCATE(IAtomicNumber,STAT=IErr)
  DEALLOCATE(IAnisoDW,STAT=IErr)
  DEALLOCATE(RAtomCoordinate,STAT=IErr)
  DEALLOCATE(CPseudoAtom,STAT=IErr)
  DEALLOCATE(CPseudoScatt,STAT=IErr)
  ! These are global variables, see smodules.f90

  IF (ISimFLAG.EQ.0) THEN
    DEALLOCATE(RImageExpi,STAT=IErr)  
  END IF  
  
  !--------------------------------------------------------------------
  ! finish off
  !--------------------------------------------------------------------

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
  ! these internal SUBROUTINEs use felixrefine's whole variable namespace

  SUBROUTINE abort
    IErr=1
    IF(l_alert(IErr,"felixrefine","ABORTING")) CONTINUE 
    CALL MPI_Abort(MPI_COMM_WORLD,1,IErr)
    STOP
  END SUBROUTINE abort

  !>
  !! Procedure-description: Refinement using the simplex method
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE SimplexRefinement

    !--------------------------------------------------------------------
    ! array allocations & make simplex matrix parallel
    !--------------------------------------------------------------------

    ALLOCATE(RSimplexVariable(INoOfVariables+1,INoOfVariables), STAT=IErr)  
    ALLOCATE(RSimplexFoM(INoOfVariables+1),STAT=IErr)  
    IF(my_rank.EQ.0) THEN !?? Since simplex not random, could be calculated by all cores
      ALLOCATE(ROnes(INoOfVariables+1,INoOfVariables), STAT=IErr) ! matrix of ones
      IF(l_alert(IErr,"SimplexRefinement","allocate ROnes")) RETURN 
      ! matrix of one +/-RSimplexLengthScale
      ALLOCATE(RSimp(INoOfVariables+1,INoOfVariables), STAT=IErr)
      IF(l_alert(IErr,"SimplexRefinement","allocate RSimp")) RETURN 
      ! diagonal matrix of variables as rows
      ALLOCATE(RVarMatrix(INoOfVariables,INoOfVariables), STAT=IErr)
      IF(l_alert(IErr,"SimplexRefinement","allocate RVarMatrix")) RETURN 
      ROnes=ONE
      RSimp=ONE
      RVarMatrix=ZERO
      FORALL(ind = 1:INoOfVariables) RSimp(ind,ind) = -1.0
      RSimp=RSimp*RSimplexLengthScale + ROnes
      FORALL(ind = 1:INoOfVariables) RVarMatrix(ind,ind) = RIndependentVariable(ind)
      RSimplexVariable=MATMUL(RSimp,RVarMatrix)
    END IF
    !===================================== ! send RSimplexVariable OUT to all cores
    CALL MPI_BCAST(RSimplexVariable,(INoOfVariables+1)*(INoOfVariables),&
          MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
    !=====================================

    !--------------------------------------------------------------------
    ! perform initial simplex simulations
    !--------------------------------------------------------------------

    DO ind = 1,(INoOfVariables+1)
      CALL message(LS,"--------------------------------")
      CALL message(LS,no_tag,"Simplex ",ind, " of ", INoOfVariables+1)
      CALL message(LS,"--------------------------------")
      CALL SimulateAndFit(RSimplexVariable(ind,:),Iter,IThicknessIndex,IErr)! Working as iteration 0?
      IF(l_alert(IErr,"SimplexRefinement","SimulateAndFit")) RETURN

      RSimplexFoM(ind)=RFigureofMerit ! RFigureofMerit returned as global variable
    END DO
    
    !--------------------------------------------------------------------
    ! Apply Simplex Method and iterate
    !--------------------------------------------------------------------

    Iter = 1  
    CALL NDimensionalDownhillSimplex(RSimplexVariable,RSimplexFoM,&
          INoOfVariables+1,INoOfVariables,INoOfVariables,&
          RExitCriteria,Iter,RStandardDeviation,RMean,IErr)
    IF(l_alert(IErr,"SimplexRefinement","NDimensionalDownhillSimplex")) RETURN

  END SUBROUTINE SimplexRefinement

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  !>
  !! Procedure-description: Refinement using a maximum gradient method
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  SUBROUTINE DownhillRefinement

    !--------------------------------------------------------------------
    ! allocations & intialise refinement variables
    !--------------------------------------------------------------------
    ALLOCATE(RVar0(INoOfVariables),STAT=IErr)! incoming set of variables
    IF(l_alert(IErr,"DownhillRefinement","allocate RVar0")) RETURN
    ! set of variables to send out for simulations
    ALLOCATE(RCurrentVar(INoOfVariables),STAT=IErr)
    IF(l_alert(IErr,"DownhillRefinement","allocate RCurrentVar")) RETURN
    ALLOCATE(RLastVar(INoOfVariables),STAT=IErr) ! set of variables updated each cycle
    IF(l_alert(IErr,"DownhillRefinement","allocate RLastVar")) RETURN
    ! the vector describing the current line in parameter space
    ALLOCATE(RPVec(INoOfVariables),STAT=IErr)
    IF(l_alert(IErr,"DownhillRefinement","allocate RPVec")) RETURN
    ! the list of fit indices resulting from small changes Rdf for each variable in RPVec
    ALLOCATE(RFitVec(INoOfVariables),STAT=IErr)
    IF(l_alert(IErr,"DownhillRefinement","allocate RFitVec")) RETURN
    
    RBestFit=RFigureofMerit
    RLastFit=RBestFit
    RLastVar=RIndependentVariable
    RCurrentVar=ONE
    Rdf=ONE
    RScale=RSimplexLengthScale
    nnd=0 ! max/min gradient flag

    !--------------------------------------------------------------------
    ! iteratively refine until improvement in fit below exit criteria
    !--------------------------------------------------------------------

    !\/------------------------------------------------------------------
    DO WHILE (Rdf.GE.RExitCriteria)

      RVar0=RIndependentVariable ! incoming point in n-dimensional parameter space
      RFit0=RFigureofMerit ! incoming fit
      !See WriteIterationOutputWrapper for what IPrintFLAG does    
      IPrintFLAG=0

      !--------------------------------------------------------------------
      ! change max gradient vector (RPVec) depending upon max/min gradient situation 
      !--------------------------------------------------------------------      
      
      IF (MOD(nnd,2).EQ.0) THEN ! even number, find the gradient
        DO ind=1,INoOfVariables ! calculate individual gradients
          Iter=Iter+1
          
          ! The type of variable being refined 
          IVariableType=MOD(IIndependentVariableType(ind),10) 
          ! print to screen
          SELECT CASE(IVariableType)
            CASE(1)
              CALL message(LS,"Changing Ug")
            CASE(2)
              CALL message(LS,"Changing atomic coordinate")
            CASE(3)
              CALL message(LS,"Changing occupancy")
            CASE(4)
              CALL message(LS,"Changing isotropic Debye-Waller factor")
            CASE(5)
              CALL message(LS,"Changing anisotropic Debye-Waller factor")
            CASE(6)
              CALL message(LS,"Changing lattice parameter")
            CASE(8)
              CALL message(LS,"Changing convergence angle")
           END SELECT
          IF (RCurrentVar(ind).LE.TINY.AND.IVariableType.EQ.4) THEN ! skip zero DW factor
            RPVec(ind)=TINY
            CYCLE
          END IF
          WRITE(SPrintString,FMT='(A18,I3,A4,I3)') "Finding gradient, ",ind," of ",INoOfVariables
          SPrintString=TRIM(ADJUSTL(SPrintString))
          CALL message(LS,SPrintString)

          ! Rdx is a small change in the current variable determined by RScale
          ! which is either RScale/10 for atomic coordinates and
          ! RScale*variable/10 for everything else 
          IF (my_rank.EQ.0) THEN
            CALL SYSTEM_CLOCK(mnd)
            RandomSign=SIGN(ONE,(REAL(MOD(mnd,10))/TEN)-0.45) ! gives random +/-1
            Rdx=0.1*RandomSign*RScale*RCurrentVar(ind)
            IF(IVariableType.EQ.2) Rdx=0.1*RandomSign*RScale
          END IF
          !=====================================! send Rdx OUT to all cores
          CALL MPI_BCAST(Rdx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
          !=====================================
          RCurrentVar=RVar0
          RCurrentVar(ind)=RCurrentVar(ind)+Rdx
          CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
          IF(l_alert(IErr,"DownhillRefinement","SimulateAndFit")) RETURN
          CALL message ( LS, ".  Writing output; difference images" )
          IF (my_rank.EQ.0) CALL WriteDifferenceImages(Iter,IThicknessIndex,ind,RCurrentVar(ind),Rdx,IErr)
          !If the fit is better, emphasise that parameter *5
          !This gives more importance to parameters that are not stuck in a
          !valley (at the expense of not using the maximum gradient when no
          !parameters are in valleys)
!          IF (RFigureofMerit.LT.RBestFit) THEN
!            REmphasis=TEN/TWO
!          ELSE
            REmphasis=ONE
!          END IF
          CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
          ! BestFitCheck copies RCurrentVar into RIndependentVariable
          ! and updates RBestFit if the fit is better
          CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
          RFitVec(ind)=RFigureofMerit*REmphasis
          RPVec(ind)=(RFit0-RFigureofMerit)/Rdx ! -df/dx: need the dx to keep track of sign
        END DO
      ELSE ! odd number: min gradient - to explore along a valley
        RandomSign=MOD(INT(nnd)+ONE,FOUR)-ONE!ok, it is not random this time!
        DO ind=1,INoOfVariables
!          ! invert gradient
!          IF (ABS(RPVec(ind)).GT.TINY) THEN ! don't invert zeros
!            RPVec(ind)=1/RPVec(ind)
!          ELSE ! keep them as zero
!            RPVec(ind)=ZERO
!          END IF 
          !swap components pairwise to give a zero dot product for even numbers of variables
          IF (MOD(ind,2).EQ.0) THEN!it's an even number, use negative of preceeding variable
            RPVec(ind)=-RLastVar(ind-1)*RandomSign
          ELSEIF (ind.NE.INoOfVariables) THEN
            RPVec(ind)=RLastVar(ind+1)*RandomSign
          END IF
        END DO
        CALL message(LS,"Checking minimum gradient")
      END IF
      !--------------------------------------------------------------------
      ! normalise the max/min gradient vector RPvec & set the first point
      !--------------------------------------------------------------------
      RPvecMag=ZERO
      DO ind=1,INoOfVariables!
        RPvecMag=RPvecMag+RPvec(ind)**2
      END DO
      RPvecMag=SQRT(RPvecMag)
      IF (ABS(RPvecMag).LT.TINY) THEN ! Zero check
        IErr=1
        WRITE(SPrintString,*) RPvec
        IF(l_alert(IErr,"DownhillRefinement",&
              "Chosen variables have no effect! Refinement vector ="//TRIM(SPrintString))) RETURN
      END IF
      IF ((RPvecMag-ONE.EQ.RPvecMag).OR.(RPvecMag.NE.RPvecMag)) THEN ! Infinity and NaN check
        IErr=1
        WRITE(SPrintString,*) RPvec
        IF(l_alert(IErr,"DownhillRefinement",&
              "Infinite or NaN gradient! Refinement vector ="//TRIM(SPrintString))) RETURN
      END IF
      RPvec=RPvec/RPvecMag ! unity vector along direction of max/min gradient
      IF(my_rank.EQ.0) THEN
        WRITE(SPrintString,*) "(A18,",SIZE(RPvec),"(F7.4,1X))"
        WRITE(SPrintString,FMT=SPrintString)"Refinement vector ",RPvec
        SPrintString=TRIM(ADJUSTL(SPrintString))
        CALL message(LS,SPrintString)
      END IF

      RVar0=RIndependentVariable ! the best point of gradient calculation
      RFigureofMerit=RBestFit ! the best fit so far
      !avoid variables that give zero change in fit
      xnd=1!index for the variable to use - we know there is one, otherwise it would have been picked up earlier
      DO WHILE (ABS(RPvec(xnd)).LT.TINY)
        xnd=xnd+1
      END DO!really need to take these variables out of the refinement, but how?
      ! First point, three points to find the minimum
      R3var(1)=RVar0(xnd)! first point is current value
      R3fit(1)=RFigureofMerit! point 1 is the incoming simulation and fit index
      RPvecMag=RVar0(xnd)*RScale ! RPvecMag gives the magnitude of vector in parameter space
      RCurrentVar=RVar0+RPvec*RPvecMag ! Update the parameters to simulate
      
      !--------------------------------------------------------------------
      ! simulate and set the 2nd and 3rd point
      !--------------------------------------------------------------------
      ! Second point
      R3var(2)=RCurrentVar(xnd) 
      CALL message(LS,"Refining, point 2 of 3")
      Iter=Iter+1
      CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
      IF(l_alert(IErr,"DownhillRefinement","SimulateAndFit")) RETURN
      CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
      CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
      R3fit(2)=RFigureofMerit

      ! Third point
      IF (R3fit(2).GT.R3fit(1)) THEN ! new 2 is not better than 1, go the other way
        RPvecMag=-RPvecMag
      ELSE ! it is better, so keep going
        RVar0=RCurrentVar
      END IF
      RCurrentVar=RVar0+RPvec*RPvecMag
      R3var(3)=RCurrentVar(xnd) ! third point
      CALL message( LS, "Refining, point 3 of 3")
      Iter=Iter+1
      CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
      IF(l_alert(IErr,"DownhillRefinement","SimulateAndFit")) RETURN
      CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
      CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
      R3fit(3)=RFigureofMerit

      !--------------------------------------------------------------------
      ! iterate until 3 points concave set
      !--------------------------------------------------------------------

      ! check the three points make a concave set
      jnd=MAXLOC(R3var,1) ! highest x
      knd=MINLOC(R3var,1) ! lowest x
      lnd=6-jnd-knd ! the mid x
      ! Rtest=0.0 would be a straight line, >0=convex, <0=concave
      Rtest=-ABS(R3fit(jnd)-R3fit(knd))
      ! Rconvex is the calculated fit index at the mid x,
      ! if there was a straight line between lowest and highest x
      Rconvex=R3fit(lnd)-(R3fit(knd)+(R3var(lnd)-R3var(knd))*&
        (R3fit(jnd)-R3fit(knd))/(R3var(jnd)-R3var(knd)))

      ! if it isn't more than 10% concave, keep going until it is sufficiently concave
      DO WHILE (Rconvex.GT.0.1*Rtest)
        CALL message( LS, "Convex, continuing")
        jnd=MAXLOC(R3fit,1) ! worst fit
        knd=MINLOC(R3fit,1) ! best fit
        lnd=6-jnd-knd ! the mid fit
        ! replace mid point with a step on from best point
        !?? RB increase the step size by the golden ratio
        !?? RPvecMag=(R3var(knd)-RVar0(1))*(0.5+SQRT(5.0)/2.0)/RPvec(1)
        RPvecMag=TWO*RPvecMag ! double the step size
        ! maximum step in Ug is RMaxUgStep
        IF (ABS(RPvecMag).GT.RMaxUgStep.AND.IRefineMode(1).EQ.1) THEN
          RPvecMag=SIGN(RMaxUgStep,RPvecMag)
        END IF
        RCurrentVar=RVar0+RPvec*RPvecMag
        R3var(lnd)=RCurrentVar(xnd)! next point
        Iter=Iter+1
        CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
        IF(l_alert(IErr,"DownhillRefinement","SimulateAndFit")) RETURN
        CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
        CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
        R3fit(lnd)=RFigureofMerit
        jnd=MAXLOC(R3var,1) ! highest x
        knd=MINLOC(R3var,1) ! lowest x
        lnd=6-jnd-knd ! the mid x
        Rconvex=R3fit(lnd)-(R3fit(knd)+(R3var(lnd)-R3var(knd))*&
              (R3fit(jnd)-R3fit(knd))/(R3var(jnd)-R3var(knd)))
        Rtest=-ABS(R3fit(jnd)-R3fit(knd))
      END DO

      !--------------------------------------------------------------------
      ! make prediction and update last fit
      !--------------------------------------------------------------------

      ! now make a prediction
      CALL Parabo3(R3var,R3fit,RvarMin,RfitMin,IErr)
      IF (my_rank.EQ.0) WRITE(SPrintString,FMT='(A32,F9.4,A16,F10.5)') &
      "Concave set, predict minimum at ",RvarMin," with fit index ",RfitMin
      SPrintString=TRIM(ADJUSTL(SPrintString))
      CALL message (LS, SPrintString)
      RCurrentVar=RVar0+RPvec*(RvarMin-RVar0(xnd))/RPvec(xnd) ! Put prediction into RCurrentVar
      Iter=Iter+1
      CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
      IF(l_alert(IErr,"DownhillRefinement","SimulateAndFit")) RETURN
      IPrintFLAG=1!Save this simulation
      CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
      CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
      ! check where we go next and update last fit etc.
      RLastVar=RPVec
      IF (MOD(nnd,2).EQ.1) THEN ! only update LastFit after a min gradient refinement
        Rdf=RLastFit-RBestFit 
        RLastFit=RBestFit
        CALL message(LS, "--------------------------------")
        WRITE(SPrintString,FMT='(A19,F8.6,A15,F8.6)') "Improvement in fit ",Rdf,", will stop at ",RExitCriteria
        CALL message (LS, SPrintString)
      END IF
      ! shrink length scale as we progress, by a smaller amount
      ! depending on the no of variables: 1->1/2; 2->3/4; 3->5/6; 4->7/8; 5->9/10;
      RScale=RScale*(ONE-ONE/(TWO*REAL(INoOfVariables)))
      nnd=nnd+1 ! do different gradient next time
    END DO
    !/\------------------------------------------------------------------
  
    !--------------------------------------------------------------------
    ! We are done, simulate and output the best fit
    !--------------------------------------------------------------------

    Iter=Iter+1
    CALL SimulateAndFit(RIndependentVariable,Iter,IThicknessIndex,IErr)
    IF(l_alert(IErr,"DownhillRefinement","SimulateAndFit")) RETURN
    IPrintFLAG=2
    CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
    IF(l_alert(IErr,"DownhillRefinement","WriteIterationOutputWrapper")) RETURN
    DEALLOCATE(RVar0,RCurrentVar,RLastVar,RPVec,RFitVec,STAT=IErr)  

  END SUBROUTINE DownhillRefinement

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  !>
  !! Procedure-description: Refinement using a maximum gradient method
  !!
  !! Major-Authors: Richard Beanland (2019)
  !!
  SUBROUTINE MaxGradientRefinement

    !--------------------------------------------------------------------
    ! allocations & intialise refinement variables
    !--------------------------------------------------------------------
    ALLOCATE(RVar0(INoOfVariables),STAT=IErr)! incoming set of variables
    IF(l_alert(IErr,"MaxGradientRefinement","allocate RVar0")) RETURN
    ! set of variables to send out for simulations
    ALLOCATE(RCurrentVar(INoOfVariables),STAT=IErr)
    IF(l_alert(IErr,"MaxGradientRefinement","allocate RCurrentVar")) RETURN
    ALLOCATE(RLastVec(INoOfVariables),STAT=IErr) ! set of variables updated each cycle
    IF(l_alert(IErr,"MaxGradientRefinement","allocate RLastVar")) RETURN
    ! the vector describing the current line in parameter space
    ALLOCATE(RPVec(INoOfVariables),STAT=IErr)
    IF(l_alert(IErr,"MaxGradientRefinement","allocate RPVec")) RETURN
    ! the list of fit indices resulting from small changes Rdf for each variable in RPVec
    
    RBestFit=RFigureofMerit
    RLastFit=RBestFit
    RLastVec=ONE
    RCurrentVar=ONE
    Rdf=ONE
    RScale=RSimplexLengthScale
    RIndependentDelta = ZERO

    !--------------------------------------------------------------------
    ! iteratively refine until improvement in fit below exit criteria
    !--------------------------------------------------------------------

    !\/------------------------------------------------------------------
    DO WHILE (Rdf.GE.RExitCriteria)
      !put current best point in n-dimensional parameter space into workspace
      RCurrentVar=RIndependentVariable
      !running best fit during this refinement cycle goes in RVar0
      RVar0=RIndependentVariable
      RFit0=RFigureofMerit ! incoming fit
      !See WriteIterationOutputWrapper for what IPrintFLAG does    
      IPrintFLAG=0
      !if all parameters have been refined, reset and restart
      IF (SUM(ABS(RLastVec)).LT.TINY) RLastVec=ONE
      DO ind=1,INoOfVariables ! calculate individual gradients
        !skip variables that previous refinements already optimised
        IF (ABS(RLastVec(ind)).LT.TINY) THEN
          RPVec(ind)=ZERO
          CYCLE
        END IF
        Iter=Iter+1
        ! The type of variable being refined 
        IVariableType=MOD(IIndependentVariableType(ind),10) 
        ! print to screen
        SELECT CASE(IVariableType)
          CASE(1)
            CALL message(LS,"Changing Ug")
          CASE(2)
            CALL message(LS,"Changing atomic coordinate")
          CASE(3)
            CALL message(LS,"Changing occupancy")
          CASE(4)
            CALL message(LS,"Changing isotropic Debye-Waller factor")
          CASE(5)
            CALL message(LS,"Changing anisotropic Debye-Waller factor")
          CASE(6)
            CALL message(LS,"Changing lattice parameter")
          CASE(8)
            CALL message(LS,"Changing convergence angle")
        END SELECT
        IF (RCurrentVar(ind).LE.TINY.AND.IVariableType.EQ.4) THEN ! skip zero DW factor
          RPVec(ind)=TINY
          CYCLE
        END IF
        WRITE(SPrintString,FMT='(A18,I3,A4,I3)') "Finding gradient, ",ind," of ",INoOfVariables
        SPrintString=TRIM(ADJUSTL(SPrintString))
        CALL message(LS,SPrintString)
        ! Rdx is a small change in the current variable determined by RScale
        ! which is either RScale for atomic coordinates and
        ! RScale*variable for everything else
        Rdx=ABS(RScale*RCurrentVar(ind))
        IF(IVariableType.EQ.2) Rdx=ABS(RScale)
        ! three point gradient measurement, + first
        RCurrentVar(ind)=RVar0(ind)+Rdx
        CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
        IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
        CALL message ( LS, ".  Writing output; difference images" )
        IF (my_rank.EQ.0) CALL WriteDifferenceImages(Iter,IThicknessIndex,ind,RCurrentVar(ind),Rdx,IErr)
        CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
        ! BestFitCheck copies RCurrentVar into RIndependentVariable
        ! and updates RBestFit if the fit is better
        CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
        RPlus=RFigureofMerit
        ! Now minus dx
        RCurrentVar(ind)=RVar0(ind)-Rdx
        CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
        IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
        CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
        Rminus=RFigureofMerit
        !Reset current point so it's correct for the next calculation
        RCurrentVar(ind)=RVar0(ind)
        !If the three points contain a minimum, predict its position using Kramer's rule
        IF (MIN(RFit0,Rplus,Rminus).EQ.RFit0) THEN
          R3var=(/ (RVar0(ind)-Rdx),RVar0(ind),(RVar0(ind)+Rdx) /)
          R3fit=(/ Rminus,RFit0,Rplus /)
          CALL Parabo3(R3var,R3fit,RvarMin,RfitMin,IErr)
          !error estimate, uses RPrecision from felix.inp
          CALL DeltaX(R3var,R3fit,RdeltaX,RPrecision,IErr)
          RIndependentDelta(ind)=RdeltaX
          !make a string for output
          CALL UncertBrak(RvarMin,RdeltaX,Sest,IErr)
          IF(l_alert(IErr,"MaxGradientRefinement","UncertBrak")) RETURN
          !We update RVar0 with the best points as we go along
          !But keep RCurrentVar the same so that the measurements of gradient
          !are accurate.  
          RVar0(ind)=RvarMin
          RPVec(ind)=ZERO !don't include this variable in the max gradient refinement
          IF (my_rank.EQ.0) WRITE(SPrintString,FMT='(A18,A,A15,F7.3,A1)') &
          "Expect minimum at ",TRIM(ADJUSTL(Sest))," with fit index",(HUNDRED*RfitMin),"%"
          SPrintString=TRIM(ADJUSTL(SPrintString))
          CALL message (LS, SPrintString)
        ELSE!this is a valid gradient descent direction
          RPVec(ind)=-(Rplus-Rminus)/(2*Rdx) ! -df/dx
          IF (MIN(RFit0,Rplus,Rminus).EQ.RPlus) RVar0(ind)=RVar0(ind)+Rdx
          IF (MIN(RFit0,Rplus,Rminus).EQ.Rminus) RVar0(ind)=RVar0(ind)-Rdx
        END IF
      END DO
      !We have not run any simulation for the predicted best point so do it now
      RCurrentVar=RVar0!
      Iter=Iter+1
!      IF (my_rank.EQ.0) THEN
!        WRITE(SPrintString,*) "(A15,",SIZE(RPvec),"(F7.4,1X),A27)"
!        WRITE(SPrintString,FMT=SPrintString) &
!        "First point at ",Rvar0," should have best fit index"
!        SPrintString=TRIM(ADJUSTL(SPrintString))
!        CALL message (LS, SPrintString)
!      END IF
      CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
      IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
      CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
      RFit0=RFigureofMerit ! should be the best fit so far
      !--------------------------------------------------------------------
      ! normalise the max/min gradient vector RPvec & set the first point
      !--------------------------------------------------------------------
      RPvecMag=ZERO
      DO ind=1,INoOfVariables!
        RPvecMag=RPvecMag+RPvec(ind)**2
      END DO
      RPvecMag=SQRT(RPvecMag)
      IF ((RPvecMag-ONE.EQ.RPvecMag).OR.(RPvecMag.NE.RPvecMag)) THEN ! Infinity and NaN check
        IErr=1
        WRITE(SPrintString,*) RPvec
        IF(l_alert(IErr,"MaxGradientRefinement",&
              "Infinite or NaN gradient! Refinement vector ="//TRIM(SPrintString))) RETURN
      END IF
      IF (ABS(RPvecMag).GT.TINY) THEN ! There are non-zero gradients, do the vector descent
        RPvec=RPvec/RPvecMag ! unity vector along direction of max/min gradient
        IF(my_rank.EQ.0) THEN
          WRITE(SPrintString,*) "(A18,",SIZE(RPvec),"(F7.4,1X))"
          WRITE(SPrintString,FMT=SPrintString)"Refinement vector ",RPvec
          SPrintString=TRIM(ADJUSTL(SPrintString))
          CALL message(LS,SPrintString)
        END IF
        !avoid variables that give zero change in fit
        !xnd=index for the variable to use (we know there is one,
        ! otherwise it would have been picked up earlier)
        xnd=1
        DO WHILE (ABS(RPvec(xnd)).LT.TINY)
          xnd=xnd+1
        END DO
        ! First point, three points to find the minimum
        R3var(1)=RVar0(xnd)! first point is current value
        R3fit(1)=RFigureofMerit! point 1 is the incoming simulation and fit index
        RPvecMag=RVar0(xnd)*RScale ! RPvecMag gives the magnitude of vector in parameter space
        RCurrentVar=RVar0+RPvec*RPvecMag ! Update the parameters to simulate

        !--------------------------------------------------------------------
        ! simulate and set the 2nd and 3rd point
        !--------------------------------------------------------------------
        ! Second point
        R3var(2)=RCurrentVar(xnd) 
        CALL message(LS,"Refining, point 2 of 3")
        Iter=Iter+1
        CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
        IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
        CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
        CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
        R3fit(2)=RFigureofMerit
  
        ! Third point
        IF (R3fit(2).GT.R3fit(1)) THEN ! new 2 is not better than 1, go the other way
          RPvecMag=-RPvecMag
        ELSE ! it is better, so keep going; move the running best fit point
          RVar0=RCurrentVar
        END IF
        RCurrentVar=RVar0+RPvec*RPvecMag
        R3var(3)=RCurrentVar(xnd) ! third point
        CALL message( LS, "Refining, point 3 of 3")
        Iter=Iter+1
        CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
        IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
        CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
        CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
        R3fit(3)=RFigureofMerit

        !--------------------------------------------------------------------
        ! iterate until 3 points concave set
        !--------------------------------------------------------------------
 
        ! check the three points make a concave set
        jnd=MAXLOC(R3var,1) ! highest x
        knd=MINLOC(R3var,1) ! lowest x
        IF (jnd.NE.knd) THEN!if j=k there is no effect on fit, so skip to end
          lnd=6-jnd-knd ! the mid x
          ! Rtest=0.0 would be a straight line, >0=convex, <0=concave
          Rtest=-ABS(R3fit(jnd)-R3fit(knd))
          ! Rconvex is the calculated fit index at the mid x,
          ! if there was a straight line between lowest and highest x
          Rconvex=R3fit(lnd)-(R3fit(knd)+(R3var(lnd)-R3var(knd))*&
           (R3fit(jnd)-R3fit(knd))/(R3var(jnd)-R3var(knd)))

          ! if it isn't more than 10% concave, keep going until it is sufficiently concave
          DO WHILE (Rconvex.GT.0.1*Rtest)
            CALL message( LS, "Convex, continuing")
            jnd=MAXLOC(R3fit,1) ! worst fit
            knd=MINLOC(R3fit,1) ! best fit
            IF (jnd.NE.knd) THEN
              lnd=6-jnd-knd ! the mid fit
              ! replace mid point with a step on from best point
              RPvecMag=TWO*RPvecMag ! double the step size
              ! maximum step in Ug is RMaxUgStep
              IF (ABS(RPvecMag).GT.RMaxUgStep.AND.IRefineMode(1).EQ.1) THEN
                RPvecMag=SIGN(RMaxUgStep,RPvecMag)
              END IF
              RCurrentVar=RVar0+RPvec*RPvecMag
              R3var(lnd)=RCurrentVar(xnd)! next point
              Iter=Iter+1
              CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
              IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
              CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
              CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
              R3fit(lnd)=RFigureofMerit
              jnd=MAXLOC(R3var,1) ! highest x
              knd=MINLOC(R3var,1) ! lowest x
              IF (jnd.NE.knd) THEN
                lnd=6-jnd-knd ! the mid x
                Rconvex=R3fit(lnd)-(R3fit(knd)+(R3var(lnd)-R3var(knd))*&
                  (R3fit(jnd)-R3fit(knd))/(R3var(jnd)-R3var(knd)))
                Rtest=-ABS(R3fit(jnd)-R3fit(knd))
              END IF
            END IF
          END DO
        END IF

        !--------------------------------------------------------------------
        ! make prediction and update last fit
        !--------------------------------------------------------------------
  
        ! now make a prediction
        CALL Parabo3(R3var,R3fit,RvarMin,RfitMin,IErr)
        IF (my_rank.EQ.0) WRITE(SPrintString,FMT='(A32,F9.4,A16,F10.5)') &
        "Concave set, predict minimum at ",RvarMin," with fit index ",RfitMin
        SPrintString=TRIM(ADJUSTL(SPrintString))
        CALL message (LS, SPrintString)
        RCurrentVar=RVar0+RPvec*(RvarMin-RVar0(xnd))/RPvec(xnd) ! Put prediction into RCurrentVar
      ELSE
        RCurrentVar=RVar0
      END IF
      Iter=Iter+1
      CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
      IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
      IPrintFLAG=1!Save this simulation
      CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
      CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
      ! check where we go next and update last fit etc.
      RLastVec=RPVec
      Rdf=RLastFit-RBestFit 
      RLastFit=RBestFit
      CALL message(LS, "--------------------------------")
      WRITE(SPrintString,FMT='(A19,F8.6,A15,F8.6)') "Improvement in fit ",Rdf,", will stop at ",RExitCriteria
      CALL message (LS, SPrintString)
      ! shrink length scale as we progress, by a smaller amount
      ! depending on the no of variables: 1->1/2; 2->3/4; 3->5/6; 4->7/8; 5->9/10;
      RScale=RScale*(ONE-ONE/(TWO*REAL(INoOfVariables)))
    END DO
    !/\------------------------------------------------------------------
  
    !--------------------------------------------------------------------
    ! We are done, simulate and output the best fit
    !--------------------------------------------------------------------

!    Iter=Iter+1
!    CALL SimulateAndFit(RIndependentVariable,Iter,IThicknessIndex,IErr)
!    IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
    IPrintFLAG=2
    CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
    IF(l_alert(IErr,"MaxGradientRefinement","WriteIterationOutputWrapper")) RETURN
 !   DEALLOCATE(RVar0,RCurrentVar,RLastVar,RPVec,STAT=IErr)  

  END SUBROUTINE MaxGradientRefinement

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Refinement using the parabola method
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  SUBROUTINE PairwiseRefinement

    !--------------------------------------------------------------------
    ! allocate & intilise refinement variables
    !--------------------------------------------------------------------

    ALLOCATE(RVar0(INoOfVariables),STAT=IErr)! incoming set of variables
    IF(l_alert(IErr,"PairwiseRefinement","allocate RVar0")) RETURN
    ! set of variables to send out for simulations
    ALLOCATE(RCurrentVar(INoOfVariables),STAT=IErr)
    IF(l_alert(IErr,"PairwiseRefinement","allocate RCurrentVar")) RETURN
    ALLOCATE(RLastVar(INoOfVariables),STAT=IErr)! set of variables updated each cycle
    IF(l_alert(IErr,"PairwiseRefinement","allocate RLastVar")) RETURN
    ! the vector describing the current line in parameter space

    ALLOCATE(RPVec(INoOfVariables),STAT=IErr)
    IF(l_alert(IErr,"PairwiseRefinement","allocate RPVec")) RETURN

    RLastFit=RFigureofMerit
    RLastVar=RIndependentVariable
    RBestFit=RFigureofMerit
    Rdf=RFigureofMerit
    ICycle=0 
    RScale=RSimplexLengthScale
    RMaxUgStep=0.005 ! maximum step in Ug is 0.5 nm^-2, 0.005 A^-2

    !--------------------------------------------------------------------
    ! iteratively refine until improvement in fit below exit criteria
    !--------------------------------------------------------------------   

    ! Each complete cycle we will look down the average refinement direction

    !\/\/------------------------------------------------------------------
    DO WHILE (Rdf.GE.RExitCriteria)
      IPrintFLAG=0
      
      !--------------------------------------------------------------------
      ! iterate over each variable to refine
      !--------------------------------------------------------------------

      !\/------------------------------------------------------------------    
      DO ind=1,INoOfVariables

        ! optional terminal output types of variables refined
        IVariableType=MOD(IIndependentVariableType(ind),10)
        SELECT CASE(IVariableType)
          CASE(1)
            CALL message(LS, "Ug refinement")
          CASE(2)
            CALL message(LS, "Atomic coordinate refinement")
          CASE(3)
            CALL message(LS, "Occupancy refinement")
          CASE(4)
            CALL message(LS, "Isotropic Debye-Waller factor refinement")
          CASE(5)
            CALL message(LS, "Convergence angle refinement")
        END SELECT
          IF (INoOfVariables.GT.1) THEN
            IF (ICycle.EQ.1) THEN
              CALL message(LS, "Refining all variables")
            ELSE
              CALL message(LS, "Refining pairs of variables")
            END IF
          END IF

        !--------------------------------------------------------------------
        ! look down average refinement direction or do pair-wise maximum gradient
         RVar0=RIndependentVariable ! incoming point in n-dimensional parameter space
        RPvec=ZERO ! Vector in n-dimensional parameter space for this refinement
        CALL message ( LL, "Current parameters=",RIndependentVariable )

        IF (ICycle.EQ.1) THEN ! look down the average refinement direction
          RPvec=(RCurrentVar-RLastVar)/(RCurrentVar(1)-RLastVar(1))
          RPvecMag=ZERO
          DO jnd=1,INoOfVariables
            RPvecMag=RPvecMag+RPvec(jnd)**2
          END DO
          RPvec=RPvec/SQRT(RPvecMag) ! unity vector
          RLastVar=RCurrentVar

        ELSE IF (INoOfVariables.GT.1) THEN ! pair-wise maximum gradient
          IF (ind.EQ.INoOfVariables) EXIT ! skip the last variable
          CALL message(LS,no_tag,&
                "Finding maximum gradient for variables",ind," and",ind+1)
          ! NB R3fit CONTAINS the three fit indices
          R3fit(1)=RFigureofMerit ! point 1 is the incoming simulation and fit index
          RPvec(ind)=RScale/5.0 ! small change in current variable for second point
          RCurrentVar=RVar0+RPvec ! second point
          CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
          IF(l_alert(IErr,"PairwiseRefinement","SimulateAndFit")) RETURN
          CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
          R3fit(2)=RFigureofMerit
          ! now look along combination of 2 parameters for third point
          RPvec(ind+1)=RScale/5.0 
          RCurrentVar=RVar0+RPvec ! Update the parameters to simulate
          CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
          IF(l_alert(IErr,"PairwiseRefinement","SimulateAndFit")) RETURN
          CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
          R3fit(3)=RFigureofMerit
          ! optimum gradient vector from the three points, magnitude unity
          RPvec=ZERO!reset
          RPvecMag=SQRT((R3fit(1)-R3fit(2))**2+(R3fit(2)-R3fit(3))**2)
          RPvec(ind)=(R3fit(1)-R3fit(2))/RPvecMag
          RPvec(ind+1)=(R3fit(2)-R3fit(3))/RPvecMag
        END IF
        CALL message( LS, "Refinement vector = ",RPvec)

        !--------------------------------------------------------------------
        ! Infinity and NaN check and set 1st point
        !--------------------------------------------------------------------

        ! Infinity and NaN check
        IF (ABS(SUM(RPvec))-1.GT.ABS(SUM(RPvec)).OR.ABS(SUM(RPvec)).NE.ABS(SUM(RPvec))) EXIT
        RVar0=RIndependentVariable ! the best point of the three
        RFigureofMerit=MINVAL(R3fit) ! the best fit of the three
        ! Small DW factor (<0.1) check
        IF (IVariableType.EQ.4.AND.RVar0(ind).LE.0.1) THEN
          CALL message( LS, "Small Debye Waller factor, resetting to 0.1")
          RVar0(ind)=0.1
          RCurrentVar=RVar0
          RPvecMag=RScale
          Iter=Iter+1
          CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
          IF(l_alert(IErr,"PairwiseRefinement","SimulateAndFit")) RETURN
          CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
          IF(l_alert(IErr,"PairwiseRefinement","WriteIterationOutputWrapper")) RETURN
          ! update RIndependentVariable if necessary
          CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
        END IF

        R3var(1)=RVar0(ind) ! first point is current value
        R3fit(1)=RFigureofMerit ! with the current fit index
        RPvecMag=RVar0(ind)*RScale  ! magnitude of vector in parameter space
        RCurrentVar=RVar0+RPvec*RPvecMag ! Update the parameters to simulate
        CALL message ( LL, "point 1",R3var(1) )
        CALL message ( LL, "fit",R3fit(1) )
        CALL message ( LL, "RPvecMag=",RPvecMag )

        !--------------------------------------------------------------------
        ! set 2nd and 3rd point
        !--------------------------------------------------------------------

        ! second point
        R3var(2)=RCurrentVar(ind) 
        CALL message( LS, "simulation 2")
        Iter=Iter+1
        CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
        IF(l_alert(IErr,"PairwiseRefinement","SimulateAndFit")) RETURN
        CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
        IF(l_alert(IErr,"PairwiseRefinement","WriteIterationOutputWrapper")) RETURN
        CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
        R3fit(2)=RFigureofMerit
        !CALL message ( LM, "point 2",R3var(2) )
        !CALL message ( LM, "fit",R3fit(2) )
        !third point
        IF (R3fit(2).GT.R3fit(1)) THEN ! new 2 is not better than 1, go the other way
          RPvecMag=-RPvecMag
          RCurrentVar=RVar0+RPvec*RPvecMag
        ELSE ! it is better, so keep going
          RCurrentVar=RVar0+TWO*RPvec*RPvecMag
        END IF
        
        ! third point
        R3var(3)=RCurrentVar(ind)
        CALL message( LS, "simulation 2") !?? JR should this be 'simulation 3' for 3rd point (remove this comment once changed)
        Iter=Iter+1
        CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
        IF(l_alert(IErr,"PairwiseRefinement","SimulateAndFit")) RETURN
        CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
        IF(l_alert(IErr,"PairwiseRefinement","WriteIterationOutputWrapper")) RETURN
        CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
        R3fit(3)=RFigureofMerit
        !CALL message ( LM, "point 3",R3var(3) )
        !CALL message ( LM, "fit",R3fit(3) )

        !?? repeated code around here like maximum gradient JR

        !--------------------------------------------------------------------
        ! iterate until 3 points concave set
        !--------------------------------------------------------------------

        ! check the three points make a concave set
        jnd=MAXLOC(R3var,1) ! highest x
        knd=MINLOC(R3var,1) ! lowest x
        lnd=6-jnd-knd ! the mid x
        ! covering the unlikely scenario that all simulations have the same fit index
        IF (jnd.EQ.knd) lnd=jnd
        ! Rtest=0.0 would be a straight line, >0=convex, <0=concave
        Rtest=-ABS(R3fit(jnd)-R3fit(knd))
        ! Rconvex is the calculated fit index at the mid x,
        ! if ithere was a straight line between lowest and highest x
        Rconvex=R3fit(lnd)-(R3fit(knd)+(R3var(lnd)-R3var(knd))*&
              (R3fit(jnd)-R3fit(knd))/(R3var(jnd)-R3var(knd)))

        ! if it isn't more than 10% concave, keep going until it is sufficiently concave
        DO WHILE (Rconvex.GT.0.1*Rtest)
          CALL message( LS, "Convex, continuing")
          jnd=MAXLOC(R3fit,1) ! worst fit
          knd=MINLOC(R3fit,1) ! best fit
          lnd=6-jnd-knd ! the mid fit
          ! replace mid point with a step on from best point
          ! increase the step size by the golden ratio
          RPvecMag=(R3var(knd)-RVar0(ind))*(0.5+SQRT(5.0)/2.0)/RPvec(ind)
          ! maximum step in Ug is RMaxUgStep
          IF (ABS(RPvecMag).GT.RMaxUgStep.AND.IRefineMode(1).EQ.1) THEN
            RPvecMag=SIGN(RMaxUgStep,RPvecMag)
          END IF
          
          RCurrentVar=RVar0+RPvec*RPvecMag     
          R3var(lnd)=RCurrentVar(ind)    ! next point
          IF (R3var(lnd).LE.ZERO.AND.IVariableType.EQ.4) THEN ! less than zero DW is requested
            R3var(lnd)=ZERO ! limit it to 0.0
            RCurrentVar=RVar0-RPvec*RVar0(ind) ! and set up for simulation outside the loop
            EXIT
          END IF

          Iter=Iter+1
          CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
          IF(l_alert(IErr,"PairwiseRefinement","SimulateAndFit")) RETURN
          CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
          IF(l_alert(IErr,"PairwiseRefinement","WriteIterationOutputWrapper")) RETURN
          CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
          R3fit(lnd)=RFigureofMerit
          jnd=MAXLOC(R3var,1) ! highest x
          knd=MINLOC(R3var,1) ! lowest x
          lnd=6-jnd-knd ! the mid x
          Rconvex=R3fit(lnd)-(R3fit(knd)+(R3var(lnd)-R3var(knd))*&
          (R3fit(jnd)-R3fit(knd))/(R3var(jnd)-R3var(knd)))
          Rtest=-ABS(R3fit(jnd)-R3fit(knd))  !?? repeated code around here
        END DO

        !--------------------------------------------------------------------
        ! make prediction and finish with this variable
        !--------------------------------------------------------------------
        IPrintFLAG=1
        IF (RCurrentVar(ind).LE.TINY.AND.(IVariableType.EQ.4.OR.IVariableType.EQ.3)) THEN
          ! We reached zero D-W factor in convexity test, skip the prediction
          CALL message( LS, "Using zero Debye Waller factor, simulate and refine next variable" )
          Iter=Iter+1
          CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
          IF(l_alert(IErr,"PairwiseRefinement","SimulateAndFit")) RETURN
          CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
          IF(l_alert(IErr,"PairwiseRefinement","WriteIterationOutputWrapper")) RETURN
          CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
        ELSE
          CALL Parabo3(R3var,R3fit,RvarMin,RfitMin,IErr)
          CALL message( LS, "Concave set, predict minimum at ",RvarMin)
          CALL message( LS, "       with fit index ",RfitMin)
          jnd=MAXLOC(R3fit,1)!worst point
          knd=MINLOC(R3fit,1)!best point
          ! replace worst point with prediction and put into RIndependentVariable
          ! check if we have reached a zero
          IF (RvarMin.LT.ZERO.AND.IVariableType.EQ.4) RvarMin=ZERO  ! DW factor
          IF (RvarMin.LT.ZERO.AND.IVariableType.EQ.3) RvarMin=ZERO  ! occupancy
          RCurrentVar=RVar0+RPvec*(RvarMin-RVar0(ind))/RPvec(ind)
          !R3var(jnd)=RCurrentVar(ind)!do I need this?
          Iter=Iter+1
          CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
          IF(l_alert(IErr,"PairwiseRefinement","SimulateAndFit")) RETURN
          CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
          IF(l_alert(IErr,"PairwiseRefinement","WriteIterationOutputWrapper")) RETURN
          CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
        END IF

        ! we have just done an average direction refinement, quit the NoOfVariables loop
        IF (ICycle.EQ.1) EXIT

      END DO
      !/\------------------------------------------------------------------

      !--------------------------------------------------------------------
      ! check where we go next and update last fit
      !--------------------------------------------------------------------

      ! We have refined all variables, check where we go next and update last fit etc.
      IF (INoOfVariables.GT.1) THEN ! refining multiple variables
        IF (ICycle.EQ.0) THEN ! we haven't done an average direction refinement, set the flag
          ICycle=1
        ELSE ! we just did an average direction refinement
          ! go back to pairwise and update RLastFit and Rdf
          ICycle=0
          Rdf=RLastFit-RBestFit 
          RLastFit=RBestFit
          CALL message(LS, "--------------------------------")
          WRITE(SPrintString,FMT='(A19,F8.6,A15,F8.6)') "Improvement in fit ",Rdf,", will stop at ",RExitCriteria
          CALL message (LS, SPrintString)
        END IF
      ELSE ! refining just one variable
        IF (RBestFit.LT.RLastFit) THEN
          Rdf=RLastFit-RBestFit 
          RLastFit=RBestFit
          CALL message(LS, "--------------------------------")
          WRITE(SPrintString,FMT='(A19,F8.6,A15,F8.6)') "Improvement in fit ",Rdf,", will stop at ",RExitCriteria
          CALL message (LS, SPrintString)
        END IF
      END IF
      ! shrink length scale as we progress, by a smaller amount
      ! depending on the no of variables: 1->1/2; 2->3/4; 3->5/6; 4->7/8; 5->9/10;
      RScale=RScale*(ONE-ONE/(TWO*REAL(INoOfVariables)))

    END DO 
    !/\/\----------------------------------------------------------------

    !--------------------------------------------------------------------
    ! We are done, finally simulate and output the best fit
    !--------------------------------------------------------------------

    IPrintFLAG=2
    Iter=Iter+1
    CALL SimulateAndFit(RIndependentVariable,Iter,IThicknessIndex,IErr)
    IF(l_alert(IErr,"PairwiseRefinement","SimulateAndFit")) RETURN
    CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)
    IF(l_alert(IErr,"PairwiseRefinement","WriteIterationOutputWrapper")) RETURN
    DEALLOCATE(RVar0,RCurrentVar,RLastVar,RPVec,STAT=IErr)
 
  END SUBROUTINE PairwiseRefinement

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Setup atomic vector movements
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  SUBROUTINE SetupAtomMovements(IErr)

    ! called once in felixrefine IF(IRefineMode(2)==1) atom coordinate refinement, B

    INTEGER(IKIND) :: IErr,knd,jnd,ind
    INTEGER(IKIND),DIMENSION(:),ALLOCATABLE :: IDegreesOfFreedom
    REAL(RKIND),DIMENSION(ITHREE,ITHREE) :: RMoveMatrix
    
    ALLOCATE(IDegreesOfFreedom(SIZE(IAtomsToRefine)),STAT=IErr)
    IF(l_alert(IErr,"SetupAtomMovements","allocate IDegreesOfFreedom")) RETURN  
    IDegreesOfFreedom=0

    !Count the degrees of freedom of movement for each atom to be refined    
    DO ind = 1,SIZE(IAtomsToRefine)
      CALL CountAllowedMovements(ISpaceGrp,SWyckoffSymbol(IAtomsToRefine(ind)),&
            IDegreesOfFreedom(ind),IErr)
      IF(l_alert(IErr,"SetupAtomMovements","CountAllowedMovements")) RETURN     
    END DO

    ALLOCATE(IAtomMoveList(SUM(IDegreesOfFreedom)),STAT=IErr)
    IF(l_alert(IErr,"SetupAtomMovements","allocate IAtomMoveList")) RETURN  
    ALLOCATE(RVector(SUM(IDegreesOfFreedom),ITHREE),STAT=IErr)
    IF(l_alert(IErr,"SetupAtomMovements","allocate RVector")) RETURN  
    
    !make a list of vectors and the atoms they move
    knd = 0
    DO ind = 1,SIZE(IAtomsToRefine)
      CALL DetermineAllowedMovements(ISpaceGrp,SWyckoffSymbol(IAtomsToRefine(ind)),&
            RMoveMatrix,IErr)
      IF(l_alert(IErr,"SetupAtomMovements","DetermineAllowedMovements")) RETURN  
      DO jnd = 1,IDegreesOfFreedom(ind)
        knd=knd+1
        RVector(knd,:)=RMoveMatrix(jnd,:)!the movement, global variable
        IAtomMoveList(knd)=IAtomsToRefine(ind)!the atom, global variable
      END DO
    END DO

  END SUBROUTINE SetupAtomMovements     

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                                                    
  !>
  !! Procedure-description:
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  SUBROUTINE BestFitCheck(RFoM,RBest,RCurrent,RIndependentVariable,IErr)
    ! called in felixrefine multiple times, utility with variable refinement and FigureOfMerit
    
    REAL(RKIND) :: RFoM,RBest ! current figure of merit and the best figure of merit
    ! current and best set of variables
    REAL(RKIND),DIMENSION(INoOfVariables) :: RCurrent,RIndependentVariable
    INTEGER(IKIND) :: IErr

    IF (RFoM.LT.RBest) THEN
      RBest=RFoM
      RIndependentVariable=RCurrent
    END IF
    WRITE(SPrintString,FMT='(A29,F7.2,A1)') "Current best figure of merit ",100*RBest,"%"
    SPrintString=TRIM(ADJUSTL(SPrintString))
    CALL message( LS, SPrintString)

          
  END SUBROUTINE BestFitCheck

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Input is a vector Rx with three x-coordinates and 
  !! Ry with three y-coordinates. Output is the x- and y-coordinate of the vertex
  !! of the fitted parabola, Rxv Ryv. Using Cramer's rules to solve the system of
  !! equations to give Ra(x^2)+Rb(x)+Rc=(y)
  !!
  !! Major-Authors: Richard Beanland (2016)
  !! 
  SUBROUTINE Parabo3(Rx,Ry,Rxv,Ryv,IErr)
    ! called twice in felixrefine, utility for MaxGradient/PairwiseRefinement
    
    USE MyNumbers
    REAL(RKIND) :: Ra,Rb,Rc,Rd,Rxv,Ryv
    REAL(RKIND),DIMENSION(3) :: Rx,Ry
    INTEGER(IKIND) :: IErr

    !y=a*x^2+b*x+c a=Ra, b=Rb, c=Rc    
    Rd = Rx(1)*Rx(1)*(Rx(2)-Rx(3)) + Rx(2)*Rx(2)*(Rx(3)-Rx(1)) + Rx(3)*Rx(3)*(Rx(1)-Rx(2))
    IF (Rd.GT.TINY) THEN  ! we get zero Rd if all three inputs are the same
      Ra =(Rx(1)*(Ry(3)-Ry(2)) + Rx(2)*(Ry(1)-Ry(3)) + Rx(3)*(Ry(2)-Ry(1)))/Rd
      Rb =( Rx(1)*Rx(1)*(Ry(2)-Ry(3)) + Rx(2)*Rx(2)*(Ry(3)-Ry(1)) + &
          Rx(3)*Rx(3)*(Ry(1)-Ry(2)) )/Rd
      Rc =(Rx(1)*Rx(1)*(Rx(2)*Ry(3)-Rx(3)*Ry(2)) + Rx(2)*Rx(2)*(Rx(3)*Ry(1)-Rx(1)*Ry(3))&
          +Rx(3)*Rx(3)*(Rx(1)*Ry(2)-Rx(2)*Ry(1)))/Rd
      Rxv = -Rb/(2*Ra);!x-coord
      Ryv = Rc-Rb*Rb/(4*Ra)!y-coord
    ELSE
      Rxv = Rx(1)
      Ryv = Ry(1)
    END IF

  END SUBROUTINE Parabo3 

  !>
  !! Procedure-description: Inputs Rx, Ry and an error estimate Rdy.
  !! Output is the corresponding error Rdx using Cramer's rules to
  !!  give Ra(x^2)+Rb(x)+Rc=(y)
  !!
  !! Major-Authors: Richard Beanland (2019)
  !! 
  SUBROUTINE DeltaX(Rx,Ry,Rdx,Rdy,IErr)

    USE MyNumbers
    REAL(RKIND) :: Ra,Rd,Rdx,Rdy
    REAL(RKIND),DIMENSION(3) :: Rx,Ry
    INTEGER(IKIND) :: IErr

    !y=a*x^2+b*x+c a=Ra, b=Rb, c=Rc    
    Rd = Rx(1)*Rx(1)*(Rx(2)-Rx(3)) + Rx(2)*Rx(2)*(Rx(3)-Rx(1)) + Rx(3)*Rx(3)*(Rx(1)-Rx(2))
    IF (Rd.GT.TINY) THEN  ! we get zero Rd if all three inputs are the same
      Ra =(Rx(1)*(Ry(3)-Ry(2)) + Rx(2)*(Ry(1)-Ry(3)) + Rx(3)*(Ry(2)-Ry(1)))/Rd
      Rdx = 0.5*SQRT(Rdy/Ra)
    ELSE
      Rdx = ZERO
    END IF

  END SUBROUTINE DeltaX


END PROGRAM Felixrefine
