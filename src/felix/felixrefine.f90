!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-18, all rights reserved
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
        INhkl,IExitFLAG,ICycle,INumInitReflections,IZerothLaueZoneLevel,&
        INumFinalReflections,IThicknessIndex,IVariableType,IArrayIndex,&
        IAnisotropicDebyeWallerFactorElementNo
  INTEGER(IKIND) :: IStartTime,IStartTime2
  REAL(RKIND) :: RHOLZAcceptanceAngle,RLaueZoneGz,RMaxGMag,RPvecMag,&
        RScale,RMaxUgStep,Rdx,RStandardDeviation,RMean,RGzUnitVec,RMinLaueZoneValue,&
        Rdf,RLastFit,RBestFit,RMaxLaueZoneValue,RMaxAcceptanceGVecMag,&
        RLaueZoneElectronWaveVectorMag,RvarMin,RfitMin,RFit0,Rconvex,Rtest
  REAL(RKIND),DIMENSION(ITHREE) :: R3var,R3fit
  INTEGER(IKIND),DIMENSION(10) :: INoOfVariablesForRefinementType

  ! allocatable arrays
  INTEGER(IKIND),DIMENSION(:),ALLOCATABLE :: IOriginGVecIdentifier
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RSimplexFoM,RIndependentVariable,&
        RCurrentVar,RVar0,RLastVar,RPvec,RFitVec
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: RSimplexVariable,RgDummyVecMat,&
        RgPoolMagLaue,RTestImage,ROnes,RVarMatrix,RSimp
  
  CHARACTER(40) :: my_rank_string
  CHARACTER(20) :: h,k,l
  CHARACTER(200) :: SPrintString

  !--------------------------------------------------------------------
  ! startup
  !--------------------------------------------------------------------

  ! initialise constants
  CALL Init_Numbers ! constants for calculations
  CALL InitialiseMessage ! constants required for formatted terminal output
  IErr=0
  IInitialSimulationFLAG = 1

  ! MPI initialization
  CALL MPI_Init(IErr) 
  IF(l_alert(IErr,"felixrefine","MPI_Init")) CALL abort
  CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,IErr) ! get rank of the current process
  IF(l_alert(IErr,"felixrefine","MPI_Comm_rank")) CALL abort
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,IErr) ! get size of the current communicator
  IF(l_alert(IErr,"felixrefine","MPI_Comm_size")) CALL abort

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

  CALL ReadInpFile(IErr) ! felix.inp
  IF(l_alert(IErr,"felixrefine","ReadInpFile")) CALL abort
  CALL SetMessageMode( IWriteFLAG, IErr )
  IF(l_alert(IErr,"felixrefine","set_message_mod_mode")) CALL abort

  CALL message(LL,'IBlochMethodFLAG =',IBlochMethodFLAG)

  CALL ReadHklFile(IErr) ! the list of hkl's to input/output
  IF(l_alert(IErr,"felixrefine","ReadHklFile")) CALL abort

  ! read experimental images (if in refinement mode)
  IF (ISimFLAG.EQ.0) THEN ! it's a refinement, so read-IN experimental images
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

  !--------------------------------------------------------------------
  ! allocate atom and Debye-Waller factor arrays
  !--------------------------------------------------------------------
  ! Reset the basis so that atomic coordinate refinement is possible 
  !CALL PreferredBasis(IErr)!A crystallography subroutine
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
  !--------------------------------------------------------------------

  ! fills unit cell from basis and symmetry, removes duplicate atoms at special positions
  CALL UniqueAtomPositions(IErr)
  IF(l_alert(IErr,"felixrefine","UniqueAtomPositions")) CALL abort
  !?? RB could re-allocate RAtomCoordinate,SAtomName,RIsoDW,ROccupancy,
  !?? IAtomicNumber,IAnisoDW to match INAtomsUnitCell?

  RHOLZAcceptanceAngle=TWODEG2RADIAN !?? RB seems way too low?
  IHKLMAXValue = 5 ! starting value, increments in loop below

  !?? RB Note the application of acceptance angle is incorrect since it uses hkl;
  !?? it should use the reciprocal lattice vectors as calculated in RgPool

  ! Count the reflections that make up the pool of g-vectors, simply returns INhkl
  CALL HKLCount(IHKLMAXValue,RZDirC,INhkl,RHOLZAcceptanceAngle,IErr)
  IF(l_alert(IErr,"felixrefine","HKLCount")) CALL abort
  DO WHILE (INhkl.LT.IMinReflectionPool) 
    IHKLMAXValue = IHKLMAXValue+1
    CALL HKLCount(IHKLMAXValue,RZDirC,INhkl,RHOLZAcceptanceAngle,IErr)
  END DO
  
  ! Fill the list of reflections Rhkl
  ! NB Rhkl are in INTEGER form [h,k,l] but are REAL to allow dot products etc.
  ALLOCATE(Rhkl(INhkl,ITHREE),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate Rhkl")) CALL abort
  CALL HKLMake(IHKLMAXValue,RZDirC,RHOLZAcceptanceAngle,IErr)
  IF(l_alert(IErr,"felixrefine","HKLMake")) CALL abort
  CALL message(LL,dbg7,"Rhkl matrix: ",NINT(Rhkl(1:INhkl,:)))

  !--------------------------------------------------------------------
  ! sort HKL, specific reflections
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
  ! allocate RgPool and dummy
  !--------------------------------------------------------------------

  ! RgPool is a list of g-vectors in the microscope ref frame,
  ! units of 1/A (NB exp(-i*q.r),  physics negative convention)
  ALLOCATE(RgPool(INhkl,ITHREE),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RgPool")) CALL abort
  ALLOCATE(RgDummyVecMat(INhkl,ITHREE),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RgDummyVecMat")) CALL abort
 
  !--------------------------------------------------------------------
  ! calculate g vector list, considering zeroth order Laue zone
  !--------------------------------------------------------------------

  ! Calculate the g vector list RgPool in reciprocal angstrom units
  ! in the microscope reference frame
  ICutOff = 1
  DO ind=1,INhkl
    DO jnd=1,ITHREE
      RgPool(ind,jnd) = Rhkl(ind,1)*RarVecM(jnd) + &
            Rhkl(ind,2)*RbrVecM(jnd) + Rhkl(ind,3)*RcrVecM(jnd)
    ENDDO
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
  
  ! sort into Laue Zones
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
  !--------------------------------------------------------------------

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
  ! calculate g vector magnitudes and components parallel to the surface
  !--------------------------------------------------------------------

  ! calculate g-vector magnitudes for the reflection pool RgPoolMag
  ! in reciprocal Angstrom units, in the Microscope reference frame
  ALLOCATE(RgPoolMag(INhkl),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RgPoolMag")) CALL abort
  DO ind=1,INhkl
     RgPoolMag(ind)= SQRT(DOT_PRODUCT(RgPool(ind,:),RgPool(ind,:)))
  END DO
  CALL message(LL,dbg7,"g-vectors and magnitude (1/A), in the microscope reference frame" )
  DO ind = 1,SIZE(Rhkl,1)
    CALL message(LL,dbg7,"hkl  :",NINT(Rhkl(ind,:)))
    CALL message(LL,dbg7,"g mag:",RgPoolMag(ind))
  END DO

  ! g-vector components parallel to the surface unit normal
  ALLOCATE(RgDotNorm(INhkl),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RgDotNorm")) CALL abort

  DO ind =1,INhkl
    RgDotNorm(ind) = DOT_PRODUCT(RgPool(ind,:),RNormDirM)
  END DO
  CALL message(LL,dbg7,"g.n list")
  DO ind = 1,SIZE(Rhkl,1)
    CALL message(LL,dbg7,"hkl :",NINT(Rhkl(ind,:)))
    CALL message(LL,dbg7,"g.n :",RgDotNorm(ind))
  END DO

  !--------------------------------------------------------------------
  ! calculate resolution in k space
  !--------------------------------------------------------------------
  RMinimumGMag = RgPoolMag(2)!since the first one is always 000
  RDeltaK = RMinimumGMag*RConvergenceAngle/REAL(IPixelCount,RKIND)

  !--------------------------------------------------------------------
  ! calculate RMaxGMag using acceptance angle, count reflections
  !--------------------------------------------------------------------
  ! acceptance angle
  IF(RAcceptanceAngle.NE.ZERO.AND.IHOLZFLAG.EQ.0) THEN
    RMaxAcceptanceGVecMag=(RElectronWaveVectorMagnitude*TAN(RAcceptanceAngle*DEG2RADIAN))
    IF(RgPoolMag(IMinReflectionPool).GT.RMaxAcceptanceGVecMag) THEN
      RMaxGMag = RMaxAcceptanceGVecMag 
    ELSE
      RMaxGMag = RgPoolMag(IMinReflectionPool)
    END IF
  ELSEIF(RAcceptanceAngle.NE.ZERO.AND.IHOLZFLAG.EQ.1) THEN
    IBSMaxLocGVecAmp=MAXVAL(IOriginGVecIdentifier)
    RMaxGMag=RgPoolMag(IBSMaxLocGVecAmp)
    IF(RgPoolMag(IBSMaxLocGVecAmp).GE.RgPoolMag(IMinreflectionPool)) THEN
      RMaxGMag = RgPoolMag(IMinReflectionPool)
    END IF
  ELSE
    RMaxGMag = RgPoolMag(IMinReflectionPool)
  END IF
  
  ! count reflections up to cutoff magnitude
  nReflections=0_IKIND
  DO ind=1,INhkl
    IF (ABS(RgPoolMag(ind)).LE.RMaxGMag) nReflections=nReflections+1
  ENDDO
  IF (nReflections.LT.INoOfLacbedPatterns) nReflections = INoOfLacbedPatterns
  
  ! deallocation
  DEALLOCATE(RgPoolMagLaue)!
  IF (RAcceptanceAngle.NE.ZERO.AND.IHOLZFLAG.EQ.1) THEN
    DEALLOCATE(IOriginGVecIdentifier)
  END IF

  !--------------------------------------------------------------------
  ! allocate Ug arrays
  !--------------------------------------------------------------------
  ! Ug matrix etc.
  ALLOCATE(CUgMatNoAbs(nReflections,nReflections),STAT=IErr)! Ug Matrix without absorption
  IF(l_alert(IErr,"felixrefine","allocate CUgMatNoAbs")) CALL abort
  ALLOCATE(CUgMatPrime(nReflections,nReflections),STAT=IErr)! U'g Matrix of just absorption
  IF(l_alert(IErr,"felixrefine","allocate CUgMatPrime")) CALL abort
  ALLOCATE(CUgMat(nReflections,nReflections),STAT=IErr)! Ug+U'g Matrix, including absorption
  IF(l_alert(IErr,"felixrefine","allocate CUgMat")) CALL abort
  ! Matrix of 2pi*g-vectors that corresponds to the CUgMatNoAbs matrix
  ALLOCATE(RgMatrix(nReflections,nReflections,ITHREE),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RgMatrix")) CALL abort
  ! Matrix of their magnitudes 
  ALLOCATE(RgMatrixMagnitude(nReflections,nReflections),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RgMatrixMagnitude")) CALL abort
  ! Matrix of sums of indices - for symmetry equivalence in the Ug matrix
  ALLOCATE(RgSumMat(nReflections,nReflections),STAT=IErr) 
  IF(l_alert(IErr,"felixrefine","allocate RgSumMat")) CALL abort
  ! Matrix with numbers marking equivalent Ug's
  ALLOCATE(ISymmetryRelations(nReflections,nReflections),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate ISymmetryRelations")) CALL abort
  
  !--------------------------------------------------------------------
  ! calculate reflection matrix & initialise structure factors
  !--------------------------------------------------------------------
  ! Calculate matrix  of g-vectors that corresponds to the Ug matrix
  IThicknessCount= NINT((RFinalThickness-RInitialThickness)/RDeltaThickness) + 1
  DO ind=1,nReflections
    DO jnd=1,nReflections
      RgMatrix(ind,jnd,:)= RgPool(ind,:)-RgPool(jnd,:)
      RgMatrixMagnitude(ind,jnd) = & 
           SQRT(DOT_PRODUCT(RgMatrix(ind,jnd,:),RgMatrix(ind,jnd,:)))
    ENDDO
  ENDDO
  CALL message(LL,dbg3,"g-vector magnitude matrix (2pi/A)", RgMatrixMagnitude(1:16,1:8)) 
  CALL message(LXL,dbg3,"first 16 g-vectors", RgMatrix(1:16,1,:)) 

  ! structure factor initialization
  ! Calculate Ug matrix for each entry in CUgMatNoAbs(1:nReflections,1:nReflections)
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
  CALL message( LM, "Initial Ug matrix, with absorption (nm^-2)" )
  DO ind = 1,16
    WRITE(SPrintString,FMT='(3(I2,1X),A2,1X,8(F7.4,1X))') NINT(Rhkl(ind,:)),": ",100*CUgMat(ind,1:4)
    CALL message( LM,dbg3, SPrintString)
  END DO
  IF(l_alert(IErr,"felixrefine","Absorption")) CALL abort
  CALL PrintEndTime(LS,IStartTime2, "Absorption" )
  CALL message(LL,dbg3,"g-vector magnitude matrix (2pi/A)", RgMatrixMagnitude(1:16,1:8)) 
  CALL SYSTEM_CLOCK( IStartTime2 )

  !--------------------------------------------------------------------
  ! INoOfVariables calculated depending upon Ug and non-Ug refinement
  !--------------------------------------------------------------------
  ! Ug refinement is a special case and must be done alone
  ! cannot do any other refinement alongside
  IF(ISimFLAG.EQ.0) THEN
    IF(IRefineMode(1).EQ.1) THEN ! It's a Ug refinement, A

      ! Count the number of Independent Variables
      IUgOffset=1  ! can skip Ug's in refinement using offset, 1 is inner potential
      jnd=1
      DO ind = 1+IUgOffset,INoofUgs+IUgOffset
        IF ( ABS(REAL(CUniqueUg(ind),RKIND)).GE.RTolerance ) jnd=jnd+1
        IF ( ABS(AIMAG(CUniqueUg(ind))).GE.RTolerance ) jnd=jnd+1
      END DO
      IF (IAbsorbFLAG.EQ.1) THEN ! proportional absorption
        INoOfVariables = jnd ! the last variable is for absorption, so included
      ELSE
        INoOfVariables = jnd-1 
      END IF
      IF ( INoOfVariables.EQ.1 ) THEN 
        CALL message(LS,"Only one independent variable")
      ELSE
        CALL message(LS,"number of independent variables = ",INoOfVariables)
      END IF

    ELSE ! It's not a Ug refinement, so count refinement variables
      ! Excluding Ug refinement, various variables can be refined together
      INoOfVariablesForRefinementType(1)=0
      ! Atom coordinate refinement, B
      IF(IRefineMode(2).EQ.1) THEN
        CALL SetupAtomMovements(IErr)
        IF(l_alert(IErr,"felixrefine","SetupAtomMovements")) CALL abort
        INoOfVariablesForRefinementType(2)=IRefineMode(2)*SIZE(IAtomMoveList)
      ELSE
        INoOfVariablesForRefinementType(2)=0
      END IF
      ! Occupancy, C
      INoOfVariablesForRefinementType(3)=IRefineMode(3)*SIZE(IAtomsToRefine)
      ! Isotropic DW, D
      INoOfVariablesForRefinementType(4)=IRefineMode(4)*SIZE(IAtomsToRefine)
      ! Anisotropic DW, E
      INoOfVariablesForRefinementType(5)=IRefineMode(5)*SIZE(IAtomsToRefine)*6
      INoOfVariablesForRefinementType(6)=IRefineMode(6)*3! Unit cell dimensions, F
      INoOfVariablesForRefinementType(7)=IRefineMode(7)*3! Unit cell angles, G
      INoOfVariablesForRefinementType(8)=IRefineMode(8)! Convergence angle, H
      INoOfVariablesForRefinementType(9)=IRefineMode(9)! Percentage Absorption, I
      INoOfVariablesForRefinementType(10)=IRefineMode(10)! kV, J
      ! Total number of independent variables
      INoOfVariables = SUM(INoOfVariablesForRefinementType)
      IF(INoOfVariables.EQ.0) THEN 
        ! there's no refinement requested, say so and quit
        IErr = 1
        IF(l_alert(IErr,"felixrefine",&
              "No refinement variables! Check IRefineModeFLAG in felix.inp. "// &
              "Valid refine modes are A,B,C,D,E,F,G,H,I,J,S")) CALL abort
      END IF
    END IF
    ALLOCATE(RIndependentVariable(INoOfVariables),STAT=IErr) 
    IF(l_alert(IErr,"felixrefine","allocate RIndependentVariable")) CALL abort

    !--------------------------------------------------------------------
    ! assign refinement variables depending upon Ug and non-Ug refinement
    !--------------------------------------------------------------------
    IF(IRefineMode(1).EQ.1) THEN ! It's a Ug refinement, A
      ! Fill up the IndependentVariable list with CUgMatNoAbs components
      jnd=1
      DO ind = 1+IUgOffset,INoofUgs+IUgOffset
        IF ( ABS(REAL(CUniqueUg(ind),RKIND)).GE.RTolerance ) THEN
          RIndependentVariable(jnd) = REAL(CUniqueUg(ind),RKIND)
          jnd=jnd+1
        END IF
        IF ( ABS(AIMAG(CUniqueUg(ind))).GE.RTolerance ) THEN 
          RIndependentVariable(jnd) = AIMAG(CUniqueUg(ind))
          jnd=jnd+1
        END IF
      END DO
      ! Proportional absorption included in structure factor refinement as last variable
      IF (IAbsorbFLAG.EQ.1) RIndependentVariable(jnd) = RAbsorptionPercentage
    ELSE ! It's not a Ug refinement 
      ! Fill up the IndependentVariable list 
      ALLOCATE(RIndependentVariable(INoOfVariables),STAT=IErr)  
      ind=1
      IF(IRefineMode(2).EQ.1) THEN ! Atomic coordinates, B
        DO jnd=1,SIZE(IAtomMoveList)
            RIndependentVariable(ind)=DOT_PRODUCT(RBasisAtomPosition(IAtomMoveList(jnd),:),RVector(jnd,:))
            ind=ind+1
        END DO
      END IF
      IF(IRefineMode(3).EQ.1) THEN ! Occupancy, C
        DO jnd=1,SIZE(IAtomsToRefine)
            RIndependentVariable(ind)=RBasisOccupancy(IAtomsToRefine(jnd))
            ind=ind+1
        END DO
      END IF
      IF(IRefineMode(4).EQ.1) THEN ! Isotropic DW, D
        DO jnd=1,SIZE(IAtomsToRefine)
            RIndependentVariable(ind)=RBasisIsoDW(IAtomsToRefine(jnd))
            ind=ind+1
        END DO
      END IF
      IF(IRefineMode(8).EQ.1) THEN ! Convergence angle, H
        RIndependentVariable(ind)=RConvergenceAngle
        ind=ind+1
      END IF
      ! Assign IDs - not needed for a Ug refinement
      ALLOCATE(IIterativeVariableUniqueIDs(INoOfVariables,2),STAT=IErr)
      IF(l_alert(IErr,"felixrefine","allocate IIterativeVariableUniqueIDs")) CALL abort
      IIterativeVariableUniqueIDs = 0 
      DO ind = 2,IRefinementVariableTypes ! Loop over iterative variables apart from Ug's
        IF(IRefineMode(ind).EQ.1) THEN
          DO jnd = 1,INoOfVariablesForRefinementType(ind)
            IArrayIndex = SUM(INoOfVariablesForRefinementType(:(ind-1))) + jnd

            SELECT CASE(ind)

            CASE(1) ! Ugs, A
              IIterativeVariableUniqueIDs(IArrayIndex,1) = ind
              IIterativeVariableUniqueIDs(IArrayIndex,2) = &
                    NINT(REAL(INoofUgs,RKIND)*(REAL(jnd/REAL(INoofUgs,RKIND),RKIND)-&
                    CEILING(REAL(jnd/REAL(INoofUgs,RKIND),RKIND)))+REAL(INoofUgs,RKIND))

            CASE(2) ! Coordinates (x,y,z), B
              IIterativeVariableUniqueIDs(IArrayIndex,1) = ind
              IIterativeVariableUniqueIDs(IArrayIndex,2) = jnd

            CASE(3) ! Occupancies, C
              IIterativeVariableUniqueIDs(IArrayIndex,1) = ind
              IIterativeVariableUniqueIDs(IArrayIndex,2) = IAtomsToRefine(jnd)

            CASE(4) ! Isotropic Debye Waller Factors , D
              IIterativeVariableUniqueIDs(IArrayIndex,1) = ind
              IIterativeVariableUniqueIDs(IArrayIndex,2) = IAtomsToRefine(jnd)

            CASE(5) ! Anisotropic Debye Waller Factors (a11-a33), E
              !?? not currently implemented
              IErr=1;IF(l_alert(IErr,"felixrefine",&
                    "Anisotropic Debye Waller Factors not implemented")) CALL abort

  !            IIterativeVariableUniqueIDs(IArrayIndex,1) = ind
  !            IIterativeVariableUniqueIDs(IArrayIndex,2) = &
  !                  IAtomsToRefine(INT(CEILING(REAL(jnd/6.0D0,RKIND))))
  !            IAnisotropicDebyeWallerFactorElementNo = &
  !                  NINT(6.D0*(REAL(jnd/6.0D0,RKIND) - &
  !                  CEILING(REAL(jnd/6.0D0,RKIND)))+6.0D0)

  !            SELECT CASE(IAnisotropicDebyeWallerFactorElementNo)
  !            CASE(1)
  !              IIterativeVariableUniqueIDs(IArrayIndex,4:5) = [1,1]
  !            CASE(2)
  !              IIterativeVariableUniqueIDs(IArrayIndex,4:5) = [2,1]
  !            CASE(3)
  !              IIterativeVariableUniqueIDs(IArrayIndex,4:5) = [2,2]
  !            CASE(4)
  !              IIterativeVariableUniqueIDs(IArrayIndex,4:5) = [3,1]
  !            CASE(5)
  !              IIterativeVariableUniqueIDs(IArrayIndex,4:5) = [3,2]
  !            CASE(6)
  !              IIterativeVariableUniqueIDs(IArrayIndex,4:5) = [3,3]
  !            END SELECT

            CASE(6) ! Lattice Parameters, E
              IIterativeVariableUniqueIDs(IArrayIndex,1) = ind
              IIterativeVariableUniqueIDs(IArrayIndex,2) = &
                    NINT(3.D0*(REAL(jnd/3.0D0,RKIND) - CEILING(REAL(jnd/3.0D0,RKIND)))+3.0D0)
               
            CASE(7) ! Lattice Angles, F
              IIterativeVariableUniqueIDs(IArrayIndex,1) = ind
              IIterativeVariableUniqueIDs(IArrayIndex,2) = &
                    NINT(3.D0*(REAL(jnd/3.0D0,RKIND) - CEILING(REAL(jnd/3.0D0,RKIND)))+3.0D0)

            CASE(8) ! Convergence angle, H
              IIterativeVariableUniqueIDs(IArrayIndex,1) = ind

            CASE(9) ! Percentage Absorption, I
              IIterativeVariableUniqueIDs(IArrayIndex,1) = ind

            CASE(10) ! kV, J
              IIterativeVariableUniqueIDs(IArrayIndex,1) = ind

            END SELECT

          END DO
        END IF
      END DO 

    END IF
  END IF

  !--------------------------------------------------------------------
  ! allocate, ImageInitialisation, ImageMaskInitialisation
  !--------------------------------------------------------------------
  ALLOCATE(RhklPositions(nReflections,2),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RhklPositions")) CALL abort
  ALLOCATE(IMask(2*IPixelCount,2*IPixelCount),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RhklPositions")) CALL abort

  CALL ImageInitialisation( IErr )
  IF(l_alert(IErr,"felixrefine","ImageInitialisation")) CALL abort

  ! creates circular or square image mask depending upon IMaskFLAG and assign 
  ! IPixelLocations ALLOCATED here
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

  IF (ICorrelationFLAG.EQ.3) THEN ! allocate images for masked correlation
    ! Baseline Images to calculate mask
    ALLOCATE(RImageBase(2*IPixelCount,2*IPixelCount,INoOfLacbedPatterns,IThicknessCount),&
          STAT=IErr)
    IF(l_alert(IErr,"felixrefine","allocate RImageBase")) CALL abort
    ! Average Images to calculate mask
    ALLOCATE(RImageAvi(2*IPixelCount,2*IPixelCount,INoOfLacbedPatterns,IThicknessCount),&
          STAT=IErr)
    IF(l_alert(IErr,"felixrefine","allocate RImageAvi")) CALL abort
    RImageAvi=ZERO
    ! Mask Images
    ALLOCATE(RImageMask(2*IPixelCount,2*IPixelCount,INoOfLacbedPatterns),STAT=IErr)
    IF(l_alert(IErr,"felixrefine","allocate RImageMask")) CALL abort
  END IF

  RSimulatedPatterns = ZERO
  ! The pixels to be calculated by this core  
  ILocalPixelCountMin= (IPixelTotal*(my_rank)/p)+1
  ILocalPixelCountMax= (IPixelTotal*(my_rank+1)/p) 
  ALLOCATE(RIndividualReflections(INoOfLacbedPatterns,IThicknessCount,&
         (ILocalPixelCountMax-ILocalPixelCountMin)+1),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RIndividualReflections")) CALL abort

  ! position of pixels calculated by this core, IDisplacements & ICount are global variables
  ALLOCATE(IDisplacements(p),ICount(p),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RhklPositions")) CALL abort
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
  IExitFLAG = 0 ! Do not exit
  IPreviousPrintedIteration = 0  ! ensures baseline simulation is printed
  IF (ISimFLAG.EQ.1) THEN ! Simulation only mode   
  
    !--------------------------------------------------------------------
    ! simulation only mode
    !--------------------------------------------------------------------
    
    ! simulate multiple thicknesses
    IF(my_rank.EQ.0) THEN
      CALL message(LS,"Writing simulations with the number of thicknesses =", IThicknessCount)
      DO IThicknessIndex = 1,IThicknessCount
        CALL WriteIterationOutput(Iter,IThicknessIndex,IExitFLAG,IErr)
        IF(l_alert(IErr,"felixrefine","WriteIterationOutput")) CALL abort 
      END DO  
    END IF 
  
  ELSE ! Refinement Mode
    IF(my_rank.EQ.0) THEN!outputs to disc come from core 0 only
      ! Figure of merit is passed back as a global variable
      CALL FigureOfMeritAndThickness(Iter,IThicknessIndex,IErr)
      IF(l_alert(IErr,"felixrefine",&
            "FigureOfMeritAndThickness")) CALL abort 
      ! Keep baseline simulation for masked correlation
      IF (ICorrelationFLAG.EQ.3) THEN
        RImageBase=RImageSimi  
        RImageAvi=RImageSimi 
      END IF
      CALL message ( LS, "Writing output; baseline simulation" )
      CALL WriteIterationOutput(Iter,IThicknessIndex,IExitFLAG,IErr)
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
    ! For single variables, their type is held in 
    ! IIterativeVariableUniqueIDs(1:INoOfVariables,2)

    SELECT CASE(IRefineMethodFLAG)

    CASE(1)
      CALL SimplexRefinement
      IF(l_alert(IErr,"felixrefine","SimplexRefinement")) CALL abort 
    
    CASE(2)
      CALL MaxGradientRefinement
      IF(l_alert(IErr,"felixrefine","MaxGradientRefinement")) CALL abort 
      
    CASE(3)
      CALL ParabolicRefinement
      IF(l_alert(IErr,"felixrefine","ParabolicRefinement")) CALL abort 
     
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
  DEALLOCATE(RgSumMat,STAT=IErr)
  DEALLOCATE(RSimulatedPatterns,STAT=IErr)
  DEALLOCATE(RAtomPosition,STAT=IErr)
  DEALLOCATE(SAtomName,STAT=IErr)
  DEALLOCATE(RIsoDW,STAT=IErr)
  DEALLOCATE(ROccupancy,STAT=IErr)
  DEALLOCATE(IAtomicNumber,STAT=IErr)
  DEALLOCATE(IAnisoDW,STAT=IErr)
  DEALLOCATE(RAtomCoordinate,STAT=IErr)
  DEALLOCATE(RgMatrixMagnitude,STAT=IErr)
  DEALLOCATE(CPseudoAtom,STAT=IErr)
  DEALLOCATE(CPseudoScatt,STAT=IErr)
  ! These are global variables, see smodules.f90

  IF (IRefineMode(1).EQ.0) THEN
    DEALLOCATE(IIterativeVariableUniqueIDs,STAT=IErr)
  END IF
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
      ! For masked correlation, add to 'average' (extreme difference from baseline)?
      IF(my_rank.EQ.0.AND.ICorrelationFLAG.EQ.3) THEN
        ! replace pixels that are the most different from baseline
        WHERE (ABS(RImageSimi-RImageBase).GT.ABS(RImageAvi-RImageBase))
          RImageAvi=RImageSimi
        END WHERE
      END IF
    END DO

    !--------------------------------------------------------------------
    ! set up masked fitting using the simplex setup
    !--------------------------------------------------------------------

    IF (my_rank.EQ.0.AND.ICorrelationFLAG.EQ.3) THEN
      ! Simple start, just take thickness 1 
      RImageMask=RImageAvi(:,:,:,1)-RImageBase(:,:,:,1)
      DO ind = 1,INoOfLacbedPatterns ! mask each pattern individually
        WHERE (ABS(RImageMask(:,:,ind)).GT.0.1*MAXVAL(ABS(RImageMask(:,:,ind))))
          RImageMask(:,:,ind)=ONE
        ELSEWHERE
          RImageMask(:,:,ind)=ZERO
        END WHERE
      END DO    
      ! debug mode output to look at the masks
      IF (dbg6%LState) THEN
        ALLOCATE(RTestImage(2*IPixelCount,2*IPixelCount),STAT=IErr)
        DO ind = 1,INoOfLacbedPatterns
          RTestImage=RImageMask(:,:,ind)
          WRITE(h,*)  NINT(Rhkl(IOutPutReflections(ind),1))
          WRITE(k,*)  NINT(Rhkl(IOutPutReflections(ind),2))
          WRITE(l,*)  NINT(Rhkl(IOutPutReflections(ind),3))
          WRITE(SPrintString,*) TRIM(ADJUSTL(h)),TRIM(ADJUSTL(k)),TRIM(ADJUSTL(l)),".mask"
          OPEN(UNIT=IChOutWIImage, STATUS= 'UNKNOWN', FILE=SPrintString,&
                FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=2*IPixelCount*8)
          IF(l_alert(IErr,"SimplexRefinement","writing .mask file")) RETURN
          DO jnd = 1,2*IPixelCount
            WRITE(IChOutWIImage,rec=jnd) RTestImage(jnd,:)
          END DO
          CLOSE(IChOutWIImage,IOSTAT=IErr)
        END DO
        DEALLOCATE(RTestImage)
      END IF
    END IF
    
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
  !! Procedure-description: Refinement using the maximum gradient method
  !!
  !! Major-Authors: Richard Beanland (2016)
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
    ALLOCATE(RLastVar(INoOfVariables),STAT=IErr) ! set of variables updated each cycle
    IF(l_alert(IErr,"MaxGradientRefinement","allocate RLastVar")) RETURN
    ! the vector describing the current line in parameter space
    ALLOCATE(RPVec(INoOfVariables),STAT=IErr)
    IF(l_alert(IErr,"MaxGradientRefinement","allocate RPVec")) RETURN
    ! the list of fit indices resulting from small changes Rdf for each variable in RPVec
    ALLOCATE(RFitVec(INoOfVariables),STAT=IErr)
    IF(l_alert(IErr,"MaxGradientRefinement","allocate RFitVec")) RETURN
    
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

      !--------------------------------------------------------------------
      ! change max gradient vector (RPVec) depending upon max/min gradient situation 
      !--------------------------------------------------------------------      
      
      IF (nnd.EQ.0) THEN ! max gradient
        DO ind=1,INoOfVariables ! calculate individual gradients
          ! The type of variable being refined 
          IVariableType=IIterativeVariableUniqueIDs(ind,1) 
          ! variable type as in what refinement mode/variables, 'type' used throughout
           
          ! print to screen
          SELECT CASE(IVariableType)
            CASE(1)
              CALL message(LS,"Ug refinement")
            CASE(2)
              CALL message(LS,"Atomic coordinate refinement")
            CASE(3)
              CALL message(LS,"Occupancy refinement")
            CASE(4)
              CALL message(LS,"Isotropic Debye-Waller factor refinement")
            CASE(5)
              CALL message(LS,"Convergence angle refinement")
          END SELECT
          IF (RCurrentVar(ind).LE.TINY.AND.IVariableType.EQ.4) THEN ! skip zero DW factor
            RPVec(ind)=TINY
            CYCLE
          END IF
          WRITE(SPrintString,FMT='(A18,I3,A4,I3)') "Finding gradient, ",ind," of ",INoOfVariables
          SPrintString=TRIM(ADJUSTL(SPrintString))
          CALL message(LS,SPrintString)

          ! Make a random number to vary the sign of dx, using system clock
          IF (my_rank.EQ.0) THEN
            CALL SYSTEM_CLOCK(mnd)
            Rdx=(REAL(MOD(mnd,10))/TEN)-0.45 ! numbers 0-4 give minus, 5-9 give plus
            Rdx=0.1*Rdx*RScale/ABS(Rdx) ! small change in current variable (RScale/10)is dx
          END IF
          !=====================================! send Rdx OUT to all cores
          CALL MPI_BCAST(Rdx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
          !=====================================
          RCurrentVar=RVar0
          RCurrentVar(ind)=RCurrentVar(ind)+Rdx
          CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
          IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
          ! Do not increment iteration here nor write iteration output
          Iter=Iter+1
          CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
          ! BestFitCheck copies RCurrentVar into RIndependentVariable
          ! and updates RBestFit if the fit is better
          CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
          RFitVec(ind)=RFigureofMerit
          RPVec(ind)=(RFit0-RFigureofMerit)/Rdx ! -df/dx: need the dx to keep track of sign
        END DO
        nnd=1 ! do min gradient next time

      ELSE ! min gradient - to explore along a valley
        DO ind=1,INoOfVariables
          ! invert gradient
          IF (ABS(RPVec(ind)).GT.TINY) THEN ! don't invert zeros
            RPVec(ind)=1/RPVec(ind)
          ELSE ! just make them quite big
            RPVec(ind)=TEN
          END IF
        END DO
        CALL message(LS,"Checking minimum gradient")
        nnd=0 ! do max gradient next time
      END IF
      
      !--------------------------------------------------------------------
      ! normalise the max gradient vector RPvec & set the first point
      !--------------------------------------------------------------------
      RPvecMag=ZERO
      DO ind=1,INoOfVariables
        RPvecMag=RPvecMag+RPvec(ind)**2
      END DO
      IF ((RPvecMag-ONE.EQ.RPvecMag).OR.(RPvecMag.NE.RPvecMag)) THEN ! Infinity and NaN check
        IErr=1
        WRITE(SPrintString,*) RPvec
        IF(l_alert(IErr,"MaxGradientRefinement",&
              "Infinity or NaN error, refinement vector ="//TRIM(SPrintString))) RETURN
      END IF
      RPvec=RPvec/SQRT(RPvecMag) ! unity vector along direction of max gradient
      IF(my_rank.EQ.0) THEN
        WRITE(SPrintString,*) "(A18,",SIZE(RPvec),"(F7.4,1X))"
        WRITE(SPrintString,FMT=SPrintString)"Refinement vector ",RPvec
        SPrintString=TRIM(ADJUSTL(SPrintString))
        CALL message(LS,SPrintString)
      END IF

      RVar0=RIndependentVariable ! the best point of gradient calculation
      RFigureofMerit=RBestFit ! the best fit so far
      ! First point, three points to find the minimum
      R3var(1)=RVar0(1)! first point is current value
      R3fit(1)=RFigureofMerit! point 1 is the incoming simulation and fit index
      RPvecMag=RVar0(1)*RScale ! RPvecMag gives the magnitude of vector in parameter space
      RCurrentVar=RVar0+RPvec*RPvecMag ! Update the parameters to simulate
      
      !--------------------------------------------------------------------
      ! simulate and set the 2nd and 3rd point
      !--------------------------------------------------------------------
      ! Second point
      R3var(2)=RCurrentVar(1) 
      CALL message(LS,"Refining, point 2 of 3")
      Iter=Iter+1
      CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
      IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
      CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
      CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
      R3fit(2)=RFigureofMerit

      ! Third point
      IF (R3fit(2).GT.R3fit(1)) THEN ! new 2 is not better than 1, go the other way
        RPvecMag=-RPvecMag
      ELSE ! it is better, so keep going
        RVar0=RCurrentVar
      END IF
      RCurrentVar=RVar0+RPvec*RPvecMag
      R3var(3)=RCurrentVar(1) ! third point
      CALL message( LS, "Refining, point 3 of 3")
      Iter=Iter+1
      CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
      IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
      CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
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
        !?? RB increase the step size by the golden ratio
        !?? RPvecMag=(R3var(knd)-RVar0(1))*(0.5+SQRT(5.0)/2.0)/RPvec(1)
        RPvecMag=TWO*RPvecMag ! double the step size
        ! maximum step in Ug is RMaxUgStep
        IF (ABS(RPvecMag).GT.RMaxUgStep.AND.IRefineMode(1).EQ.1) THEN
          RPvecMag=SIGN(RMaxUgStep,RPvecMag)
        END IF
        RCurrentVar=RVar0+RPvec*RPvecMag
        R3var(lnd)=RCurrentVar(1)! next point
        Iter=Iter+1
        CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
        IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
        CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
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
      WRITE(SPrintString,FMT='(A32,F7.4,A16,F7.4)') &
      "Concave set, predict minimum at ",RvarMin," with fit index ",RfitMin
      SPrintString=TRIM(ADJUSTL(SPrintString))
      CALL message (LS, SPrintString)
      RCurrentVar=RVar0+RPvec*(RvarMin-RVar0(1))/RPvec(1) ! Put prediction into RCurrentVar
      Iter=Iter+1
      CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
      IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
      CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
      CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)

      ! check where we go next and update last fit etc.
      IF (nnd.EQ.0) THEN ! only update LastFit after a max gradient refinement
        Rdf=RLastFit-RBestFit 
        RLastFit=RBestFit
        CALL message(LS, "--------------------------------")
        WRITE(SPrintString,FMT='(A19,F8.6,A15,F8.6)') "Improvement in fit ",Rdf,", will stop at ",RExitCriteria
        CALL message (LS, SPrintString)
      END IF
      ! shrink length scale as we progress, by a smaller amount
      ! depending on the no of variables: 1->1/2; 2->3/4; 3->5/6; 4->7/8; 5->9/10;
      RScale=RScale*(ONE-ONE/(TWO*REAL(INoOfVariables)))

    END DO
    !/\------------------------------------------------------------------
  
    !--------------------------------------------------------------------
    ! We are done, simulate and output the best fit
    !--------------------------------------------------------------------

    IExitFLAG=1
    Iter=Iter+1
    CALL SimulateAndFit(RIndependentVariable,Iter,IThicknessIndex,IErr)
    IF(l_alert(IErr,"MaxGradientRefinement","SimulateAndFit")) RETURN
    CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
    IF(l_alert(IErr,"MaxGradientRefinement","WriteIterationOutputWrapper")) RETURN
    DEALLOCATE(RVar0,RCurrentVar,RLastVar,RPVec,RFitVec,STAT=IErr)  

  END SUBROUTINE MaxGradientRefinement

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Refinement using the parabola method
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  SUBROUTINE ParabolicRefinement

    !--------------------------------------------------------------------
    ! allocate & intilise refinement variables
    !--------------------------------------------------------------------

    ALLOCATE(RVar0(INoOfVariables),STAT=IErr)! incoming set of variables
    IF(l_alert(IErr,"ParabolicRefinement","allocate RVar0")) RETURN
    ! set of variables to send out for simulations
    ALLOCATE(RCurrentVar(INoOfVariables),STAT=IErr)
    IF(l_alert(IErr,"ParabolicRefinement","allocate RCurrentVar")) RETURN
    ALLOCATE(RLastVar(INoOfVariables),STAT=IErr)! set of variables updated each cycle
    IF(l_alert(IErr,"ParabolicRefinement","allocate RLastVar")) RETURN
    ! the vector describing the current line in parameter space

    ALLOCATE(RPVec(INoOfVariables),STAT=IErr)
    IF(l_alert(IErr,"ParabolicRefinement","allocate RPVec")) RETURN

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
      
      !--------------------------------------------------------------------
      ! iterate over each variable to refine
      !--------------------------------------------------------------------

      !\/------------------------------------------------------------------    
      DO ind=1,INoOfVariables

        ! optional terminal output types of variables refined
        IVariableType=IIterativeVariableUniqueIDs(ind,1)
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
          IF(l_alert(IErr,"ParabolicRefinement","SimulateAndFit")) RETURN
          CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
          R3fit(2)=RFigureofMerit
          ! now look along combination of 2 parameters for third point
          RPvec(ind+1)=RScale/5.0 
          RCurrentVar=RVar0+RPvec ! Update the parameters to simulate
          CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
          IF(l_alert(IErr,"ParabolicRefinement","SimulateAndFit")) RETURN
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
          IF(l_alert(IErr,"ParabolicRefinement","SimulateAndFit")) RETURN
          CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
          IF(l_alert(IErr,"ParabolicRefinement","WriteIterationOutputWrapper")) RETURN
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
        IF(l_alert(IErr,"ParabolicRefinement","SimulateAndFit")) RETURN
        CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
        IF(l_alert(IErr,"ParabolicRefinement","WriteIterationOutputWrapper")) RETURN
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
        IF(l_alert(IErr,"ParabolicRefinement","SimulateAndFit")) RETURN
        CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
        IF(l_alert(IErr,"ParabolicRefinement","WriteIterationOutputWrapper")) RETURN
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
          IF(l_alert(IErr,"ParabolicRefinement","SimulateAndFit")) RETURN
          CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
          IF(l_alert(IErr,"ParabolicRefinement","WriteIterationOutputWrapper")) RETURN
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

        IF (RCurrentVar(ind).LE.TINY.AND.IVariableType.EQ.4) THEN
          ! We reached zero D-W factor in convexity test, skip the prediction
          CALL message( LS, "Using zero Debye Waller factor, simulate and refine next variable" )
          Iter=Iter+1
          CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
          IF(l_alert(IErr,"ParabolicRefinement","SimulateAndFit")) RETURN
          CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
          IF(l_alert(IErr,"ParabolicRefinement","WriteIterationOutputWrapper")) RETURN
          CALL BestFitCheck(RFigureofMerit,RBestFit,RCurrentVar,RIndependentVariable,IErr)
        ELSE
          CALL Parabo3(R3var,R3fit,RvarMin,RfitMin,IErr)
          CALL message( LS, "Concave set, predict minimum at ",RvarMin)
          CALL message( LS, "       with fit index ",RfitMin)
          jnd=MAXLOC(R3fit,1)!worst point
          knd=MINLOC(R3fit,1)!best point
          ! replace worst point with prediction and put into RIndependentVariable
          ! checnk if we have reached zero D-W factor
          IF (RvarMin.LT.ZERO.AND.IVariableType.EQ.4) RvarMin=ZERO
          RCurrentVar=RVar0+RPvec*(RvarMin-RVar0(ind))/RPvec(ind)
          !R3var(jnd)=RCurrentVar(ind)!do I need this?
          Iter=Iter+1
          CALL SimulateAndFit(RCurrentVar,Iter,IThicknessIndex,IErr)
          IF(l_alert(IErr,"ParabolicRefinement","SimulateAndFit")) RETURN
          CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
          IF(l_alert(IErr,"ParabolicRefinement","WriteIterationOutputWrapper")) RETURN
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
      IF (INoOfVariables.GT.2) THEN ! refining multiple variables
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

    IExitFLAG=1
    Iter=Iter+1
    CALL SimulateAndFit(RIndependentVariable,Iter,IThicknessIndex,IErr)
    IF(l_alert(IErr,"ParabolicRefinement","SimulateAndFit")) RETURN
    CALL WriteIterationOutputWrapper(Iter,IThicknessIndex,IExitFLAG,IErr)
    IF(l_alert(IErr,"ParabolicRefinement","WriteIterationOutputWrapper")) RETURN
    DEALLOCATE(RVar0,RCurrentVar,RLastVar,RPVec,STAT=IErr)
 
  END SUBROUTINE ParabolicRefinement

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Setup atomic vector movements
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  SUBROUTINE SetupAtomMovements(IErr)

    ! called once in felixrefine IF(IRefineMode(2)==1) atom coordinate refinement, B

    INTEGER(IKIND) :: IErr,knd,jnd,ind,ISpaceGrp
    INTEGER(IKIND),DIMENSION(:),ALLOCATABLE :: IDegreesOfFreedom
    REAL(RKIND),DIMENSION(ITHREE,ITHREE) :: RMoveMatrix
    
    CALL ConvertSpaceGroupToNumber(ISpaceGrp,IErr)
    IF(l_alert(IErr,"SetupAtomMovements","ConvertSpaceGroupToNumber")) RETURN  

    ALLOCATE(IDegreesOfFreedom(SIZE(IAtomsToRefine)),STAT=IErr)
    IF(l_alert(IErr,"SetupAtomMovements","allocate IDegreesOfFreedom")) RETURN  

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
    ! called twice in felixrefine, utility for MaxGradient/ParabolicRefinement
    
    REAL(RKIND) :: Ra,Rb,Rc,Rd,Rxv,Ryv
    REAL(RKIND),DIMENSION(3) :: Rx,Ry
    INTEGER(IKIND) :: IErr
    
    Rd = Rx(1)*Rx(1)*(Rx(2)-Rx(3)) + Rx(2)*Rx(2)*(Rx(3)-Rx(1)) + Rx(3)*Rx(3)*(Rx(1)-Rx(2))
    Ra =(Rx(1)*(Ry(3)-Ry(2)) + Rx(2)*(Ry(1)-Ry(3)) + Rx(3)*(Ry(2)-Ry(1)))/Rd
    Rb =( Rx(1)*Rx(1)*(Ry(2)-Ry(3)) + Rx(2)*Rx(2)*(Ry(3)-Ry(1)) + &
          Rx(3)*Rx(3)*(Ry(1)-Ry(2)) )/Rd
    Rc =(Rx(1)*Rx(1)*(Rx(2)*Ry(3)-Rx(3)*Ry(2)) + Rx(2)*Rx(2)*(Rx(3)*Ry(1)-Rx(1)*Ry(3))&
        +Rx(3)*Rx(3)*(Rx(1)*Ry(2)-Rx(2)*Ry(1)))/Rd
    Rxv = -Rb/(2*Ra);!x-coord
    Ryv = Rc-Rb*Rb/(4*Ra)!y-coord

  END SUBROUTINE Parabo3 

END PROGRAM Felixrefine
