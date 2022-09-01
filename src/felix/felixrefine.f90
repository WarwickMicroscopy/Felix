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
  USE refinementcontrol_mod       ! CONTAINS Simulate
  USE write_output_mod

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
        IAnisotropicDebyeWallerFactorElementNo,IStartTime,IStartTime2,&
        IXPixelIndex,IYPixelIndex
  INTEGER(4) :: IErr4
  REAL(RKIND) :: REmphasis,RGlimit,RLaueZoneGz,RMaxGMag,&
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
  
  CALL ReadInpFile(IErr) ! felix.inp
  IF(l_alert(IErr,"felixrefine","ReadInpFile")) CALL abort
  CALL SetMessageMode( IWriteFLAG, IErr )
  IF(l_alert(IErr,"felixrefine","set_message_mod_mode")) CALL abort

  CALL read_cif(IErr) ! felix.cif ! some allocations are here
  IF(l_alert(IErr,"felixrefine","ReadCif")) CALL abort

  CALL ReadHklFile(IErr) ! the list of hkl's to input/output
  IF(l_alert(IErr,"felixrefine","ReadHklFile")) CALL abort

  WRITE(SPrintString, FMT='(A11,I6,1x,A1,I6)') "Simulation ",ISizeX,"x",ISizeY 
  CALL message(LS,SPrintString)
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
  !resolution in k-space N.B. in cRED we define convergence angle as half the y-size
  RDeltaK = FOURPI*RConvergenceAngle/REAL(ISizeY,RKIND)
!DBG   IF(my_rank.EQ.0)PRINT*, "Delta K", RDeltaK  

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
  ALLOCATE(Rhkl(INhkl,ITHREE),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate Rhkl")) CALL abort
  RGlimit = 10.0*TWOPI    
  CALL HKLMake(RGlimit,IErr)
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
  ! allocate, ImageInitialisation, ImageMaskInitialisation
  !--------------------------------------------------------------------
  ALLOCATE(RhklPositions(INhkl,2),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RhklPositions")) CALL abort

  IPixelTotal = ISizeX*ISizeY
  ALLOCATE(IPixelLocations(IPixelTotal,2),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate IPixelLocations")) CALL abort
  ! we keep track of where a calculation goes in the image using two
  ! 1D IPixelLocations arrays.  Remember fortran indexing is [row,col]=[y,x]
  knd=0
  DO IYPixelIndex = 1,ISizeY
    DO IXPixelIndex = 1,ISizeX
      knd = knd+1
      IPixelLocations(knd,1) = IYPixelIndex
      IPixelLocations(knd,2) = IXPixelIndex
    END DO
  END DO


  !--------------------------------------------------------------------
  ! allocate & setup image arrays for pixel-parallel simulations
  !--------------------------------------------------------------------
  ! All the individual calculations go into RSimulatedPatterns later with MPI_GATHERV
  ! NB RSimulatedPatterns is a vector with respect to pixels, not a 2D image
  ALLOCATE(RSimulatedPatterns(INoOfLacbedPatterns,IThicknessCount,IPixelTotal),STAT=IErr)
  IF(l_alert(IErr,"felixrefine","allocate RSimulatedPatterns")) CALL abort
  ! Images, NB Fortan arrays are [row,column]=[y,x]
  ALLOCATE(RImageSimi(ISizeY,ISizeX,INoOfLacbedPatterns,IThicknessCount),&
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
  ! baseline simulation with timer
  CALL Simulate(IErr)
  IF(l_alert(IErr,"felixrefine","Simulate")) CALL abort
  CALL PrintEndTime(LS,IStartTime2, "Simulation" )

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

  SUBROUTINE abort
    IErr=1
    IF(l_alert(IErr,"felixrefine","ABORTING")) CONTINUE
    CALL MPI_Abort(MPI_COMM_WORLD,1,IErr)
    STOP
  END SUBROUTINE abort

END PROGRAM Felixrefine

