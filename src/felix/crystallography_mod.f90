!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: 2.0
! Date: 19-12-2022
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
!! Module-description: This defines lattice vectors as well as the fractional atomic coordinates
!!
MODULE crystallography_mod
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: ReciprocalLattice, CrystalOrientation, UniqueAtomPositions, gVectors

  CONTAINS

  !>
  !! Procedure-description: Calculates g-vector matrices, global variables
  !!
  !! Author:  r.beanland@warwick.ac.uk
  !!
  SUBROUTINE gVectors(IErr)

    USE MyNumbers
    USE message_mod
    USE MyMPI 
 
    ! global inputs
    USE RPARA, ONLY : RarVecM,RbrVecM,RcrVecM,RNormDirM,Rhkl,RConvergenceAngle
    USE IPARA, ONLY : INhkl
    USE SPARA, ONLY : SPrintString    
    ! global outputs
    USE RPARA, ONLY : RgPool,RgPoolMag,RgDotNorm,RDeltaK,RgMatrix
    
    IMPLICIT NONE
    
    INTEGER(IKIND) :: IErr,ind,jnd

    IErr=0!No route to throw an error here in fact
    
    !calculate g-vector pool, the magnitudes and component parallel to specimen surface
    DO ind=1,INhkl
      DO jnd=1,ITHREE
        RgPool(ind,jnd) = Rhkl(ind,1)*RarVecM(jnd) + &
            Rhkl(ind,2)*RbrVecM(jnd) + Rhkl(ind,3)*RcrVecM(jnd)
      END DO
      RgPoolMag(ind)= SQRT(DOT_PRODUCT(RgPool(ind,:),RgPool(ind,:)))
      RgDotNorm(ind) = DOT_PRODUCT(RgPool(ind,:),RNormDirM)
    END DO
    
    ! Calculate matrix  of g-vectors that corresponds to the Ug matrix
    DO ind=1,INhkl
      DO jnd=1,INhkl
        RgMatrix(ind,jnd,:)= RgPool(ind,:)-RgPool(jnd,:)
      END DO
    END DO

    !outputs if requested    
    CALL message(LL,dbg3,"first 16 g-vectors", RgMatrix(1:16,1,:)) 
    CALL message(LL,dbg7,"g-vectors and magnitude (1/A), in the microscope reference frame" )
    DO ind = 1,INhkl
      CALL message(LL,dbg7,"hkl  :",NINT(Rhkl(ind,:)))
      CALL message(LL,dbg7,"g mag:",RgPoolMag(ind))
    END DO
    CALL message(LL,dbg7,"g.n list")
    DO ind = 1,INhkl
      CALL message(LL,dbg7,"hkl :",NINT(Rhkl(ind,:)))
      CALL message(LL,dbg7,"g.n :",RgDotNorm(ind))
    END DO
    
    END SUBROUTINE gVectors

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Creates list of reciprocal lattice vectors covering a 3D
  !! volume in an orthogonal reference frame and sorts them in magnitude.  
  !!
  !! Major-Authors: Richard Beanland (2023)
  !!
  SUBROUTINE ReciprocalLattice(RLatticeLimit, IErr)

    USE MyNumbers
    USE MyMPI
    USE message_mod

    ! global inputs
    USE IPARA, ONLY : IVolumeFLAG
    USE RPARA, ONLY : RAlpha,RBeta,RGamma,RCellA,RCellB,RCellC,RXDirC_0,RZDirC_0
    USE SPARA, ONLY : SSpaceGroupName,SPrintString
    USE CPARA, ONLY : CFg

    ! global outputs
    USE RPARA, ONLY : RaVecO,RbVecO,RcVecO,RVolume,RarVecO,RbrVecO,RcrVecO,&
            RXDirO,RYDirO,RZDirO,RarMag,RbrMag,RcrMag,RgLatticeO,RLatMag
    USE IPARA, ONLY : IhklLattice,InLattice

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr,ind,jnd,knd,lnd,mnd,nnd,inda,indb,indc,IlogNB2,Id3(ITHREE),ISel
    REAL(RKIND) :: Rt,ALN2I,LocalTINY,Rg(ITHREE),RxAngle
    REAL(RKIND),INTENT(IN) :: RLatticeLimit
    PARAMETER (ALN2I=1.4426950D0, LocalTINY=1.D-5)
    
    !direct lattice vectors in an orthogonal reference frame, Angstrom units 
    ! a is parallel to [100]
    RaVecO(1)= RCellA
    RaVecO(2)= ZERO
    RaVecO(3)= ZERO
    ! b lies in the x-y plane 
    RbVecO(1)= RCellB*COS(RGamma)
    RbVecO(2)= RCellB*SIN(RGamma)
    RbVecO(3)= ZERO
    ! c lies... wherever
    RcVecO(1)= RCellC*COS(RBeta)
    RcVecO(2)= RCellC*(COS(RAlpha)-COS(RBeta)*COS(RGamma))/SIN(RGamma)
    RcVecO(3)= RCellC*(SQRT(1.D0-COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)&
      -COS(RGamma)*COS(RGamma)+TWO*COS(RAlpha)*COS(RBeta)*COS(RGamma)) / SIN(RGamma))

    !calculate cell volume if required
    IF(IVolumeFLAG .EQ. 0) THEN
       RVolume= RCellA*RCellB*RCellC* &
            SQRT(1.0D0 - &
            COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)-COS(RGamma)*COS(RGamma) + &
            TWO*COS(RAlpha)*COS(RBeta)*COS(RGamma))
    END IF

    !Some checks for rhombohedral cells?
    Rt = DOT_PRODUCT(RaVecO/DOT_PRODUCT(RaVecO,RaVecO),RbVecO/DOT_PRODUCT(RbVecO,RbVecO))*&
         DOT_PRODUCT(RbVecO/DOT_PRODUCT(RbVecO,RbVecO),RcVecO/DOT_PRODUCT(RcVecO,RcVecO))*&
         DOT_PRODUCT(RcVecO/DOT_PRODUCT(RcVecO,RcVecO),RaVecO/DOT_PRODUCT(RaVecO,RaVecO))
       
    IF(SCAN(SSpaceGroupName,'rR').NE.0) THEN
       IF(ABS(Rt).LT.TINY) THEN
        SSpaceGroupName = TRIM(ADJUSTL("V"))
        ! Crystal is either Obverse or Reverse
        ! Selection Rules are not in place to determine the difference, 
        ! assume the crystal is Obverse")
      ELSE
        SSpaceGroupName=TRIM(ADJUSTL('P'))
        ! Primitive setting (Rhombohedral axes)
      END IF
    END IF

    ! Set up Reciprocal Lattice Vectors: orthogonal reference frame in 1/Angstrom units
    ! RarDirO,RbrDirO,RcrDirO vectors are reciprocal lattice vectors 
    ! 2pi/a, 2pi/b, 2pi/c in an orthogonal frame
    ! Note that reciprocal lattice vectors have two pi included,
    ! we are using the physics convention exp(i*g.r)
    RarVecO= TWOPI*CROSS(RbVecO,RcVecO)/DOT_PRODUCT(RbVecO,CROSS(RcVecO,RaVecO))
    RbrVecO= TWOPI*CROSS(RcVecO,RaVecO)/DOT_PRODUCT(RcVecO,CROSS(RaVecO,RbVecO))
    RcrVecO= TWOPI*CROSS(RaVecO,RbVecO)/DOT_PRODUCT(RaVecO,CROSS(RbVecO,RcVecO))
    DO ind=1,ITHREE
       IF (abs(RarVecO(ind)).LT.TINY) THEN
          RarVecO(ind) = ZERO
       END IF
       IF (abs(RbrVecO(ind)).LT.TINY) THEN
          RbrVecO(ind) = ZERO
       END IF
       IF (abs(RcrVecO(ind)).LT.TINY) THEN
          RcrVecO(ind) = ZERO
       END IF
    ENDDO
    ! their magnitudes
    RarMag=SQRT(DOT_PRODUCT(RarVecO,RarVecO))!magnitude of a*
    RbrMag=SQRT(DOT_PRODUCT(RbrVecO,RbrVecO))!magnitude of b*
    RcrMag=SQRT(DOT_PRODUCT(RcrVecO,RcrVecO))!magnitude of c*

    ! Now build the reciprocal lattice, from which we will take slices for
    ! each beam pool using HKLMake
    ! IhklLattice is the list of Miller indices for the 3D lattice
    ! RgLatticeO is the corresponding list of coordinates in reciprocal space
    ! maximum a*,b*,c* limit is determined by the G magnitude limit
    inda=NINT(RLatticeLimit/RarMag)
    indb=NINT(RLatticeLimit/RbrMag)
    indc=NINT(RLatticeLimit/RcrMag)
    InLattice = (2*inda+1)*(2*indb+1)*(2*indc+1)
    ALLOCATE(IhklLattice(InLattice, ITHREE), STAT=IErr)! Miller indices
    ALLOCATE(RgLatticeO(InLattice, ITHREE), STAT=IErr)! g-vector
    ALLOCATE(RLatMag(InLattice), STAT=IErr)! magnitude
    ALLOCATE(CFg(InLattice), STAT=IErr)! Structure factors
    
    ! populate the lists
    lnd = 0
    DO ind = -inda,inda
      DO jnd = -indb,indb
        DO knd = -indc,indc
          lnd = lnd + 1
          CALL SelectionRules(ind, jnd, knd, ISel, IErr)! Systematic absences - lattice
          IF (ISel.EQ.1) THEN
            IhklLattice(lnd,:) = (/ ind, jnd, knd /) !Miller indices
            Rg = ind*RarVecO + jnd*RbrVecO + knd*RcrVecO !g-vector
            RgLatticeO(lnd,:) = Rg  ! in the O reference frame
            RLatMag(lnd) = SQRT(DOT_PRODUCT(Rg,Rg)) !g-magnitude
            ! Calculate structure factor
            DO knd=1,INAtomsUnitCell
              CALL AtomicScatteringFactor(RScatteringFactor,IErr)
              CFg(lnd) = RScatteringFactor*EXP(-CIMAGONE*DOT_PRODUCT(Rg,RAtomCoordinate(knd,:)) )
            END DO
          END IF
          END DO
      END DO
    END DO
    
    ! Sort them in ascending order of magnitude (re-purposed HKLSort routine)
    ! Based on ShellSort from "Numerical Recipes", routine SHELL()
    IlogNB2=INT(LOG(REAL(InLattice))*ALN2I+LocalTINY)
    mnd = InLattice
    DO nnd=1,IlogNB2
      mnd=mnd/2
      knd=InLattice-mnd
      DO jnd=1,knd
        ind=jnd
3       CONTINUE
        lnd=ind+mnd
        IF( RLatMag(lnd) .LT. RLatMag(ind)) THEN
          Id3 = IhklLattice(ind,:) ! swap indices
          IhklLattice(ind,:) = IhklLattice(lnd,:)
          IhklLattice(lnd,:) = Id3
          Rg = RgLatticeO(ind,:) ! swap g-vectors
          RgLatticeO(ind,:) = RgLatticeO(lnd,:)
          RgLatticeO(lnd,:) = Rg
          Rt = RLatMag(ind) ! swap magnitudes
          RLatMag(ind) = RLatMag(lnd)
          RLatMag(lnd) = Rt
          ind=ind-mnd
          IF(ind.GE.1) GOTO 3
        ENDIF
      ENDDO
    ENDDO

    ! Make the final g-vector something big, to fill up incomplete matrices later
!    IhklLattice(InLattice,:) = (/ 666,666,666 /)
!    RgLatticeO(InLattice,:) = REAL( (/ 666,666,666 /),RKIND )
!    RLatMag(InLattice) = 666666.666

    ! Set up initial microscope reference frame
    ! X, Y and Z are orthogonal vectors that defines the simulation
    ! Also referred to as the microscope reference frame M.
    ! The electron beam propagates along +Zm.
    ! The alpha rotation axis is along Ym.  Positive alpha rotation moves the field
    ! of view of the simulation along +Xm.
    ! In the crystal reference frame we read in reciprocal vectors RXDirC_0 and RZDirC_0
    ! These define Xm & Zm in the inital reference frame
    ! RXDirO,RYDirO,RZDirO are UNIT reciprocal lattice vectors parallel to X,Y,Z
    RXDirO = RXDirC_0(1)*RarVecO + RXDirC_0(2)*RbrVecO + RXDirC_0(3)*RcrVecO
    RXDirO = RXDirO/SQRT(DOT_PRODUCT(RXDirO,RXDirO))
    RZDirO = RZDirC_0(1)*RaVecO + RZDirC_0(2)*RbVecO + RZDirC_0(3)*RcVecO
    RZDirO = RZDirO/SQRT(DOT_PRODUCT(RZDirO,RZDirO))
    ! Check the input is sensible, i.e. Xm is perpendicular to Zm
    RxAngle = ABS(180.0D0*ACOS(DOT_PRODUCT(RXDirO,RZDirO))/PI)
    IF(ABS(RxAngle-90.0D0).GT.0.1)THEN! with a tolerance of 0.1 degrees
      WRITE(SPrintString,"(A15,F5.1,A27)") "Error: X is at ",RxAngle," degrees to Z, should be 90"
      CALL message(LS,SPrintString)
      IErr = 1
    ELSE!fine correction of x
      ! take off any component parallel to z & renormalise
      RXDirO = RXDirO - DOT_PRODUCT(RXDirO,RZDirO)*RZDirO
      RXDirO = RXDirO/SQRT(DOT_PRODUCT(RXDirO,RXDirO))
    END IF
    RYDirO = CROSS(RZDirO,RXDirO)  ! the rotation axis

  END SUBROUTINE ReciprocalLattice

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE CrystalOrientation(IErr)
    USE MyNumbers
    USE MyMPI
    USE message_mod

    USE RPARA, ONLY : RarVecM,RbrVecM,RcrVecM,RaVecM,RbVecM,RcVecM,RNormDirM,RaVecO,RbVecO,&
          RcVecO,RarVecO,RbrVecO,RcrVecO,RXDirO,RYDirO,RZDirO,RarMag,RbrMag,RcrMag
    USE SPARA, ONLY : SSpaceGroupName

    ! global inputs
    USE IPARA, ONLY : IVolumeFLAG,INAtomsUnitCell
    USE RPARA, ONLY : RAlpha,RBeta,RGamma,RCellA,RCellB,RCellC,RNormDirC,RXDirC,&
          RZDirC,RXDirM,RAtomCoordinate,RAtomPosition
    USE SPARA, ONLY : SPrintString

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr,ind,jnd
    REAL(RKIND), DIMENSION(ITHREE,ITHREE) :: RTMatC2O,RTMatO2M
    REAL(RKIND) :: RNormAngle

    ! RTmatC2O transforms from crystal (implicit units) 
    ! to orthogonal reference frame (Angstrom units)
    RTMatC2O(:,1) = RaVecO(:)
    RTMatC2O(:,2) = RbVecO(:)
    RTMatC2O(:,3) = RcVecO(:)

    ! RTmatO2M transforms from orthogonal to microscope reference frame
    RTMatO2M(1,:) = RXDirO(:)
    RTMatO2M(2,:) = RYDirO(:)
    RTMatO2M(3,:) = RZDirO(:)

    ! Unit normal to the specimen in REAL space
    ! This is used for all g-vectors as a boundary condition
    RNormDirM = MATMUL(RTMatO2M,MATMUL(RTMatC2O,RNormDirC))
    RNormDirM = RNormDirM/SQRT(DOT_PRODUCT(RNormDirM,RNormDirM)) 

    ! Check to see if the normal is close to perpendicular to the incident beam
    ! remember in the microscope frame z = [001] so n.z is just n(3)
    RNormAngle = ABS(180.0D0*ACOS(RNormDirM(3))/PI)
    IF(RNormAngle.GT.75.0D0 .AND. my_rank.EQ.0) THEN
      WRITE(SPrintString,"(A30,F5.1,A8)") "Warning: surface normal is at ",RNormAngle," degrees"
      CALL message(LS,SPrintString) 
    END IF

    ! Transform atomic coordinates to microscope frame to give RAtomCoordinate
    ! RaVecM, RbVecM, RbVecM unit cell vectors in Angstrom units in the microscope frame
    RaVecM= MATMUL(RTMatO2M,RaVecO)
    RbVecM= MATMUL(RTMatO2M,RbVecO)
    RcVecM= MATMUL(RTMatO2M,RcVecO)
    ! Calculate atomic position vectors RAtomCoordinate
    ! In microscope reference frame, in Angstrom units (NB RAtomPosition=crystal frame, in .cif)
    DO ind=1,INAtomsUnitCell
      DO jnd=1,ITHREE
        RAtomCoordinate(ind,jnd)= RAtomPosition(ind,1)*RaVecM(jnd) + &
              RAtomPosition(ind,2)*RbVecM(jnd)+RAtomPosition(ind,3)*RcVecM(jnd)
      END DO
    END DO
    
    ! create reciprocal lattice vectors in Microscope reference frame
    ! Note that reciprocal lattice vectors have two pi included,
    ! we are using the optical convention exp(i*g.r)
    RarVecM= TWOPI*CROSS(RbVecM,RcVecM)/DOT_PRODUCT(RbVecM,CROSS(RcVecM,RaVecM))  ! a*
    RbrVecM= TWOPI*CROSS(RcVecM,RaVecM)/DOT_PRODUCT(RcVecM,CROSS(RaVecM,RbVecM))  ! b*
    RcrVecM= TWOPI*CROSS(RaVecM,RbVecM)/DOT_PRODUCT(RaVecM,CROSS(RbVecM,RcVecM))  ! c*
    ! their magnitudes
    RarMag=SQRT(DOT_PRODUCT(RarVecM,RarVecM))!magnitude of a*
    RbrMag=SQRT(DOT_PRODUCT(RbrVecM,RbrVecM))!magnitude of b*
    RcrMag=SQRT(DOT_PRODUCT(RcrVecM,RcrVecM))!magnitude of c*

  END SUBROUTINE CrystalOrientation

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Calculates the full set of possible fractional atomic positions,
  !! the mean inner potential and wavevector in the material K
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016, 2017)
  !!
  SUBROUTINE UniqueAtomPositions(IErr)

    ! updates full crystal arrays from basis atom refinement
    
    USE MyNumbers
    USE message_mod
    USE ug_matrix_mod

    ! global outputs
    USE RPARA, ONLY : RAtomCoordinate,ROccupancy,RIsoDW,RAtomPosition,RMeanInnerPotential,RBigK
    USE IPARA, ONLY : IAtomicNumber,IAnisoDW
    USE SPARA, ONLY : SAtomLabel, SAtomName

    ! global inputs
    USE RPARA, ONLY : RBasisOccupancy,RBasisIsoDW,RSymVec,RBasisAtomPosition,RSymMat, &
          RcVecM,RbVecM,RaVecM,RCurrentGMagnitude,RScattFacToVolts,RElectronWaveVectorMagnitude
    USE SPARA, ONLY : SBasisAtomLabel, SBasisAtomName
    USE IPARA, ONLY : IBasisAtomicNumber, IBasisAnisoDW, IMaxPossibleNAtomsUnitCell, &
         INAtomsUnitCell,ICurrentZ
    USE SPARA, ONLY : SPrintString
    
    IMPLICIT NONE
    
    INTEGER(IKIND) :: IErr,ind,jnd,knd
    INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: IAllAtomicNumber, RAllAnisoDW
    REAL(RKIND),ALLOCATABLE :: RAllAtomPosition(:,:), RAllOccupancy(:), RAllIsoDW(:)
    REAL(RKIND) :: RScatteringFactor
    LOGICAL :: Lunique
    CHARACTER(2), DIMENSION(:), ALLOCATABLE :: SAllAtomName
    CHARACTER(5), DIMENSION(:), ALLOCATABLE :: SAllAtomLabel

    ind=SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)
    ! All atom positions generated by symmetry, including duplicates: local variables
    ALLOCATE( RAllAtomPosition(ind,ITHREE),SAllAtomName(ind),SAllAtomLabel(ind),&
        RAllOccupancy(ind),RAllIsoDW(ind),IAllAtomicNumber(ind),RAllAnisoDW(ind),STAT=IErr )
    IF(l_alert(IErr,"UniqueAtomPositions","allocations")) RETURN
   
    !--------------------------------------------------------------------  
    ! apply symmetry elements to generate all equivalent positions 
    !--------------------------------------------------------------------  
    knd=1
    DO ind=1, SIZE(RSymVec,1)  
      DO jnd=1, SIZE(RBasisAtomPosition,1)
        RAllAtomPosition(knd,:)= MATMUL(RSymMat(ind,:,:),RBasisAtomPosition(jnd,:)) &
              + RSymVec(ind,:)
        SAllAtomLabel(knd) = SBasisAtomLabel(jnd)
        SAllAtomName(knd) = SBasisAtomName(jnd)
        RAllOccupancy(knd) = RBasisOccupancy(jnd)
        RAllIsoDW(knd) = RBasisIsoDW(jnd)
        IAllAtomicNumber(knd) = IBasisAtomicNumber(jnd)
        RAllAnisoDW(knd) = IBasisAnisoDW(jnd)
        knd=knd+1
      END DO
!WRITE(SPrintString,'(F3.0,1X,F3.0,1X,F3.0,2X,F3.0,1X,F3.0,1X,F3.0,2X,F3.0,1X,F3.0,1X,F3.0)') RSymMat(ind,1,:),RSymMat(ind,2,:),RSymMat(ind,3,:)
!IF(my_rank.EQ.0) PRINT*, ind,"RSymMat:  ", SPrintString             
    END DO
    RAllAtomPosition = MODULO(RAllAtomPosition,ONE)
    WHERE(ABS(RAllAtomPosition).LT.TINY) RAllAtomPosition = ZERO 

    !--------------------------------------------------------------------  
    ! Reduce to the set of unique fractional atomic positions
    !--------------------------------------------------------------------
    
    ! first atom has to be in this set
    RAtomPosition(1,:)= RAllAtomPosition(1,:)
    SAtomLabel(1)= SAllAtomLabel(1)
    SAtomName(1)= SAllAtomName(1)
    RIsoDW(1) = RAllIsoDW(1)
    ROccupancy(1) = RAllOccupancy(1)
    IAtomicNumber(1) = IAllAtomicNumber(1)
    IAnisoDW(1) = RAllAnisoDW(1)
    jnd=2
    ! work through all possible atom coords and check for duplicates
    DO ind=2,IMaxPossibleNAtomsUnitCell
      Lunique=.TRUE.
      DO knd=1,jnd-1 ! check against the unique ones found so far
        IF (SUM(ABS(RAllAtomPosition(ind,:)-RAtomPosition(knd,:))).LE.TINY) THEN ! position same
          IF (SAllAtomLabel(ind).EQ.SAtomLabel(knd)) THEN ! Label is the same too, so not unique
            Lunique=.FALSE.
            EXIT
          END IF
        END IF
      END DO
      IF (Lunique .EQV. .TRUE.) THEN
        RAtomPosition(jnd,:)= RAllAtomPosition(ind,:)
        SAtomLabel(jnd)= SAllAtomLabel(ind)
        SAtomName(jnd)= SAllAtomName(ind)
        RIsoDW(jnd) = RAllIsoDW(ind)
        ROccupancy(jnd) = RAllOccupancy(ind)
        IAtomicNumber(jnd) = IAllAtomicNumber(ind)!
        IAnisoDW(jnd) = RAllAnisoDW(ind)
        jnd=jnd+1
      END IF
    END DO
    INAtomsUnitCell = jnd-1 ! this is how many unique atoms there are in the unit cell

    DO ind=1,INAtomsUnitCell    
      CALL message( LL, dbg7, "Atom ",ind)
      WRITE(SPrintString,"(A18,F8.4,F8.4,F8.4)") ": Atom position = ", RAtomPosition(ind,:)
      CALL message( LL, dbg7, SAtomName(ind)//SPrintString )
      CALL message( LL, dbg7, "(DWF, occupancy) = ",(/ RIsoDW(ind), ROccupancy(ind) /) )
    END DO

    !--------------------------------------------------------------------
    ! calculate mean inner potential and wave vector magnitude
    !--------------------------------------------------------------------
    ! calculate the mean inner potential as the sum of scattering factors
    ! at g=0 multiplied by h^2/(2pi*m0*e*CellVolume)
    RMeanInnerPotential=ZERO
    RCurrentGMagnitude=ZERO  ! this is a global variable, sets g=0
    DO ind=1,INAtomsUnitCell
      ICurrentZ = IAtomicNumber(ind)
      CALL AtomicScatteringFactor(RScatteringFactor,IErr)
      CALL message( LL, dbg3, "Atom ",ind)
      CALL message( LL, dbg3, "f(theta) at g=0 ",RScatteringFactor)
      RMeanInnerPotential = RMeanInnerPotential+RScatteringFactor
    END DO
    RMeanInnerPotential = RMeanInnerPotential*RScattFacToVolts
    WRITE(SPrintString,FMT='(A21,F6.2,A6)') "Mean inner potential ",RMeanInnerPotential," Volts"
    SPrintString=TRIM(ADJUSTL(SPrintString))
    CALL message(LS,SPrintString)

    ! Wave vector magnitude in crystal
    ! high-energy approximation (not HOLZ compatible)
    ! K^2=k^2+U0
    RBigK= SQRT(RElectronWaveVectorMagnitude**2)!-RMeanInnerPotential)
    CALL message ( LM, dbg3, "K (Angstroms) = ",RBigK )

    ! Finished with these variables now
    DEALLOCATE( &
         RAllAtomPosition, SAllAtomName, RAllOccupancy, RAllIsoDW, &
         IAllAtomicNumber, RAllAnisoDW, STAT=IErr)
    IF(l_alert(IErr,"UniqueAtomPositions","deallocations")) RETURN

  END SUBROUTINE UniqueAtomPositions

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Checks a g-vector Ih,Ik,Il
  !! against the selection rules for the global variable SSpaceGroupName
  !! IFlag comes in as zero and goes out as 1 if it is an allowed reflection
  !!
  !! Major-Authors: Richard Beanland (2021)
  !!  
  SUBROUTINE SelectionRules(Ih, Ik, Il, ISel, IErr)

    ! This procedure is called from ReciprocalLattice
    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE SPARA, ONLY : SSpaceGroupName
    
    IMPLICIT NONE

    INTEGER (IKIND),INTENT(IN) :: Ih, Ik, Il
    INTEGER (IKIND),INTENT(INOUT) :: ISel, IErr

    SELECT CASE(SSpaceGroupName)

      CASE("F") !Face Centred, all odd or all even
        IF( (MOD(Ih+Ik,2).EQ.0).AND.&
          (MOD(Ik+Il,2).EQ.0).AND.&
          (MOD(Il+Ih,2).EQ.0) ) ISel=1

      CASE("I")! Body Centred
        IF(MOD(Ih+Ik+Il,2).EQ.0) ISel=1

      CASE("A")! A-Face Centred
        IF(MOD(Ik+Il,2).EQ.0) ISel=1

      CASE("B")! B-Face Centred
        IF(MOD(Ih+Il,2).EQ.0) ISel=1

      CASE("C")! C-Face Centred
        IF(MOD(Ih+Ik,2).EQ.0) ISel=1

      CASE("R")! Rhombohedral Reverse
        IF(MOD(Ih-Ik+Il,3).EQ.0) ISel=1

      CASE("V")! Rhombohedral Obverse
        IF(MOD(-Ih+Ik+Il,3).EQ.0) ISel=1

      CASE("P")! Primitive
        ISel=1

      CASE DEFAULT
      IErr=1
      IF(l_alert(IErr,"SelectionRules",&
          "Space Group Name unrecognised")) RETURN
          
    END SELECT
     
  END SUBROUTINE SelectionRules
  
END MODULE crystallography_mod
