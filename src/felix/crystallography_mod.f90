!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: :VERSION: RB_coord / 1.14 /
! Date:    :DATE: 15-01-2019
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

!>
!! Module-description: This defines lattice vectors as well as the fractional atomic coordinates
!!
MODULE crystallography_mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ReciprocalLattice, UniqueAtomPositions, gVectors

  CONTAINS

  !>
  !! Procedure-description: Calculates g-vector matrices, global variables
  !!
  !! Author: Richard Beanland (2019)
  !!
  SUBROUTINE gVectors(IErr)

    USE MyNumbers
    USE message_mod
    USE MyMPI 
 
    ! global inputs
    USE RPARA, ONLY : RarVecM,RbrVecM,RcrVecM,RNormDirM,Rhkl,RConvergenceAngle
    USE IPARA, ONLY : INhkl,IPixelCount
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
!    !resolution in k-space
!    RDeltaK = RConvergenceAngle/REAL(IPixelCount,RKIND)
    
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
  !! Procedure-description: Creates reciprocal lattice vectors in Microscope
  !! reference frame. This involves transforms between the orthogonal, crystal
  !! and microscopic frame. 
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE ReciprocalLattice(IErr)

    USE MyNumbers
    USE MyMPI

    USE RPARA, ONLY : RarVecM,RbrVecM,RcrVecM,RaVecM,RbVecM,RcVecM,RNormDirM,RaVecO,RbVecO,&
          RcVecO,RVolume,RarVecO,RbrVecO,RcrVecO
    USE SPARA, ONLY : SSpaceGroupName
    
    ! global inputs
    USE IPARA, ONLY : IDiffractionFLAG,IVolumeFLAG
    USE RPARA, ONLY : RAlpha,RBeta,RGamma,RLengthX,RLengthY,RLengthZ,RNormDirC,RXDirC,&
          RZDirC
    USE SPARA, ONLY : SPrintString
    
    IMPLICIT NONE

    INTEGER(IKIND) :: IErr,ind
    REAL(RKIND) :: RTTest
    REAL(RKIND), DIMENSION(ITHREE) :: RXDirO, RYDirO, RZDirO, RYDirC
    REAL(RKIND), DIMENSION(ITHREE,ITHREE) :: RTMatC2O,RTMatO2M
    CHARACTER(50) :: indString
    CHARACTER(400) :: RTMatString

    !direct lattice vectors in an orthogonal reference frame, Angstrom units 
    RaVecO(1)= RLengthX
    RaVecO(2)= ZERO
    RaVecO(3)= ZERO

    RbVecO(1)= RLengthY*COS(RGamma)
    RbVecO(2)= RLengthY*SIN(RGamma)
    RbVecO(3)= ZERO

    RcVecO(1)= RLengthZ*COS(RBeta)
    RcVecO(2)= RLengthZ*(COS(RAlpha)-COS(RBeta)*COS(RGamma))/SIN(RGamma)
    RcVecO(3)= RLengthZ*(SQRT(1.D0-COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)&
      -COS(RGamma)*COS(RGamma)+TWO*COS(RAlpha)*COS(RBeta)*COS(RGamma)) / SIN(RGamma))

    !calculate cell volume if required
    IF(IVolumeFLAG .EQ. 0) THEN
       RVolume= RLengthX*RLengthY*RLengthZ* &
            SQRT(1.0D0 - &
            COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)-COS(RGamma)*COS(RGamma) + &
            TWO*COS(RAlpha)*COS(RBeta)*COS(RGamma))
    END IF

    !Some checks for rhombohedral cells?
    IF(IDiffractionFLAG.EQ.0) THEN  
      RTTest = &
            DOT_PRODUCT(RaVecO/DOT_PRODUCT(RaVecO,RaVecO),RbVecO/DOT_PRODUCT(RbVecO,RbVecO))*&
            DOT_PRODUCT(RbVecO/DOT_PRODUCT(RbVecO,RbVecO),RcVecO/DOT_PRODUCT(RcVecO,RcVecO))*&
            DOT_PRODUCT(RcVecO/DOT_PRODUCT(RcVecO,RcVecO),RaVecO/DOT_PRODUCT(RaVecO,RaVecO))
       
       IF(SCAN(SSpaceGroupName,'rR').NE.0) THEN
          IF(ABS(RTTest).LT.TINY) THEN
             SSpaceGroupName = TRIM(ADJUSTL("V"))
             ! Crystal is either Obverse or Reverse
             ! Selection Rules are not in place to determine the difference, 
             ! assume the crystal is Obverse")
          ELSE
             SSpaceGroupName=TRIM(ADJUSTL('P'))
             ! Primitive setting (Rhombohedral axes)
          END IF
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

    ! RTmat transforms from crystal (implicit units) 
    ! to orthogonal reference frame (Angstrom units)
    RTMatC2O(:,1)= RaVecO(:)
    RTMatC2O(:,2)= RbVecO(:)
    RTMatC2O(:,3)= RcVecO(:)
    
    ! RXDirC is the reciprocal lattice vector that defines the x-axis of the
    ! diffraction pattern and RZDirC the beam direction, coming from felix.inp
    ! No check has been made to ensure that they are perpendicular, it is
    ! assumed
    ! RXDirO,RYDirO,RZDirO vectors are UNIT reciprocal lattice vectors parallel 
    ! to the above in an orthogonal frame
    RXDirO= RXDirC(1)*RarVecO + RXDirC(2)*RbrVecO + RXDirC(3)*RcrVecO
    RXDirO= RXDirO/SQRT(DOT_PRODUCT(RXDirO,RXDirO))
    RZDirO= RZDirC(1)*RaVecO + RZDirC(2)*RbVecO + RZDirC(3)*RcVecO
    RZDirO= RZDirO/SQRT(DOT_PRODUCT(RZDirO,RZDirO))
    RYDirO= CROSS(RZDirO,RXDirO)

    ! RTmatO2M transforms from orthogonal to microscope reference frame
    RTMatO2M(1,:)= RXDirO(:)
    RTMatO2M(2,:)= RYDirO(:)
    RTMatO2M(3,:)= RZDirO(:)
    
    ! Unit normal to the specimen in REAL space
    ! This is used in diffraction pattern calculation
    RNormDirM = MATMUL(RTMatO2M,MATMUL(RTMatC2O,RNormDirC))
    RNormDirM = RNormDirM/SQRT(DOT_PRODUCT(RNormDirM,RNormDirM)) 

    ! now transform from crystal reference frame to orthogonal and then to microscope frame

    ! RaVecM, RbVecM, RbVecM unit cell vectors in Angstrom units in the microscope frame
    RaVecM= MATMUL(RTMatO2M,RaVecO)
    RbVecM= MATMUL(RTMatO2M,RbVecO)
    RcVecM= MATMUL(RTMatO2M,RcVecO)
    
    ! create new set of reciprocal lattice vectors in Microscope reference frame
    ! Note that reciprocal lattice vectors have two pi included,
    ! we are using the optical convention exp(i*g.r)
    RarVecM= TWOPI*CROSS(RbVecM,RcVecM)/DOT_PRODUCT(RbVecM,CROSS(RcVecM,RaVecM))
    RbrVecM= TWOPI*CROSS(RcVecM,RaVecM)/DOT_PRODUCT(RcVecM,CROSS(RaVecM,RbVecM))
    RcrVecM= TWOPI*CROSS(RaVecM,RbVecM)/DOT_PRODUCT(RaVecM,CROSS(RbVecM,RcVecM))
!debugging output, can be deleted
!WRITE(SPrintString,FMT='(A4,3(F7.4,1X))') "a*: ",RarVecM
!IF(my_rank.EQ.0)PRINT*,TRIM(ADJUSTL(SPrintString))
!WRITE(SPrintString,FMT='(A4,3(F7.4,1X))') "b*: ",RbrVecM
!IF(my_rank.EQ.0)PRINT*,TRIM(ADJUSTL(SPrintString))
!WRITE(SPrintString,FMT='(A4,3(F7.4,1X))') "c*: ",RcrVecM
!IF(my_rank.EQ.0)PRINT*,TRIM(ADJUSTL(SPrintString))
    
  END SUBROUTINE ReciprocalLattice

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Calculates the full set of possible fractional atomic positions
  !! and then gets rid of duplicates
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016, 2017)
  !!
  SUBROUTINE UniqueAtomPositions(IErr)

    ! updates full crystal arrays from basis atom refinement
    
    USE MyNumbers
    USE message_mod

    ! global outputs
    USE RPARA, ONLY : RAtomCoordinate,ROccupancy,RIsoDW,RAtomPosition
    USE IPARA, ONLY : IAtomicNumber,IAnisoDW
    USE SPARA, ONLY : SAtomLabel, SAtomName

    ! global inputs
    USE RPARA, ONLY : RBasisOccupancy,RBasisIsoDW,RSymVec,RBasisAtomPosition,RSymMat, &
          RcVecM,RbVecM,RaVecM
    USE SPARA, ONLY : SBasisAtomLabel, SBasisAtomName
    USE IPARA, ONLY : IBasisAtomicNumber, IBasisAnisoDW, IMaxPossibleNAtomsUnitCell, &
         INAtomsUnitCell
    USE SPARA, ONLY : SPrintString
    
    IMPLICIT NONE
    
    INTEGER(IKIND) :: IErr,ind,jnd,knd
    INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: IAllAtomicNumber, RAllAnisoDW
    REAL(RKIND),ALLOCATABLE :: RAllAtomPosition(:,:), RAllOccupancy(:), RAllIsoDW(:)
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
!DBG      IF (my_rank.EQ.0) PRINT*, SAtomName(ind)//SPrintString
      CALL message( LL, dbg7, SAtomName(ind)//SPrintString )
      CALL message( LL, dbg7, "(DWF, occupancy) = ",(/ RIsoDW(ind), ROccupancy(ind) /) )
    END DO
    
    ! Finished with these variables now
    DEALLOCATE( &
         RAllAtomPosition, SAllAtomName, RAllOccupancy, RAllIsoDW, &
         IAllAtomicNumber, RAllAnisoDW, STAT=IErr)
    IF(l_alert(IErr,"UniqueAtomPositions","deallocations")) RETURN
      
    !--------------------------------------------------------------------
    ! Calculate atomic position vectors RAtomCoordinate
    !--------------------------------------------------------------------

    ! calculate from Fractional Coordinates and Lattice Vectors
    ! In microscope reference frame, in Angstrom units
    DO ind=1,INAtomsUnitCell
      DO jnd=1,ITHREE
        RAtomCoordinate(ind,jnd)= RAtomPosition(ind,1)*RaVecM(jnd) + &
              RAtomPosition(ind,2)*RbVecM(jnd)+RAtomPosition(ind,3)*RcVecM(jnd)
      END DO
    END DO

  END SUBROUTINE UniqueAtomPositions

END MODULE crystallography_mod
