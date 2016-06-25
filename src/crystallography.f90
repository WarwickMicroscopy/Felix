!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all rights reserved
!
! Version: $VERSION$
! Date:    $DATE$
! Time:    $TIME$
! Status:  $RLSTATUS$
! Build:   $BUILD$
! Author:  $AUTHOR$
! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This file is part of felixsim.
!
!  felixsim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  felixsim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with felixsim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ReciprocalLattice(IErr)
  
  USE MyNumbers
  USE WriteToScreen
  
  USE RPara; USE IPara; USE SPara
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind
  REAL(RKIND) :: RTTest
  REAL(RKIND), DIMENSION(ITHREE) :: RXDirO, RYDirO, RZDirO, RYDirC
  REAL(RKIND), DIMENSION(ITHREE,ITHREE) :: RTMatC2O,RTMatO2M
  CHARACTER*50 indString
  CHARACTER*400  RTMatString

  CALL Message("ReciprocalLattice",IMust,IErr)

  ! Setup Crystal Lattice Vectors: orthogonal reference frame in Angstrom units

  RaVecO(1)= RLengthX
  RaVecO(2)= ZERO
  RaVecO(3)= ZERO

  RbVecO(1)= RLengthY*COS(RGamma)
  RbVecO(2)= RLengthY*SIN(RGamma)
  RbVecO(3)= ZERO

  RcVecO(1)= RLengthZ*COS(RBeta)
  RcVecO(2)= RLengthZ*(COS(RAlpha)-COS(RBeta)*COS(RGamma))/SIN(RGamma)
  RcVecO(3)= RLengthZ*( &
       SQRT(1.D0- &
        COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)-COS(RGamma)*COS(RGamma) + &
        TWO*COS(RAlpha)*COS(RBeta)*COS(RGamma)) / &
        SIN(RGamma))

  IF(IVolumeFLAG .EQ. 0) THEN
     RVolume= RLengthX*RLengthY*RLengthZ* &
          SQRT(1.0D0 - &
          COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)-COS(RGamma)*COS(RGamma) + &
          TWO*COS(RAlpha)*COS(RBeta)*COS(RGamma))
  END IF

  IF(IDiffractionFLAG.EQ.0) THEN  
     RTTest = DOT_PRODUCT(RaVecO/DOT_PRODUCT(RaVecO,RaVecO),RbVecO/DOT_PRODUCT(RbVecO,RbVecO))*&
          DOT_PRODUCT(RbVecO/DOT_PRODUCT(RbVecO,RbVecO),RcVecO/DOT_PRODUCT(RcVecO,RcVecO))*&
          DOT_PRODUCT(RcVecO/DOT_PRODUCT(RcVecO,RcVecO),RaVecO/DOT_PRODUCT(RaVecO,RaVecO))
     
     IF(SCAN(SSpaceGroupName,'rR').NE.0) THEN
        IF(ABS(RTTest).LT.TINY) THEN
           SSpaceGroupName = TRIM(ADJUSTL("V"))
           CALL Message("ReciprocalLattice",IMust,IErr, &
                MessageString = "Warning: Crystal is either Obverse or Reverse,")
           CALL Message("ReciprocalLattice",IMust,IErr, &
                MessageString = "Selection Rules are not Currently In place to determine the difference,")
           CALL Message("ReciprocalLattice",IMust,IErr, &
                MessageString = "felix will assume the crystal is Obverse")
        ELSE
           SSpaceGroupName=TRIM(ADJUSTL('P'))
           CALL Message("ReciprocalLattice",IMust,IErr, &
                MessageString = "Crystal is in Primitive setting (Rhombohedral axes)")
        END IF
     END IF
  END IF
  ! Set up Reciprocal Lattice Vectors: orthogonal reference frame in 1/Angstrom units
  ! RarDirO,RbrDirO,RcrDirO vectors are reciprocal lattice vectors 2pi/a, 2pi/b, 2pi/c in an orthogonal frame
  ! Note that reciprocal lattice vectors have two pi included, we are using the physics convention exp(i*g.r)
  RarVecO= TWOPI*CROSS(RbVecO,RcVecO)/DOT_PRODUCT(RbVecO,CROSS(RcVecO,RaVecO))
  RbrVecO= TWOPI*CROSS(RcVecO,RaVecO)/DOT_PRODUCT(RcVecO,CROSS(RaVecO,RbVecO))
  RcrVecO= TWOPI*CROSS(RaVecO,RbVecO)/DOT_PRODUCT(RaVecO,CROSS(RbVecO,RcVecO))

  DO ind=1,ITHREE
     IF (abs(RarVecO(ind)).lt.TINY) THEN
        RarVecO(ind) = ZERO
     END IF
     IF (abs(RbrVecO(ind)).lt.TINY) THEN
        RbrVecO(ind) = ZERO
     END IF
     IF (abs(RcrVecO(ind)).lt.TINY) THEN
        RcrVecO(ind) = ZERO
     END IF
  ENDDO

  ! RTmat transforms from crystal (implicit units)to orthogonal reference frame (Angstrom units)
  RTMatC2O(:,1)= RaVecO(:)
  RTMatC2O(:,2)= RbVecO(:)
  RTMatC2O(:,3)= RcVecO(:)

  CALL Message("ReciprocalLattice",IMoreInfo,IErr,MessageString = "RTMatC2O")
  DO ind=1,ITHREE
    WRITE(indString,*)ind
    WRITE(RTMatString,'(3(F8.3,1X))') RTMatC2O(:,ind)
    CALL Message("ReciprocalLattice",IMoreInfo,IErr, &
         MessageVariable = "RTMatC2O(:,"//TRIM(ADJUSTL(indString))//")", &
         MessageString = TRIM(ADJUSTL(RTMatString)))
  END DO
  
  ! RXDirC,RYDirC,RZDirC vectors are the reciprocal lattice vectors that define the x-axis of the
  ! diffraction pattern and the beam direction, they are given in felix.inp
  ! RXDirO,RYDirO,RZDirO vectors are UNIT reciprocal lattice vectors parallel to the above in an orthogonal frame
  RXDirO= RXDirC(1)*RarVecO + RXDirC(2)*RbrVecO + RXDirC(3)*RcrVecO
  RXDirO= RXDirO/SQRT(DOT_PRODUCT(RXDirO,RXDirO))
  RZDirO= RZDirC(1)*RaVecO + RZDirC(2)*RbVecO + RZDirC(3)*RcVecO
  RZDirO= RZDirO/SQRT(DOT_PRODUCT(RZDirO,RZDirO))
  RYDirO= CROSS(RZDirO,RXDirO)

  ! RTmatO2M transforms from orthogonal to microscope reference frame
  RTMatO2M(1,:)= RXDirO(:)
  RTMatO2M(2,:)= RYDirO(:)
  RTMatO2M(3,:)= RZDirO(:)

  CALL Message("ReciprocalLattice",IMoreInfo,IErr,MessageString = "RTMatO2M")
  DO ind=1,ITHREE
     WRITE(indString,*)ind
     WRITE(RTMatString,'(3(F8.3,1X))') RTMatO2M(:,ind)
     CALL Message("ReciprocalLattice",IMoreInfo,IErr, &
          MessageVariable = "RTMatO2M(:,"//TRIM(ADJUSTL(indString))//")", &
          MessageString = TRIM(ADJUSTL(RTMatString)))
  END DO
  
  !Unit normal to the specimen in real space
  !This is used in diffraction pattern calculation subroutine
  RNormDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RNormDirC))
  RNormDirM = RNormDirM/SQRT(DOT_PRODUCT(RNormDirM,RNormDirM)) 
  
  ! now transform from crystal reference frame to orthogonal and then to microscope frame

  ! RaVecM, RbVecM, RbVecM unit cell vectors in Angstrom units in the microscope frame
  RaVecM= MATMUL(RTMatO2M,RaVecO)
  RbVecM= MATMUL(RTMatO2M,RbVecO)
  RcVecM= MATMUL(RTMatO2M,RcVecO)
  
  ! create new set of reciprocal lattice vectors in Microscope reference frame
  ! Note that reciprocal lattice vectors dot have two pi included, we are using the optical convention exp(i*g.r)
  RarVecM= TWOPI*CROSS(RbVecM,RcVecM)/DOT_PRODUCT(RbVecM,CROSS(RcVecM,RaVecM))
  RbrVecM= TWOPI*CROSS(RcVecM,RaVecM)/DOT_PRODUCT(RcVecM,CROSS(RaVecM,RbVecM))
  RcrVecM= TWOPI*CROSS(RaVecM,RbVecM)/DOT_PRODUCT(RaVecM,CROSS(RbVecM,RcVecM))
  
END SUBROUTINE ReciprocalLattice

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE UniqueAtomPositions(IErr)
!---------------------------------------------------------------
!Calculates the FULL set of possible fractional atomic positions
!And then gets rid of duplicates
!---------------------------------------------------------------
  
  USE WriteToScreen
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,ind,jnd,knd
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: IAllAtomicNumber, RAllAnisoDW
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RAllAtomPosition
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: RAllOccupancy, RAllIsoDW
  LOGICAL :: Lunique
  CHARACTER*2, DIMENSION(:), ALLOCATABLE :: SAllAtomName 
  CHARACTER*100 RAtomPositionString
  CHARACTER*100 indString

  CALL Message("UniqueAtomPositions",IMust,IErr)
  
  !All atom positions generated by symmetry, including duplicates
  ALLOCATE(RAllAtomPosition(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1),ITHREE),STAT=IErr)  
  ALLOCATE(SAllAtomName(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)),STAT=IErr)  
  ALLOCATE(RAllOccupancy(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)),STAT=IErr)  
  ALLOCATE(RAllIsoDW(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)),STAT=IErr)  
  ALLOCATE(IAllAtomicNumber(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)),STAT=IErr)  
  ALLOCATE(RAllAnisoDW(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)),STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"UniqueAtomPositions(",my_rank,")error in allocations"
     RETURN
  END IF
 
  !Apply symmetry elements to generate all equivalent positions 
  knd=1
  DO ind=1, SIZE(RSymVec,1)  
    DO jnd=1, SIZE(RBasisAtomPosition,1)
      RAllAtomPosition(knd,:)= MATMUL(RSymMat(ind,:,:),RBasisAtomPosition(jnd,:)) &
             + RSymVec(ind,:)
      SAllAtomName(knd) = SBasisAtomName(jnd)
      RAllOccupancy(knd) = RBasisOccupancy(jnd)
      RAllIsoDW(knd) = RBasisIsoDW(jnd)
      IAllAtomicNumber(knd) = IBasisAtomicNumber(jnd)
      RAllAnisoDW(knd) = IBasisAnisoDW(jnd)
	  knd=knd+1
    ENDDO
  ENDDO
  RAllAtomPosition=MODULO(RAllAtomPosition,ONE)
  WHERE(ABS(RAllAtomPosition).LT.TINY) RAllAtomPosition = ZERO 

  !--------------------------------------------------------------------  
  ! Now reduce to the set of unique fractional atomic positions, used to be subroutine CrystalUniqueFractionalAtomicPostitionsCalculation
  
  !first atom has to be in this set
  RAtomPosition(1,:)= RAllAtomPosition(1,:)
  SAtomName(1)= SAllAtomName(1)
  RIsoDW(1) = RAllIsoDW(1)
  ROccupancy(1) = RAllOccupancy(1)
  IAtomicNumber(1) = IAllAtomicNumber(1)
  RAnisoDW(1) = RAllAnisoDW(1)
  jnd=2
  !work through all possible atom coords and check for duplicates
  DO ind=2,IMaxPossibleNAtomsUnitCell
    Lunique=.TRUE.
    DO knd=1,jnd-1!check against the unique ones found so far
	  !PRINT*,ind,jnd,SUM(ABS(RAllAtomPosition(ind,:)-RAtomPosition(knd,:)))
      IF (SUM(ABS(RAllAtomPosition(ind,:)-RAtomPosition(knd,:))).LE.TINY) THEN  !position is the same
        IF (SAllAtomName(ind).EQ.SAtomName(knd)) THEN !name is the same too, so not unique
		  Lunique=.FALSE.
		  EXIT
		ENDIF
	  ENDIF
    ENDDO
	IF (Lunique .EQV. .TRUE.) THEN
      RAtomPosition(jnd,:)= RAllAtomPosition(ind,:)
      SAtomName(jnd)= SAllAtomName(ind)!never used
      RIsoDW(jnd) = RAllIsoDW(ind)
      ROccupancy(jnd) = RAllOccupancy(ind)
      IAtomicNumber(jnd) = IAllAtomicNumber(ind)!
      RAnisoDW(jnd) = RAllAnisoDW(ind)
      jnd=jnd+1
    ENDIF
  ENDDO
  INAtomsUnitCell= jnd-1

  !Finished with these variables now
  DEALLOCATE(RAllAtomPosition,SAllAtomName,RAllOccupancy,RAllIsoDW,IAllAtomicNumber,RAllAnisoDW)
  IF( IErr.NE.0 ) THEN
     PRINT*,"UniqueAtomPositions(",my_rank,")error deallocating RAllAtomPositions etc"
     RETURN
  ENDIF
    
!!$  Display the above variables to the user - index stored as string for formatting using message subroutine
  CALL Message("UniqueAtomPositions",IInfo,IErr, &
       MessageVariable = "NAtomsUnitCell",IVariable = INAtomsUnitCell) 
  DO ind=1,INAtomsUnitCell     
        WRITE(RAtomPositionString,'(3(F5.3,1X))') RAtomPosition(ind,:)
        WRITE(indString,*)ind
        CALL Message("UniqueAtomPositions",IMoreInfo,IErr, &
              MessageVariable = "RAtomPosition("//TRIM(ADJUSTL(indString))//",:)",MessageString = RAtomPositionString)
        CALL Message("UniqueAtomPositions",IMoreInfo,IErr, &
              MessageVariable = "SAtomName("//TRIM(ADJUSTL(indString))//")",MessageString = SAtomName(ind))
        CALL Message("UniqueAtomPositions",IMoreInfo,IErr, &
              MessageVariable = "RIsoDW("//TRIM(ADJUSTL(indString))//")",RVariable = RIsoDW(ind))
        CALL Message("UniqueAtomPositions",IMoreInfo,IErr, &
              MessageVariable = "ROccupancy("//TRIM(ADJUSTL(indString))//")",RVariable = ROccupancy(ind))
        CALL Message("UniqueAtomPositions",IMoreInfo,IErr, &
              MessageVariable = "IAtomicNumber("//TRIM(ADJUSTL(indString))//")",IVariable = IAtomicNumber(ind))
  ENDDO
  DO ind=1,INAtomsUnitCell
     WRITE(indString,*)ind
     CALL Message("UniqueAtomPositions",IDebug+IMoreInfo,IErr, &
              MessageVariable = "RAnisoDW("//TRIM(ADJUSTL(indString))//")",IVariable =RAnisoDW(ind))
  END DO

  ! Calculate atomic position vectors RAtomCoordinate from Fractional Coordinates and Lattice Vectors
  ! In microscope reference frame, in Angstrom units
  DO ind=1,INAtomsUnitCell
    DO jnd=1,ITHREE
       RAtomCoordinate(ind,jnd)= RAtomPosition(ind,1)*RaVecM(jnd)+RAtomPosition(ind,2)*RbVecM(jnd)+RAtomPosition(ind,3)*RcVecM(jnd)
    ENDDO
  ENDDO

END SUBROUTINE UniqueAtomPositions