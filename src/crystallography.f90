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
  ! Note that reciprocal lattice vectors have two pi included!!	
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

  ! Transform from orthogonal reference frame to microscope reference frame
  
  ! RTmat transforms from crystal reference to orthogonal frame
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
  
  ! R?DirO vectors are reference vectors in orthogonal frame
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
    ! now transform from crystal reference frame to orthogonal and then to microscope frame

  !RXDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RXDirC))
  !RYDirM = MATMUL(RTMatO2M,RYDirO)
  !RZDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RZDirC))
  
  !RB This is never used, but should be?
  RNormDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RNormDirC))
  RNormDirM = RNormDirM/SQRT(DOT_PRODUCT(RNormDirM,RNormDirM)) 
  
  ! now transform from crystal reference frame to orthogonal and then to microscope frame

  ! since RVec is already in orthogonal frame, only transform into microscope frame needed
  RaVecM= MATMUL(RTMatO2M,RaVecO)
  RbVecM= MATMUL(RTMatO2M,RbVecO)
  RcVecM= MATMUL(RTMatO2M,RcVecO)
  
  ! create new set of reciprocal lattice vectors
  RarVecM= TWOPI*CROSS(RbVecM,RcVecM)/DOT_PRODUCT(RbVecM,CROSS(RcVecM,RaVecM))
  RbrVecM= TWOPI*CROSS(RcVecM,RaVecM)/DOT_PRODUCT(RcVecM,CROSS(RaVecM,RbVecM))
  RcrVecM= TWOPI*CROSS(RaVecM,RbVecM)/DOT_PRODUCT(RaVecM,CROSS(RbVecM,RcVecM))
  
END SUBROUTINE ReciprocalLattice

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE AllAtomPositions(IErr)
!---------------------------------------------------------------
!Calculates the FULL set of possible fractional atomic positions
!---------------------------------------------------------------
  
  USE WriteToScreen
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,ind,jnd,knd,Ifullind,Iuniind
  REAL(RKIND):: norm
  LOGICAL :: Lunique
  CHARACTER*100 MNPString
  CHARACTER*100 indString

  CALL Message("AllAtomPositions",IMust,IErr)
  
  ALLOCATE(RAtomPosition(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1),ITHREE),STAT=IErr)  
  ALLOCATE(SAtomName(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)),STAT=IErr)  
  ALLOCATE(ROccupancy(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)),STAT=IErr)  
  ALLOCATE(RIsoDW(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)),STAT=IErr)  
  ALLOCATE(IAtomicNumber(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)),STAT=IErr)  
  ALLOCATE(IAnisoDW(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)),STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"AllAtomPositions(",my_rank,")error in allocations"
     RETURN
  END IF
 
  !Apply symmetry elements to generate equivalent positions 
  knd=1
  DO ind=1, SIZE(RSymVec,1)  
    DO jnd=1, SIZE(RBasisAtomPosition,1)
      RAtomPosition(knd,:)= MATMUL(RSymMat(ind,:,:),RBasisAtomPosition(jnd,:)) &
             + RSymVec(ind,:)
      SAtomName(knd) = SBasisAtomName(jnd)
      ROccupancy(knd) = RBasisOccupancy(jnd)
      RIsoDW(knd) = RBasisIsoDW(jnd)
      IAtomicNumber(knd) = IBasisAtomicNumber(jnd)
      IAnisoDW(knd) = IBasisAnisoDW(jnd)
	  knd=knd+1
    ENDDO
  ENDDO
  RAtomPosition=MODULO(RAtomPosition,ONE)
  WHERE(ABS(RAtomPosition).LT.TINY) RAtomPosition = ZERO 

  !--------------------------------------------------------------------  
  ! Now reduce to the set of unique fractional atomic positions, used to be subroutine CrystalUniqueFractionalAtomicPostitionsCalculation
  CALL Message("AllAtomPositions",IMust,IErr)
  
  !Initial value of total atoms/unit cell, set to the maximum and reduce to correct value at the end
  IMaxPossibleNAtomsUnitCell=SIZE(RBasisAtomPosition,1)*SIZE(RSymVec,1)

  ALLOCATE(MNP(IMaxPossibleNAtomsUnitCell,ITHREE),STAT=IErr)
  ALLOCATE(SMNP(IMaxPossibleNAtomsUnitCell),STAT=IErr)
  ALLOCATE(RDWF(IMaxPossibleNAtomsUnitCell),STAT=IErr)
  ALLOCATE(ROcc(IMaxPossibleNAtomsUnitCell),STAT=IErr)
  ALLOCATE(IAtoms(IMaxPossibleNAtomsUnitCell),STAT=IErr)
  ALLOCATE(IAnisoDWFT(IMaxPossibleNAtomsUnitCell),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"AllAtomPositions(",my_rank,")error in allocations"
     RETURN
  ENDIF

  !first atom
  MNP(1,:)= RAtomPosition(1,:)
  SMNP(1)= SAtomName(1)
  RDWF(1) = RIsoDW(1)
  ROcc(1) = ROccupancy(1)
  IAtoms(1) = IAtomicNumber(1)
  IAnisoDWFT(1) = IAnisoDW(1)
  jnd=2
  DO ind=2,IMaxPossibleNAtomsUnitCell!work through all possible atom coords
    Lunique=.TRUE.
    DO knd=1,jnd-1!check against the unique ones found so far
	  !PRINT*,ind,jnd,SUM(ABS(RAtomPosition(ind,:)-MNP(knd,:)))
      IF (SUM(ABS(RAtomPosition(ind,:)-MNP(knd,:))).LE.TINY) THEN  !position is the same
        IF (SAtomName(ind).EQ.SMNP(knd)) THEN !name is the same too, so not unique
		  Lunique=.FALSE.
		  EXIT
		ENDIF
	  ENDIF
    ENDDO
	IF (Lunique .EQV. .TRUE.) THEN
      MNP(jnd,:)= RAtomPosition(ind,:)
      SMNP(jnd)= SAtomName(ind)!never used
      RDWF(jnd) = RIsoDW(ind)
      ROcc(jnd) = ROccupancy(ind)
      IAtoms(jnd) = IAtomicNumber(ind)!
      IAnisoDWFT(jnd) = IAnisoDW(ind)
      jnd=jnd+1
    ENDIF
  ENDDO
  INAtomsUnitCell= jnd-1
  
!!$  Display the above variables to the user - index stored as string for formatting using message subroutine
  CALL Message("AllAtomPositions",IInfo,IErr, &
       MessageVariable = "NAtomsUnitCell",IVariable = INAtomsUnitCell) 
  DO ind=1,INAtomsUnitCell     
        WRITE(MNPString,'(3(F5.3,1X))') MNP(ind,:)
        WRITE(indString,*)ind
        CALL Message("AllAtomPositions",IMoreInfo,IErr, &
              MessageVariable = "MNP("//TRIM(ADJUSTL(indString))//",:)",MessageString = MNPString)
        CALL Message("AllAtomPositions",IMoreInfo,IErr, &
              MessageVariable = "SMNP("//TRIM(ADJUSTL(indString))//")",MessageString = SMNP(ind))
        CALL Message("AllAtomPositions",IMoreInfo,IErr, &
              MessageVariable = "RDWF("//TRIM(ADJUSTL(indString))//")",RVariable = RDWF(ind))
        CALL Message("AllAtomPositions",IMoreInfo,IErr, &
              MessageVariable = "ROcc("//TRIM(ADJUSTL(indString))//")",RVariable = ROcc(ind))
        CALL Message("AllAtomPositions",IMoreInfo,IErr, &
              MessageVariable = "IAtoms("//TRIM(ADJUSTL(indString))//")",IVariable = IAtoms(ind))
  ENDDO
  DO ind=1,INAtomsUnitCell
     WRITE(indString,*)ind
     CALL Message("AllAtomPositions",IDebug+IMoreInfo,IErr, &
              MessageVariable = "IAnisoDWFT("//TRIM(ADJUSTL(indString))//")",IVariable =IAnisoDWFT(ind))
  END DO

  ALLOCATE(RAtomCoordinate(INAtomsUnitCell,ITHREE),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"AllAtomPositions(",my_rank,")error allocating RAtomCoordinate"
     RETURN
  ENDIF
  ! Calculate atomic position vectors r from Fractional Coordinates and Lattice Vectors
  ! in microscopy reference frame in Angstrom units
  DO ind=1,INAtomsUnitCell
    DO jnd=1,ITHREE
       RAtomCoordinate(ind,jnd)= MNP(ind,1)*RaVecM(jnd)+MNP(ind,2)*RbVecM(jnd)+MNP(ind,3)*RcVecM(jnd)
    ENDDO
  ENDDO

!zz an ambition to deallocate some local variables here
  DEALLOCATE(SMNP,MNP,STAT=IErr)
  DEALLOCATE(ROccupancy,RAtomPosition,SAtomName,RIsoDW,IAtomicNumber,IAnisoDW)
  IF( IErr.NE.0 ) THEN
     PRINT*,"AllAtomPositions(",my_rank,")error allocating RAtomCoordinate"
     RETURN
  ENDIF
  
END SUBROUTINE AllAtomPositions