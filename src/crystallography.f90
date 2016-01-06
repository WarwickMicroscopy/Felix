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

SUBROUTINE CrystalLatticeVectorDetermination(IErr)
  
  USE MyNumbers
  USE WriteToScreen
  
  USE RPara; USE IPara; USE SPara
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind

  REAL(RKIND) :: RTTest
  REAL(RKIND), DIMENSION(THREEDIM) :: &
       RXDirO, RYDirO, RZDirO, RYDirC
  REAL(RKIND), DIMENSION(THREEDIM,THREEDIM) :: &
       RTMatC2O,RTMatO2M

  CHARACTER*50 indString
  CHARACTER*400  RTMatString

  CALL Message("CrystalLatticeVectorDetermination",IMust,IErr)

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
  ENDIF

  IF(IDiffractionFLAG.EQ.0) THEN
     
     RTTest = DOT_PRODUCT(RaVecO/DOT_PRODUCT(RaVecO,RaVecO),RbVecO/DOT_PRODUCT(RbVecO,RbVecO))*&
          DOT_PRODUCT(RbVecO/DOT_PRODUCT(RbVecO,RbVecO),RcVecO/DOT_PRODUCT(RcVecO,RcVecO))*&
          DOT_PRODUCT(RcVecO/DOT_PRODUCT(RcVecO,RcVecO),RaVecO/DOT_PRODUCT(RaVecO,RaVecO))
     
     IF(SCAN(SSpaceGroupName,'rR').NE.0) THEN
        IF(ABS(RTTest).LT.TINY) THEN
           SSpaceGroupName = TRIM(ADJUSTL("V"))
           CALL Message("CrystalLatticeVectorDetermination",IMust,IErr, &
                MessageString = "Warning: Crystal is either Obverse or Reverse,")
           CALL Message("CrystalLatticeVectorDetermination",IMust,IErr, &
                MessageString = "Selection Rules are not Currently In place to determine the difference,")
           CALL Message("CrystalLatticeVectorDetermination",IMust,IErr, &
                MessageString = "felix will assume the crystal is Obverse")
        ELSE
           SSpaceGroupName=TRIM(ADJUSTL('P'))
           CALL Message("CrystalLatticeVectorDetermination",IMust,IErr, &
                MessageString = "Crystal is in Primitive setting (Rhombohedral axes)")
        END IF
     END IF
  END IF
  ! Set up Reciprocal Lattice Vectors: orthogonal reference frame in 1/Angstrom units
  ! Note that reciprocal lattice vectors have two pi included!!	
  RarVecO= TWOPI*CROSS(RbVecO,RcVecO)/DOT_PRODUCT(RbVecO,CROSS(RcVecO,RaVecO))
  RbrVecO= TWOPI*CROSS(RcVecO,RaVecO)/DOT_PRODUCT(RcVecO,CROSS(RaVecO,RbVecO))
  RcrVecO= TWOPI*CROSS(RaVecO,RbVecO)/DOT_PRODUCT(RaVecO,CROSS(RbVecO,RcVecO))

  DO ind=1,THREEDIM
     IF (abs(RarVecO(ind)).lt.TINY) THEN
        RarVecO(ind) = ZERO
     ENDIF
     IF (abs(RbrVecO(ind)).lt.TINY) THEN
        RbrVecO(ind) = ZERO
     ENDIF
     IF (abs(RcrVecO(ind)).lt.TINY) THEN
        RcrVecO(ind) = ZERO
     ENDIF
  ENDDO

  ! Transform from orthogonal reference frame to microscope reference frame
  
  ! RTmat transforms from crystal reference to orthogonal frame
  RTMatC2O(:,1)= RaVecO(:)
  RTMatC2O(:,2)= RbVecO(:)
  RTMatC2O(:,3)= RcVecO(:)

  CALL Message("CrystalLatticeVectorDetermination",IMoreInfo,IErr,MessageString = "RTMatC2O")
  
  DO ind=1,THREEDIM
     WRITE(indString,*)ind
     WRITE(RTMatString,'(3(F8.3,1X))') RTMatC2O(:,ind)
     CALL Message("CrystalLatticeVectorDetermination",IMoreInfo,IErr, &
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

  CALL Message("CrystalLatticeVectorDetermination",IMoreInfo,IErr,MessageString = "RTMatO2M")

  DO ind=1,THREEDIM
     WRITE(indString,*)ind
     WRITE(RTMatString,'(3(F8.3,1X))') RTMatO2M(:,ind)
     CALL Message("CrystalLatticeVectorDetermination",IMoreInfo,IErr, &
          MessageVariable = "RTMatO2M(:,"//TRIM(ADJUSTL(indString))//")", &
          MessageString = TRIM(ADJUSTL(RTMatString)))
  END DO
    ! now transform from crystal reference frame to orthogonal and then to microscope frame

  RXDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RXDirC))
  RYDirM = MATMUL(RTMatO2M,RYDirO)
  RZDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RZDirC))
  ! now transform from crystal reference frame to orthogonal and then to microscope frame

  RNormDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RNormDirC))
  RNormDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RNormDirC))
  RNormDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RNormDirC))

  ! since RVec is already in orthogonal frame, only transform into microscope frame needed
  RaVecM= MATMUL(RTMatO2M,RaVecO)
  RbVecM= MATMUL(RTMatO2M,RbVecO)
  RcVecM= MATMUL(RTMatO2M,RcVecO)
  
  ! create new set of reciprocal lattice vectors

  RarVecM= TWOPI*CROSS(RbVecM,RcVecM)/DOT_PRODUCT(RbVecM,CROSS(RcVecM,RaVecM))
  RbrVecM= TWOPI*CROSS(RcVecM,RaVecM)/DOT_PRODUCT(RcVecM,CROSS(RaVecM,RbVecM))
  RcrVecM= TWOPI*CROSS(RaVecM,RbVecM)/DOT_PRODUCT(RaVecM,CROSS(RbVecM,RcVecM))
  
END SUBROUTINE CrystalLatticeVectorDetermination

!---------------------------------------------------------------
!Calculates the FULL set of possible fractional atomic positions
!---------------------------------------------------------------
SUBROUTINE CrystalFullFractionalAtomicPostitionsCalculation(IErr)
  
  USE WriteToScreen
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       IErr, ind,jnd,knd, Ifullind

  CALL Message("CrystalFullFractionalAtomicPostitionsCalculation",IMust,IErr)
  
  ALLOCATE( &
       RFullAtomicFracCoordVec( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1),&
       THREEDIM), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystalFullFractionalAtomicPostitionsCalculation(",&
          my_rank, ") error ", IErr, " in ALLOCATE RFullAtomicFracCoordVec"
     RETURN
  ENDIF
  
  ALLOCATE( &
       SFullAtomicNameVec( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystalFullFractionalAtomicPostitionsCalculation(",&
          my_rank, ") error ", IErr, " in ALLOCATE SFullAtomicNameVec"
     RETURN
  ENDIF
  
  ALLOCATE( &
       RFullPartialOccupancy( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystalFullFractionalAtomicPostitionsCalculation(",&
          my_rank, ") error ", IErr, " in ALLOCATE RFullPartialOccupancy"
     RETURN
  ENDIF
  
  ALLOCATE( &
       RFullIsotropicDebyeWallerFactor( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystalFullFractionalAtomicPostitionsCalculation(",&
          my_rank, ") error ", IErr, " in ALLOCATE RFullIsotropicDebyeWallerFactor"
     RETURN
  ENDIF
  ALLOCATE( &
       IFullAtomNumber( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystalFullFractionalAtomicPostitionsCalculation(",&
          my_rank, ") error ", IErr, " in ALLOCATE IFullAtomNumber"
     RETURN
  ENDIF
  
  ALLOCATE( &
       IFullAnisotropicDWFTensor( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystalFullFractionalAtomicPostitionsCalculation(",&
          my_rank, ") error ", IErr, " in ALLOCATE IFullAnisotropicDWFTensor"
     RETURN
  ENDIF
  
  DO ind=1, SIZE(RSymVec,DIM=1)
     
     DO jnd=1, SIZE(RAtomSiteFracCoordVec,DIM=1)
        Ifullind= SIZE(RSymVec,1)*(jnd-1) + ind
        
        RFullAtomicFracCoordVec(Ifullind,:)= &
             MATMUL(RSymMat(ind,:,:),RAtomSiteFracCoordVec(jnd,:)) &
             + RSymVec(ind,:)
        SFullAtomicNameVec(Ifullind) = SAtomName(jnd)
        RFullPartialOccupancy(Ifullind) = RAtomicSitePartialOccupancy(jnd)
        RFullIsotropicDebyeWallerFactor(Ifullind) = RIsotropicDebyeWallerFactors(jnd)
        IFullAtomNumber(Ifullind) = IAtomNumber(jnd)
        IFullAnisotropicDWFTensor(Ifullind) = IAnisotropicDWFTensor(jnd)
        
        ! renormalize such that all values are non-negative
        DO knd=1,THREEDIM
           IF( RFullAtomicFracCoordVec(Ifullind,knd) .LT. ZERO) THEN
              RFullAtomicFracCoordVec(Ifullind,knd)= &
                   RFullAtomicFracCoordVec(Ifullind,knd)+1.D0
           ENDIF
        ENDDO
        
     ENDDO
     
  ENDDO

  DO ind = 1,SIZE(RFullAtomicFracCoordVec,DIM=1)
     DO jnd = 1,SIZE(RFullAtomicFracCoordVec,DIM=2)
        IF (RFullAtomicFracCoordVec(ind,jnd).LT.ZERO) THEN
           RFullAtomicFracCoordVec(ind,jnd) = &
                RFullAtomicFracCoordVec(ind,jnd) + ONE
        END IF
        IF (RFullAtomicFracCoordVec(ind,jnd).GE.ONE) THEN
           RFullAtomicFracCoordVec(ind,jnd) = &
                RFullAtomicFracCoordVec(ind,jnd) - ONE
        END IF
     END DO
  END DO

  WHERE(ABS(RFullAtomicFracCoordVec).LT.TINY) RFullAtomicFracCoordVec = ZERO 

END SUBROUTINE CrystalFullFractionalAtomicPostitionsCalculation


SUBROUTINE CrystalUniqueFractionalAtomicPostitionsCalculation (IErr)


  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  
  IMPLICIT NONE
  
  REAL(RKIND):: &
       norm

  INTEGER(IKIND) :: &
       IErr, ind,jnd, Iuniind
  LOGICAL :: &
       Lunique

  CHARACTER*100 MNPString
  CHARACTER*100 indString
  
  CALL Message("CrystalUniqueFractionalAtomicPostitionsCalculation",IMust,IErr)

  ! Calculate the set of unique fractional atomic positions

  ALLOCATE( &
       MNP(ITotalAtoms,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystalUniqueFractionalAtomicPostitionsCalculation(", my_rank, ") error ", IErr, " in ALLOCATE MNP"
     RETURN
  ENDIF

  ALLOCATE( &
       SMNP(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystalUniqueFractionalAtomicPostitionsCalculation(", my_rank, ") error ", IErr, " in ALLOCATE SMNP"
     RETURN
  ENDIF

  ALLOCATE( &
       RDWF(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystalUniqueFractionalAtomicPostitionsCalculation(", my_rank, ") error ", IErr, " in ALLOCATE RDWF"
     RETURN
  ENDIF

  ALLOCATE( &
       ROcc(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystalUniqueFractionalAtomicPostitionsCalculation(", my_rank, ") error ", IErr, " in ALLOCATE ROcc"
     RETURN
  ENDIF

  ALLOCATE( &
       IAtoms(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystalUniqueFractionalAtomicPostitionsCalculation(", my_rank, ") error ", IErr, " in ALLOCATE IAtoms"
     RETURN
  ENDIF

  ALLOCATE( &
       IAnisoDWFT(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystalUniqueFractionalAtomicPostitionsCalculation(", my_rank, ") error ", IErr, " in ALLOCATE IAnisoDWFT"
     RETURN
  ENDIF

  MNP(1,:)= RFullAtomicFracCoordVec(1,:)
  SMNP(1)= SFullAtomicNameVec(1)
  RDWF(1) = RFullIsotropicDebyeWallerFactor(1)
  ROcc(1) = RFullPartialOccupancy(1)
  IAtoms(1) = IFullAtomNumber(1)
  IAnisoDWFT(1) = IFullAnisotropicDWFTensor(1)
  Iuniind=1

  DO ind=2,SIZE(RFullAtomicFracCoordVec,1)
     DO jnd=1,Iuniind
        IF ( RFullAtomicFracCoordVec(ind,1) .EQ. MNP(jnd,1) .AND. &
             RFullAtomicFracCoordVec(ind,2) .EQ. MNP(jnd,2) .AND. &
             RFullAtomicFracCoordVec(ind,3) .EQ. MNP(jnd,3) .AND. &
             SFullAtomicNameVec(ind) .EQ. SMNP(jnd) ) THEN
           !this seems NOT a unique coordinate
           Lunique=.FALSE.
           EXIT
        ENDIF
        Lunique=.TRUE.
     ENDDO

     IF(Lunique .EQV. .TRUE.) THEN
        Iuniind=Iuniind+1
        MNP(Iuniind,:)= RFullAtomicFracCoordVec(ind,:)
        SMNP(Iuniind)= SFullAtomicNameVec(ind)
        RDWF(Iuniind) = RFullIsotropicDebyeWallerFactor(ind)
        ROcc(Iuniind) = RFullPartialOccupancy(ind)
        IAtoms(Iuniind) = IFullAtomNumber(ind)
        IAnisoDWFT(Iuniind) = IFullAnisotropicDWFTensor(ind)
        
     ENDIF
     
  ENDDO

!!$  Display the above variables to the user - index stored as string for formatting using message subroutine
  DO ind=1,Iuniind     

        WRITE(MNPString,'(3(F5.3,1X))') MNP(ind,:)
        WRITE(indString,*)ind
        CALL Message("CrystalUniqueFractionalAtomicPostitionsCalculation",IMoreInfo,IErr, &
              MessageVariable = "MNP("//TRIM(ADJUSTL(indString))//",:)",MessageString = MNPString)
        CALL Message("CrystalUniqueFractionalAtomicPostitionsCalculation",IMoreInfo,IErr, &
              MessageVariable = "SMNP("//TRIM(ADJUSTL(indString))//")",MessageString = SMNP(ind))
        CALL Message("CrystalUniqueFractionalAtomicPostitionsCalculation",IMoreInfo,IErr, &
              MessageVariable = "RDWF("//TRIM(ADJUSTL(indString))//")",RVariable = RDWF(ind))
        CALL Message("CrystalUniqueFractionalAtomicPostitionsCalculation",IMoreInfo,IErr, &
              MessageVariable = "ROcc("//TRIM(ADJUSTL(indString))//")",RVariable = ROcc(ind))
        CALL Message("CrystalUniqueFractionalAtomicPostitionsCalculation",IMoreInfo,IErr, &
              MessageVariable = "IAtoms("//TRIM(ADJUSTL(indString))//")",IVariable = IAtoms(ind))
    
        
        !PRINT *,"Crystallography : ",MNP(ind,:),SMNP(ind),RDWF(ind),ROcc(ind),IAtoms(ind),ind
   !  END IF
  ENDDO

  !IAnisoDWFT now contains a list of indices referring to the correct Anisotropic Debye Wall Factor Tensor for each atom in the unit cell

  DO ind=1,ITotalAtoms

     WRITE(indString,*)ind
     CALL Message("CrystalUniqueFractionalAtomicPostitionsCalculation",IDebug+IMoreInfo,IErr, &
              MessageVariable = "IAnisoDWFT("//TRIM(ADJUSTL(indString))//")",IVariable =IAnisoDWFT(ind))
  END DO

  ! Calculate atomic position vectors r from Fractional Coordinates and Lattice Vectors
  ! in microscopy reference frame in Angstrom units
  
  DO ind=1,SIZE(MNP,DIM=1)
     DO jnd=1,THREEDIM
        RrVecMat(ind,jnd)= MNP(ind,1)*RaVecM(jnd) + MNP(ind,2)*RbVecM(jnd) + MNP(ind,3)*RcVecM(jnd)
     ENDDO
  ENDDO
  
  INAtomsUnitCell= SIZE(MNP,DIM=1)
  
  CALL Message("CrystalUniqueFractionalAtomicPostitionsCalculation",IInfo,IErr, &
       MessageVariable = "NAtomsUnitCell",IVariable = INAtomsUnitCell)
  
  DO ind = 1,INAtomsUnitCell
     DO jnd = 1,THREEDIM
        IF (ABS(RrVecMat(ind,jnd)).LE.TINY) THEN
           RrVecMat(ind,jnd) = ZERO
        ENDIF
     ENDDO
  ENDDO

  ! select the bases in the reciprocal space

  RETURN

END SUBROUTINE  CrystalUniqueFractionalAtomicPostitionsCalculation

