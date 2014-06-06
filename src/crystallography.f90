!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! FelixSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
!
! (C) 2013/14, all rights reserved
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This file is part of FelixSim.
!
!  FelixSim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  FelixSim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with FelixSim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE Crystallography( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
   USE MyMPI

  IMPLICIT NONE

  REAL(RKIND) norm
  REAL(RKIND), DIMENSION(THREEDIM) :: &
       RXDirO, RYDirO, RZDirO, RYDirC
  REAL(RKIND), DIMENSION(THREEDIM,THREEDIM) :: &
       RTMatC2O,RTMatO2M

  INTEGER IErr, ind,jnd,knd, Ifullind, Iuniind
  LOGICAL Lunique
  REAL(RKIND) :: RTTEst

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Crystallography(",my_rank,")"
  END IF
  
  IF(my_rank.EQ.0) THEN
     PRINT*,RZDirC,RXDirC
  END IF
  
  IF(IDiffractionFLAG.EQ.1) THEN
     DEALLOCATE(&
          RFullAtomicFracCoordVec,SFullAtomicNameVec,&
          RFullPartialOccupancy,RFullIsotropicDebyeWallerFactor, &
          IFullAtomNumber, IFullAnisotropicDWFTensor, &
          MNP,SMNP,RDWF,ROcc,IAtoms,IAnisoDWFT,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
        RETURN
     ENDIF
  END IF
  
  IF(IPseudoCubicFLAG.EQ.1.AND.IDiffractionFLAG.EQ.0) THEN
     RZDirC(1) = (IIncidentBeamDirectionX/3.0D0)+&
          (IIncidentBeamDirectionY/3.0D0)-&
          (2.0D0*(IIncidentBeamDirectionZ/3.0D0))
     RZDirC(2) = (-IIncidentBeamDirectionX/3.0D0)+&
          (2.0D0*(IIncidentBeamDirectionY/3.0D0))-&
          (IIncidentBeamDirectionZ/3.0D0)
     RZDirC(3) = (IIncidentBeamDirectionX/3.0D0)+&
          (IIncidentBeamDirectionY/3.0D0)+&
          (IIncidentBeamDirectionZ/3.0D0)
     RXDirC(1) = (IXDirectionX/3.0D0)+&
          (IXDirectionY/3.0D0)-&
          (2.0D0*(IXDirectionZ/3.0D0))
     RXDirC(2) = (-IXDirectionX/3.0D0)+&
          (2.0D0*(IXDirectionY/3.0D0))-&
          (IXDirectionZ/3.0D0)
     RXDirC(3) = (IXDirectionX/3.0D0)+&
          (IXDirectionY/3.0D0)+&
          (IXDirectionZ/3.0D0)
     RNormDirC(1) = (INormalDirectionX/3.0D0)+&
          (INormalDirectionY/3.0D0)-&
          (2.0D0*(INormalDirectionZ/3.0D0))
     RNormDirC(2) = (-INormalDirectionX/3.0D0)+&
          (2.0D0*(INormalDirectionY/3.0D0))-&
          (INormalDirectionZ/3.0D0)
     RNormDirC(3) = (INormalDirectionX/3.0D0)+&
          (INormalDirectionY/3.0D0)+&
          (INormAlDirectionZ/3.0D0)

     
     IF(my_rank.EQ.0) THEN
        PRINT*,RZDirC,RXDirC
     END IF
     
     DO ind =1,3
        IF(ABS(RZDirC(ind)).LE.TINY) THEN
           RZDirC(ind) = 100000000.0D0 ! A large number
        END IF
        IF(ABS(RXDirC(ind)).LE.TINY) THEN
           RXDirC(ind) = 100000000.0D0 ! A large number
        END IF
        IF(ABS(RNormDirC(ind)).LE.TINY) THEN
           RNormDirC(ind) = 100000000.0D0 ! A large number
        END IF
     END DO
        RZDirC = RZDirC/MINVAL(ABS(RZDirC))
        RXDirC = RXDirC/MINVAL(ABS(RXDirC))

     DO ind =1,3
        IF(RZDirC(ind).GT.10000000.0D0) THEN
           RZDirC(ind) = ZERO ! A large number
        END IF
        IF(RXDirC(ind).GT.10000000.0D0) THEN
           RXDirC(ind) = ZERO ! A large number
        END IF
        IF(RNormDirC(ind).GT.10000000.0D0) THEN
           RNormDirC(ind) = ZERO ! A large number
        END IF
     END DO
     
     RZDirC = REAL(NINT(RZDirC))
     RXDirC = REAL(NINT(RXDirC))
     RNormDirC = REAL(NINT(RNormDirC))
     
  END IF

  IF(my_rank.EQ.0) THEN
     PRINT*,RZDirC,RXDirC
  END IF

  !RYDirC = CROSS(RZDirC,RYDirC)

  ! Setup Crystal Lattice Vectors: orthogonal reference frame in Angstrom units

  RaVecO(1)= RLengthX
  RaVecO(2)= 0.D0
  RaVecO(3)= 0.D0

  RbVecO(1)= RLengthY*COS(RGamma)
  RbVecO(2)= RLengthY*SIN(RGamma)
  RbVecO(3)= 0.D0

  RcVecO(1)= RLengthZ*COS(RBeta)
  RcVecO(2)= RLengthZ*(COS(RAlpha)-COS(RBeta)*COS(RGamma))/SIN(RGamma)
  RcVecO(3)= RLengthZ*( &
       SQRT(1.D0- &
        COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)-COS(RGamma)*COS(RGamma) + &
        2.D0*COS(RAlpha)*COS(RBeta)*COS(RGamma)) / &
       SIN(RGamma))

  IF(IVolumeFLAG .EQ. 0) THEN
     RVolume= RLengthX*RLengthY*RLengthZ* &
          SQRT(1.0D0 - &
          COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)-COS(RGamma)*COS(RGamma) + &
          2.0D0*COS(RAlpha)*COS(RBeta)*COS(RGamma))
  ENDIF

  IF(IDiffractionFLAG.EQ.0) THEN
     
     RTTest = DOT_PRODUCT(RaVecO/DOT_PRODUCT(RaVecO,RaVecO),RbVecO/DOT_PRODUCT(RbVecO,RbVecO))*&
          DOT_PRODUCT(RbVecO/DOT_PRODUCT(RbVecO,RbVecO),RcVecO/DOT_PRODUCT(RcVecO,RcVecO))*&
          DOT_PRODUCT(RcVecO/DOT_PRODUCT(RcVecO,RcVecO),RaVecO/DOT_PRODUCT(RaVecO,RaVecO))
     
     IF(ABS(RTTest).LT.TINY.AND.SCAN(SSpaceGroupName,'rR').NE.0) THEN
        SSpaceGroupName = TRIM(ADJUSTL("V"))
        IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"Crystal is Obverse"
        END IF
     ELSE
        IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"Crystal is Reverse"
        END IF
     END IF
  END IF

  ! Set up Reciprocal Lattice Vectors: orthogonal reference frame in 1/Angstrom units

  RarVecO= TWOPI*CROSS(RbVecO,RcVecO)/DOT(RbVecO,CROSS(RcVecO,RaVecO))
  RbrVecO= TWOPI*CROSS(RcVecO,RaVecO)/DOT(RcVecO,CROSS(RaVecO,RbVecO))
  RcrVecO= TWOPI*CROSS(RaVecO,RbVecO)/DOT(RaVecO,CROSS(RbVecO,RcVecO))

  DO ind=1,THREEDIM
     IF (abs(RarVecO(ind)).lt.1.D-3) THEN
        RarVecO(ind) = 0.0D0
     ENDIF
     IF (abs(RbrVecO(ind)).lt.1.D-3) THEN
        RbrVecO(ind) = 0.0D0
     ENDIF
     IF (abs(RcrVecO(ind)).lt.1.D-3) THEN
        RcrVecO(ind) = 0.0D0
     ENDIF
  ENDDO

  ! Transform from orthogonal reference frame to microscope reference frame
  
  ! RTmat transforms from crystal reference to orthogonal frame
  RTMatC2O(:,1)= RaVecO(:)
  RTMatC2O(:,2)= RbVecO(:)
  RTMatC2O(:,3)= RcVecO(:)

  IF((IWriteFLAG.GE.3.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Crystallography(",my_rank,") : RTMatC2O=", RTMatC2O
  END IF

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
  
  IF((IWriteFLAG.GE.3.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Crystallography(",my_rank,"): RTMatO2M=", RTMatO2M
  END IF

  ! now transform from crystal reference frame to orthogonal and then to microscope frame

  RXDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RXDirC))
  RYDirM = MATMUL(RTMatO2M,RYDirO)
  RZDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RZDirC))
  ! now transform from crystal reference frame to orthogonal and then to microscope frame

  RNormDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RNormDirC))
  RNormDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RNormDirC))
  RNormDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RNormDirC))

  ! since R?Vec is already in orthogonal frame, only transform into microscope frame needed
  RaVecM= MATMUL(RTMatO2M,RaVecO)
  RbVecM= MATMUL(RTMatO2M,RbVecO)
  RcVecM= MATMUL(RTMatO2M,RcVecO)
  
  ! create new set of reciprocal lattice vectors

  RarVecM= TWOPI*CROSS(RbVecM,RcVecM)/DOT(RbVecM,CROSS(RcVecM,RaVecM))
  RbrVecM= TWOPI*CROSS(RcVecM,RaVecM)/DOT(RcVecM,CROSS(RaVecM,RbVecM))
  RcrVecM= TWOPI*CROSS(RaVecM,RbVecM)/DOT(RaVecM,CROSS(RbVecM,RcVecM))

  ! Calculate the FULL set of possible fractional atomic positions

  ALLOCATE( &
       RFullAtomicFracCoordVec( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1),&
       THREEDIM), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       SFullAtomicNameVec( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       RFullPartialOccupancy( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       RFullIsotropicDebyeWallerFactor( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE( &
       IFullAtomNumber( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       IFullAnisotropicDWFTensor( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
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
        RFullIsotropicDebyeWallerFactor = RIsotropicDebyeWallerFactors(jnd)
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

  ! Calculate the set of unique fractional atomic positions
  
  ALLOCATE( &
       MNP(ITotalAtoms,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       SMNP(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       RDWF(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       ROcc(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       IAtoms(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       IAnisoDWFT(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
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
!!$        IF(IWriteFLAG.GE.1) THEN
!!$           PRINT*,"MNP(Iuniind,:) = ",MNP(Iuniind,:),SMNP(Iuniind)
!!$        END IF
        RDWF(Iuniind) = RFullIsotropicDebyeWallerFactor(ind)
        ROcc(Iuniind) = RFullPartialOccupancy(ind)
        IAtoms(Iuniind) = IFullAtomNumber(ind)
        IAnisoDWFT(Iuniind) = IFullAnisotropicDWFTensor(ind)
        
     ENDIF
     
  ENDDO
  
  DO ind=1,Iuniind     
     IF((IWriteFLAG.GE.3.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        IF (ind.EQ.1) THEN
           PRINT*,"DBG : Crystallography : ","MNP ","SMNP ","RDWF ","ROcc "
        END IF
        
        PRINT *,"DBG: Crystallography : ",MNP(ind,:),SMNP(ind),RDWF(ind),ROcc(ind),IAtoms(ind),ind
     END IF
  ENDDO

  !IAnisoDWFT now contains a list of indices referring to the correct Anisotropic Debye Wall Factor Tensor for each atom in the unit cell

  DO ind=1,ITotalAtoms
     
     IF((IWriteFLAG.GE.3.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"DBG: IAnisoDWFT",IAnisoDWFT(ind)
     END IF
     
  END DO
  
  ! Calculate atomic position vectors r from Fractional Coordinates and Lattice Vectors
  ! in microscopy reference frame in Angstrom units
  
  DO ind=1,SIZE(MNP,DIM=1)
     DO jnd=1,THREEDIM
        RrVecMat(ind,jnd)= MNP(ind,1)*RaVecM(jnd) + MNP(ind,2)*RbVecM(jnd) + MNP(ind,3)*RcVecM(jnd)
     ENDDO
  ENDDO
  
  INAtomsUnitCell= SIZE(MNP,DIM=1)
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: main(", my_rank, ") NAtomsUnitCell=", INAtomsUnitCell
  END IF
  
  DO ind = 1,INAtomsUnitCell
     DO jnd = 1,THREEDIM
        IF (ABS(RrVecMat(ind,jnd)).LE.TINY) THEN
           RrVecMat(ind,jnd) = ZERO
        ENDIF
     ENDDO
  ENDDO

  ! select the bases in the reciprocal space

  RETURN

END SUBROUTINE Crystallography
