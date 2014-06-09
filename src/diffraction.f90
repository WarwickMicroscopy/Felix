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

SUBROUTINE DiffractionPatternDefinitions( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE SPara
  USE CPara
  USE IChannels
  USE BlochPara
  
  USE MyMPI

  IMPLICIT NONE

  REAL(RKIND) :: &
       norm, dummyVec(THREEDIM),dummy
  INTEGER IErr, ind,jnd,icheck,ihklrun
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DiffractionPatternDefinitions()"
  END IF

  icheck = 0
  
  IHKLMAXValue = 15
  ihklrun = 0
  DO WHILE (icheck.EQ.0)
     ihklrun = ihklrun+1
     
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"DiffractionPatternDefinitions(",my_rank,") hklrun = ",ihklrun
     END IF
     
     CALL NewHKLMake(IHKLMAXValue,RZDirC,TWOPI/180.0D0,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error in NewHKLMake()"
        RETURN
     ENDIF
     
     IF(SIZE(RHKl,DIM=1).LT.IMinReflectionPool) THEN
        IHKLMAXValue = IHKLMAXValue*2
        Deallocate(RHKL,STAT=ierr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error ", IErr, &
                " in DEALLOCATE() of DYNAMIC variables RHKL"
           RETURN
        ENDIF
        
        CYCLE
        
     ELSE
        icheck = 1
     END IF
     
     CALL ReSortHKL( RHKL, SIZE(RHKL,1))
     IF( IErr.NE.0 ) THEN
        PRINT*,"DiffractionPatternDefintions(): error in ReSortHKL()"
        RETURN
     ENDIF     
     
     IF(IXDirectionFLAG.EQ.0) THEN
        IDiffractionFLAG = 1
        RXDirC(1) = RHKL(2,1)
        RXDirC(2) = RHKL(2,2)
        RXDirC(3) = RHKL(2,3)
        CALL Crystallography( IErr )
        IF( IErr.NE.0 ) THEN
           PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error in Crystallography()"
           RETURN
        ENDIF
        
     END IF
     
     ALLOCATE(&
          RgVecMatT(SIZE(RHKL,DIM=1),THREEDIM), &
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables RgVecMatT(HKL)"
        RETURN
     ENDIF
     
     ALLOCATE(&
          RgVecMag(SIZE(RHKL,DIM=1)), &
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables RgVecMag(HKL)"
        RETURN
     ENDIF
     
     DO ind=1,SIZE(RHKL,DIM=1)
        DO jnd=1,THREEDIM
           RgVecMatT(ind,jnd)= &
                RHKL(ind,1)*RarVecM(jnd) + &
                RHKL(ind,2)*RbrVecM(jnd) + &
                RHKL(ind,3)*RcrVecM(jnd)
        ENDDO
     ENDDO
     
     ! G vector magnitudes in 1/Angstrom units
     
     DO ind=1,SIZE(RHKL,DIM=1)
        RgVecMag(ind)= SQRT(DOT_PRODUCT(RgVecMatT(ind,:),RgVecMatT(ind,:)))
     ENDDO
  
     RBSMaxGVecAmp = RgVecMag(IMinReflectionPool)
     
     nReflections = 0
     nStrongBeams = 0
     nWeakBeams = 0
     
     DO ind=1,SIZE(RHKL,DIM=1)
        IF (ABS(RgVecMag(ind)).LE.RBSMaxGVecAmp) THEN
           nReflections = nReflections + 1
        ENDIF
     ENDDO
     
  END DO

  IF((IWriteFLAG.GE.4.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     DO ind=1,SIZE(RHKL,DIM=1)
        PRINT*,RHKL(ind,:)
     END DO
  END IF

  ALLOCATE(&
       RGn(SIZE(RHKL,DIM=1)), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RGn(HKL)"
     RETURN
  ENDIF

  RNormDirM = RNormDirM/sqrt(DOT_PRODUCT(RNormDirM,RNormDirM))

  DO ind =1,SIZE(RHKL,DIM=1)
     RGn(ind) = DOT_PRODUCT(RgVecMatT(ind,:),RNormDirM)
  END DO
  
  ALLOCATE(&
       RSg(SIZE(RHKL,DIM=1)), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RSg(HKL)"
     RETURN
  ENDIF

  IF(ICentralBeamFLAG.EQ.0) THEN
     !no central beam, 1st beam is used
     dummy= RHKL(1,1)*RHKL(1,1)+RHKL(1,2)*RHKL(1,2)+RHKL(1,3)*RHKL(1,3)
     dummy= SQRT(dummy)
  ELSE
     !there is a central beam, 2nd beam is used
     dummy= RHKL(2,1)*RHKL(2,1)+RHKL(2,2)*RHKL(2,2)+RHKL(2,3)*RHKL(2,3)
     dummy= SQRT(dummy)
  ENDIF

  RBraggCentral= RElectronWaveLength/(2.0D0*RLengthX)*dummy

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DiffractionPatternDefinitions (",my_rank,") BraggCentral=", RBraggCentral
  END IF


  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DiffractionPatternDefinitions (",my_rank,") GMax = ",RBSMaxGVecAmp
  END IF

  ! smallest g is gmag(2) IF 000 beam is included !!!add error catch here
  
  IF (ICentralBeamFLAG.EQ.1) THEN
     RMinimumGMag = RgVecMag(2)
  ELSE
     RMinimumGMag = RgVecMag(1)
  ENDIF

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DiffractionPatternDefinitions (",my_rank,") MinimumGMag = ",RMinimumGMag
  END IF

  ! Calculate the Angles to the centre of all the disks from the central disk
  
  DO ind=1,SIZE(RHKL,DIM=1)
     RSg(ind) = 2*RElectronWaveVectorMagnitude* &
          ((-2*RElectronWaveVectorMagnitude +&
          SQRT((2*RElectronWaveVectorMagnitude)**2 + &
          4*RgVecMag(ind)**2))/2)
  ENDDO
  
  ! Determine the number of Gs within GMax (input variable) 
  ! Furthermore, determine which Gs have Sg < SgMax (Strong Beams)
  
  IF (nReflections.LT.IReflectOut) THEN

     IReflectOut = nReflections

  END IF

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"DiffractionPatternDefinitions (",my_rank,") No. of Reflections = ",nReflections

  END IF
  ! resolution in k space
  RDeltaK = RMinimumGMag*RConvergenceAngle/IPixelCount

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DiffractionPatternDefinitions (",my_rank,") DeltaK = ", RDeltaK
  END IF

  RETURN

END SUBROUTINE DiffractionPatternDefinitions

SUBROUTINE NewHKLmake(Ihklmax,Rhkl0Vec,RAcceptanceAngle,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  USE MPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) IErr, Ihklmax,ind,jnd,knd,INhkl
  REAL(RKIND) RAcceptanceAngle
  REAL(RKIND), DIMENSION(THREEDIM) :: Rhkl0Vec,RhklDummyUnitVec,RhklDummyVec,Rhkl0UnitVec

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"NewHKLmake(",my_rank,")"
  END IF
  INhkl = 0

  Rhkl0UnitVec= Rhkl0Vec/SQRT(DOT_PRODUCT(REAL(Rhkl0Vec,RKIND),REAL(Rhkl0Vec,RKIND)))
  
  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1
           
           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)
           
           RhklDummyUnitVec= RhklDummyVec / &
                SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))
           !-NINT(MOD(RhklDummyVec(1)+RhklDummyVec(2),2.D0))
           SELECT CASE(SSpaceGroupName)
           CASE("F") !Face Centred
              IF(((ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),2.D0)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),2.D0)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),2.D0)).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),2.D0))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),2.D0))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),2.D0))).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),2.D0))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),2.D0))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),2.D0))).GT.TINY))) THEN
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                    
                 ENDIF
                 ! INhkl = INhkl + 1
              END IF              
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),2.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),2.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                    
                 ENDIF
              END IF
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),2.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                    
                 ENDIF
              END IF
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),2.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),3.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),3.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
           CASE("P")! Primitive
              !INhkl = INhkl + 1
              
              IF(IZolzFLAG.EQ.1) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                   .LE.SIN(RAcceptanceAngle)) THEN
                 INhkl = INhkl +1       
              ENDIF
           CASE DEFAULT
              PRINT*,"HKLMake(): unknown space group", SSpaceGroupName, "--- aborting"
              IErr=1
           END SELECT

        END DO
     END DO
  END DO
  
  Allocate(&
       RHKL((INhkl+1),THREEDIM),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"hklMake(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Rhkl"
     RETURN
  ENDIF
  
  INhkl = 0

  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1

           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)

           RhklDummyUnitVec= RhklDummyVec / &
                SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))
           
           SELECT CASE(SSpaceGroupName)
           CASE("F") !Face Centred
              IF(((ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),2.D0)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),2.D0)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),2.D0)).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),2.D0))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),2.D0))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),2.D0))).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),2.D0))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),2.D0))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),2.D0))).GT.TINY))) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),2.0D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),2.0D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),2.0D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),2.0D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),3.0D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),3.0D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("P")! Primitive
              !INhkl = INhkl + 1
              !RHKL(INhkl,:) = RhklDummyVec
              IF(IZolzFLAG.EQ.1) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                    RHKL(INhkl,:)= RhklDummyVec
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                 INhkl =  INhkl + 1
                 RHKL(INhkl,:) = RhklDummyVec                 
              END IF
           CASE DEFAULT
              PRINT*,"HKLMake(): unknown space group", SSpaceGroupName, "--- aborting"
              IErr=1
              RETURN
           END SELECT

        END DO
     END DO
  END DO

  
  RHKL(INhkl+1,:)= (/ 0.0D0,0.0D0,0.0D0 /)

END SUBROUTINE NewHKLmake
