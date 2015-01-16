!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all right reserved
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: diffractionpatterndefinitions.f90,v 1.11 2014/03/25 15:37:30 phsht Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ReflectionDetermination( IErr )

  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE SPara
  USE CPara
  USE IChannels
  USE BlochPara
  
  USE MyMPI

  IMPLICIT NONE

  REAL(RKIND) :: &
       dummy
  !norm,dummyVec(THREEDIM) -not used?
  INTEGER(IKIND) IErr, ind,jnd,icheck,ihklrun,IFind,IFound,knd
  
  CALL Message("ReflectionDetermination",IMust,IErr)
!!$  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
!!$     PRINT*,"DiffractionPatternDefinitions()"
!!$  END IF

  icheck = 0
  
  IHKLMAXValue = 15
  ihklrun = 0

  DO WHILE (icheck.EQ.0)
     ihklrun = ihklrun+1
     

     CALL Message("ReflectionDetermination",IInfo,IErr,MessageVariable = "IHklrun", &
          IVariable = IHklrun)


     
     CALL NewHKLMake(IHKLMAXValue,RZDirC,TWOPI/180.0D0,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Reflectiondetermination(", my_rank, ") error in NewHKLMake()"
        RETURN
     ENDIF
     
     IF(SIZE(RHKl,DIM=1).LT.IMinReflectionPool) THEN
        IHKLMAXValue = IHKLMAXValue*2
        Deallocate(RHKL,STAT=ierr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"ReflectionDetermination(", my_rank, ") error ", IErr, &
                " in DEALLOCATE() of DYNAMIC variables RHKL"
           RETURN
        ENDIF
        
        CYCLE
        
     ELSE
        icheck = 1
     END IF
     
     CALL ReSortHKL( RHKL, SIZE(RHKL,1),IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReflectionDetermination(): error in ReSortHKL()"
        RETURN
     ENDIF    
     
     IF(IXDirectionFLAG.EQ.0) THEN
        IDiffractionFLAG = 1
        RXDirC(1) = RHKL(2,1)
        RXDirC(2) = RHKL(2,2)
        RXDirC(3) = RHKL(2,3)
        CALL CrystallographyInitialisation( IErr )
        IF( IErr.NE.0 ) THEN
           PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error",IErr, &
                "in CrystallographyInitialisation()"
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
     
!!$        PRINT*,RHKL(:,1)
     ! G vector magnitudes in 1/Angstrom units
     
     DO ind=1,SIZE(RHKL,DIM=1)
        RgVecMag(ind)= SQRT(DOT_PRODUCT(RgVecMatT(ind,:),RgVecMatT(ind,:)))
!!$        PRINT*,RgVecMatT(ind,:)
     ENDDO
  
     RBSMaxGVecAmp = RgVecMag(IMinReflectionPool)
     
     nReflections = 0
     nStrongBeams = 0
     nWeakBeams = 0

!!$     PRINT*,SIZE(RHKL,DIM=1),SIZE(RgVecMag,DIM=1),SIZE(RgVecMatT,DIM=1)   
!!$     PRINT*,RBSMaxGVecAmp
     DO ind=1,SIZE(RHKL,DIM=1)
        IF (ABS(RgVecMag(ind)).LE.RBSMaxGVecAmp) THEN
           nReflections = nReflections + 1
        ENDIF
     ENDDO

!!$     PRINT*,"nReflections =",nReflections
     
  END DO

  END SUBROUTINE ReflectionDetermination
  
  SUBROUTINE SpecificReflectionDetermination (IErr)
    
    USE MyNumbers
    USE WriteToScreen

    USE IPara; USE RPara

    USE MyMPI

    IMPLICIT NONE

    INTEGER(IKIND) :: &
         IFind,IFound,ind,jnd,knd,IErr
    
  IFind = 0
  
  IF(IHKLSelectFLAG.EQ.1) THEN
     DO ind = 1,IReflectOut
        DO jnd = 1,SIZE(RHKL,DIM=1)
           IF(ABS(RHKL(jnd,1)-RInputHKLs(ind,1)).LE.RTolerance.AND.&
                ABS(RHKL(jnd,2)-RInputHKLs(ind,2)).LE.RTolerance.AND.&
                ABS(RHKL(jnd,3)-RInputHKLs(ind,3)).LE.RTolerance) THEN
              IFound = 0
              DO knd = 1,IFind
                 IF(ABS(RHKL(IOutputReflections(knd),1)-RInputHKLs(ind,1)).LE.RTolerance.AND.&
                      ABS(RHKL(IOutputReflections(knd),2)-RInputHKLs(ind,2)).LE.RTolerance.AND.&
                      ABS(RHKL(IOutputReflections(knd),3)-RInputHKLs(ind,3)).LE.RTolerance) THEN
                    IFound = 1
                    EXIT
                 END IF
              END DO

              IF(IFound.EQ.0) THEN
                 IFind = IFind +1
                 IOutputReflections(IFind) = jnd
                 
                 CALL Message("SpecificReflectionDetermination",IMoreInfo,IErr, &
                      MessageVariable = "Found HKL",RVariable = RInputHKLs(ind,jnd))
                 CALL Message("SpecificReflectionDetermination",IMoreInfo,IErr, &
                      MessageVariable = "At",IVariable = jnd)
                 
!!$                 IF((IWriteFLAG.GE.3.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
!!$                    PRINT*,"DiffractionPatternDefinitions(",my_rank,&
!!$                         ") Found HKL ",RInputHKLs(ind,:)," at ",jnd
!!$                 END IF
              ELSE 

                 CALL Message("SpecificReflectionDetermination",IMoreInfo,IErr, &
                      MessageVariable = "Found HKL",RVariable = RInputHKLs(ind,jnd))
                 CALL Message("SpecificReflectionDetermination",IMoreInfo,IErr, &
                      MessageVariable = "At",IVariable = jnd)
                 CALL Message("SpecificReflectionDetermination",IMoreInfo,IErr, &
                      MessageString = "However it is not unique")

!!$                 IF((IWriteFLAG.GE.3.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
!!$                    PRINT*,"DiffractionPatternDefinitions(",my_rank,&
!!$                         ") Found HKL ",RInputHKLs(ind,:)," at ",jnd,&
!!$                         "However it is not unique"
!!$                 END IF
                 
              END IF
              EXIT
           ELSE
              IF((jnd.EQ.SIZE(RHKL,DIM=1).AND.IWriteFLAG.GE.3.AND.my_rank.EQ.0).or.&
                   (jnd.EQ.SIZE(RHKL,DIM=1).AND.IWriteFLAG.GE.10)) THEN
                 PRINT*,"DiffractionPatternDefinitions(",my_rank,&
                      ") Could Not Find Requested HKL ",&
                      RInputHKLs(ind,:)," Will Ignore and Continue"
              END IF
              CYCLE
           END IF
        END DO
     END DO
     
     IF(IFind.LE.0) THEN
        IErr = 1
        IF( IErr.NE.0 ) THEN
           PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error ", IErr, &
                " No requested HKLs are allowed using the purposed geometry"
           RETURN
        ENDIF
     END IF
     IF(IReflectOut.NE.IFind) THEN
        IReflectOut = IFind
     END IF
  END IF
  
  IF((IWriteFLAG.GE.4.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     DO ind=1,SIZE(RHKL,DIM=1)
        PRINT*,RHKL(ind,:)
     END DO
  END IF

  END SUBROUTINE SpecificReflectionDetermination
 
  SUBROUTINE DiffractionPatternCalculation (IErr)

    USE MyNumbers
    USE WriteToScreen

    USE IPara; USE RPara;

    USE MyMPI
    
    IMPLICIT NONE

    INTEGER(IKIND) ind,IErr

    REAL(RKIND):: dummy

    CALL Message("DiffractionPatternDefinitions",IMust,IErr)

    ALLOCATE(&
         RGn(SIZE(RHKL,DIM=1)), &
         STAT=IErr)
    IF( IErr.NE.0 ) THEN
       PRINT*,"DiffractionPatternCalculation(", my_rank, ") error ", IErr, &
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
     PRINT*,"DiffractionPatternCalculation(", my_rank, ") error ", IErr, &
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

  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "BraggCentral", RVariable = RBraggCentral)

  
  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "GMax", RVariable = RBSMaxGVecAmp)     

  ! smallest g is gmag(2) IF 000 beam is included !!!add error catch here
  
  IF (ICentralBeamFLAG.EQ.1) THEN
     RMinimumGMag = RgVecMag(2)
  ELSE
     RMinimumGMag = RgVecMag(1)
  ENDIF

  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "MinimumGMag", RVariable = RMinimumGMag)

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

  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "No. of Reflections", IVariable = nReflections)

  ! resolution in k space
  RDeltaK = RMinimumGMag*RConvergenceAngle/IPixelCount

  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "RDeltaK", RVariable = RDeltaK)

  RETURN

END SUBROUTINE DiffractionPatternCalculation


SUBROUTINE NewHKLmake(Ihklmax,Rhkl0Vec,RAcceptanceAngle,IErr)
  
  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  USE MPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) IErr, Ihklmax,ind,jnd,knd,INhkl
  REAL(RKIND) RAcceptanceAngle
  REAL(RKIND), DIMENSION(THREEDIM) :: Rhkl0Vec,RhklDummyUnitVec,RhklDummyVec,Rhkl0UnitVec

  CALL Message("NewHKLMake",IMust,IErr)

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
