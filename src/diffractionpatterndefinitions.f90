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
       dummy,RMaxAcceptanceGVecMag
  INTEGER(IKIND) :: &
       IErr, ind,jnd,icheck,ihklrun,IFind,IFound,knd,IMaxLaueZoneLevel
  CHARACTER*20 :: &
       Sind
  
  CALL Message("ReflectionDetermination",IMust,IErr)

  icheck = 0
  
  IHKLMAXValue = 15
  ihklrun = 0

  DO WHILE (icheck.EQ.0)
     ihklrun = ihklrun+1     

     CALL Message("ReflectionDetermination",IInfo,IErr,MessageVariable = "IHklrn", &
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
  END DO
  
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
  
  ALLOCATE(&
       RgVecMagLaueZone(SIZE(RHKL,DIM=1),2), & !TWO dimension here from smodules to add
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
     WRITE(Sind,'(I10.1)')ind
     CALL Message("ReflectionDetermination", IMoreInfo,IErr, &
          MessageVariable="Reciprocal Vector" &
          //",RgVecMatT"//"("//ADJUSTL(TRIM(Sind))//":)", &
          RVariable=RgVecMatT(ind,3))
  ENDDO
  !produce vector from Acceptance angle - mrad to 1/Angstrom max,
!!$  any values outside this can be cut...
  !have to think about boundaries, if equal to or less than, then allow.
  !restrict MinReflectionPool on this basis, ie. change MinReflectionPool 
  ! G vector magnitudes in 1/Angstrom units

  IMaxLaueZoneLevel=INT(MAXLOC(RgVecMatT(:,3),DIM=1))

!!$     DO ind=1,SIZE(RHKL,DIM=1)
!!$        DO jnd=0,IMaxLaueZoneLevel
!!$           WHERE(RgVecMatT(ind,3).EQ.RMaxLaueZoneLevel)
!!$              RgVecMagLaueZone(ind,2)= &
!!$                   SQRT(DOT_PRODUCT(RgVecMatT(ind,:),RgVecMatT(ind,:)))
!!$           END WHERE
!!$        END DO
!!$     END DO
         
  
  DO ind=1,SIZE(RHKL,DIM=1)
     RgVecMag(ind)= SQRT(DOT_PRODUCT(RgVecMatT(ind,:),RgVecMatT(ind,:)))
  ENDDO

!!$If the user specifies an acceptance angle, felix uses it to restrict the number
!!$of reflections in IMinReflectionPool & nReflections
!!$check to ensure everything is only in magnitudes (there are negative g-vectors) 

  IF(RAcceptanceAngle.NE.ZERO.AND.IZOLZFLAG.EQ.1) THEN
     CALL Message("ReflectionDetermination",IInfo,IErr, &
          MessageString="Determining restriction effect of Acceptance Angle" &
          // " (in k-space for ZOLZ only). If unintended, please cancel the simulation and" &
          // " set the RAcceptance angle to 0.0, and/or switch IZOLZFLAG to 0")
     RMaxAcceptanceGVecMag=(RElectronWaveVectorMagnitude*TAN(RAcceptanceAngle*1E-3))

     CALL Message("ReflectionDetermination",IInfo,IErr, &
          MessageVariable="RMaxGVecMAG",RVariable=RMaxAcceptanceGVecMag)

     CALL Message("ReflectionDetermination",IInfo,IErr, &
          MessageVariable="RElectronWaveVector",RVariable=RElectronWaveVectorMagnitude)

     IF(RgVecMag(IMinReflectionPool).GT.RMaxAcceptanceGVecMag) THEN
        CALL Message("ReflectionDetermination",IInfo,IErr, &
             MessageString="Number of Reflections (IMinReflectionPool) Exceeds cut-off from" &
             //" Acceptance angle, calculating new cut-off value (reciprocal angstroms)")
!!$     New max Gvector amplitude is the G-vector specified by the acceptance angle
        RBSMaxGVecAmp = RMaxAcceptanceGVecMag 
     ELSE
        CALL Message("ReflectionDetermination",IInfo,IErr, &
             MessageString="Warning: Number of Reflections in Reflection pool does not" &
             //" exceed the Acceptance angle (reciprocal angstroms)" &
             //" continuing in normal mode, for full range increase IMinReflectionPool")
!!$     Normal Reflection Pool value
        RBSMaxGVecAmp = RgVecMag(IMinReflectionPool)
     END IF
     
  ELSEIF(RAcceptanceAngle.NE.ZERO.AND.IZOLZFLAG.EQ.0) THEN
     CALL Message("ReflectionDetermination",IInfo,IErr, &
          MessageString="Determining restriction effect of Acceptance Angle" &
          // " (in k-space with HOLZ). If unintended, please cancel the simulation and" &
          // " set the RAcceptance angle to 0.0 and/or switch IZOLZFLAG to 1")
  ELSE
     RBSMaxGVecAmp = RgVecMag(IMinReflectionPool)
  END IF

  
  nReflections = 0
  nStrongBeams = 0
  nWeakBeams = 0

  

  DO ind=1,SIZE(RHKL,DIM=1)
     IF (ABS(RgVecMag(ind)).LE.RBSMaxGVecAmp) THEN
        WRITE(Sind,'(I10.1)')ind
        CALL Message("ReflectionDetermination", IMoreInfo,IErr, &
             MessageVariable="Allowed Reciprocal g-vector magnitude" &
             //",RgVecMag" // "("//ADJUSTL(TRIM(Sind))//")", &
             RVariable=RgVecMag(ind))
        nReflections = nReflections + 1
     ENDIF
  ENDDO
  
  CALL Message("ReflectionDetermination", IInfo,IErr, &
       MessageVariable="Reflection Pool reduced from, IMinReflectionpool", &
       IVariable=IMinReflectionPool)
  CALL Message("ReflectionDetermination", IInfo,IErr, &
       MessageVariable="to nReflections",IVariable=nReflections)
  
  END SUBROUTINE ReflectionDetermination
  
  SUBROUTINE SpecificReflectionDetermination (IErr)
    
    USE MyNumbers
    USE WriteToScreen
    USE IConst

    USE IPara; USE RPara

    USE MyMPI

    IMPLICIT NONE

    INTEGER(IKIND) :: &
         IFind,IFound,ind,jnd,knd,IErr

    CALL Message("SpecificReflectionDetermination",IMust,IErr)
    
  IFind = 0
  
  IF(IHKLSelectFLAG.EQ.1) THEN
     DO ind = 1,SIZE(RInputHKLs,DIM=1)
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
                      MessageVariable = "At",IVariable = jnd)
              ELSE 
                 CALL Message("SpecificReflectionDetermination",IMoreInfo,IErr, &
                      MessageVariable = "At",IVariable = jnd)
                 CALL Message("SpecificReflectionDetermination",IMoreInfo,IErr, &
                      MessageString = "However it is not unique")
                 
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
  
END SUBROUTINE SpecificReflectionDetermination

SUBROUTINE DiffractionPatternCalculation (IErr)
  
  USE MyNumbers
  USE WriteToScreen
  USE IConst
  
  USE IPara; USE RPara;
  
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       ind,IErr
  
  CALL Message("DiffractionPatternCalculation",IMust,IErr)
  
  ALLOCATE(&
       RgVecVec(SIZE(RHKL,DIM=1)), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RgVecVec(HKL)"
     RETURN
  ENDIF
  
  RNormDirM = RNormDirM/sqrt(DOT_PRODUCT(RNormDirM,RNormDirM))
  
  DO ind =1,SIZE(RHKL,DIM=1)
     RgVecVec(ind) = DOT_PRODUCT(RgVecMatT(ind,:),RNormDirM)
  END DO
  
  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "GMax", RVariable = RBSMaxGVecAmp)     
  
  ! smallest g is gmag(2) IF 000 beam is included !!!add error catch here
  
  RMinimumGMag = RgVecMag(2)
  
  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "MinimumGMag", RVariable = RMinimumGMag)
  
  IF (nReflections.LT.IReflectOut) THEN
     IReflectOut = nReflections
  END IF
  
  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "No. of Reflections", IVariable = nReflections)
  
  ! resolution in k space
  RDeltaK = RMinimumGMag*RConvergenceAngle/REAL(IPixelCount,RKIND)

  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "RDeltaK", RVariable = RDeltaK)

  RETURN

END SUBROUTINE DiffractionPatternCalculation


SUBROUTINE NewHKLMake(Ihklmax,Rhkl0Vec,RHOLZAcceptanceAngle,IErr)
  
  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  USE MPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       IErr, Ihklmax,ind,jnd,knd,INhkl
  REAL(RKIND) :: &
       RHOLZAcceptanceAngle
  REAL(RKIND), DIMENSION(THREEDIM) :: &
       Rhkl0Vec,RhklDummyUnitVec,RhklDummyVec,Rhkl0UnitVec

  CALL Message("NewHKLMake",IMust,IErr)

  INhkl = 0

  Rhkl0UnitVec= Rhkl0Vec/SQRT(DOT_PRODUCT(REAL(Rhkl0Vec,RKIND),REAL(Rhkl0Vec,RKIND)))
  
  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1
           
           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)

           IF(ind.NE.0.AND.jnd.NE.0.AND.knd.NE.0) THEN
              RhklDummyUnitVec= RhklDummyVec / &
                   SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))
           ELSE
              RhklDummyUnitVec = RhklDummyVec
           END IF

           SELECT CASE(SSpaceGroupName)
           CASE("F") !Face Centred
              IF(((ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),TWO))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),TWO))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),TWO))).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),TWO))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),TWO))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),TWO))).GT.TINY))) THEN
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                    
                 ENDIF
                 ! INhkl = INhkl + 1
              END IF              
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                    
                 ENDIF
              END IF
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                    
                 ENDIF
              END IF
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                !INhkl = INhkl + 1
                
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
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
                   .LE.SIN(RHOLZAcceptanceAngle)) THEN
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
       RHKL((INhkl),THREEDIM),&
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

           IF(ind.NE.0.AND.jnd.NE.0.AND.knd.NE.0) THEN
              RhklDummyUnitVec= RhklDummyVec / &
                   SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))
           ELSE
              RhklDummyUnitVec = RhklDummyVec
           END IF
           
           SELECT CASE(SSpaceGroupName)
           CASE("F") !Face Centred
              IF(((ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),TWO))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),TWO))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),TWO))).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),TWO))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),TWO))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),TWO))).GT.TINY))) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
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
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
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

END SUBROUTINE NewHKLmake
