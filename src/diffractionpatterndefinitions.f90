!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-17, all rights reserved
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: diffractionpatterndefinitions.f90,v 1.11 2014/03/25 15:37:30 phsht Exp $
! $Id: diffractionpatterndefinitions.f90,v 2 2016/02/12 R.Beanland
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! for future removal
SUBROUTINE ReflectionDetermination( IErr )
!this is redundant for felixrefine
  USE MyNumbers
  USE WriteToScreen
  USE ErrorCodes

  USE CConst; USE IConst
  USE IPara; USE RPara
  USE SPara
  USE CPara
  USE IChannels
  USE BlochPara

  USE MyMPI

  IMPLICIT NONE
  INTEGER(IKIND) :: IErr
  PRINT*,"Attempt to use subroutine ReflectionDetermination!"
  IErr=1

END SUBROUTINE ReflectionDetermination

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
SUBROUTINE SpecificReflectionDetermination (IErr)
!What does this do    
  USE MyNumbers
  USE WriteToScreen
  USE IConst
  USE IPara; USE RPara
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IFind,IFound,ind,jnd,knd,IErr
    
  IFind = 0

  IF(IHKLSelectFLAG.EQ.1) THEN
    DO ind = 1,SIZE(RInputHKLs,DIM=1)
      DO jnd = 1,SIZE(Rhkl,DIM=1)
        IF(ABS(Rhkl(jnd,1)-RInputHKLs(ind,1)).LE.RTolerance.AND.&
           ABS(Rhkl(jnd,2)-RInputHKLs(ind,2)).LE.RTolerance.AND.&
           ABS(Rhkl(jnd,3)-RInputHKLs(ind,3)).LE.RTolerance) THEN
          IFound = 0
          DO knd = 1,IFind
            IF(ABS(Rhkl(IOutputReflections(knd),1)-RInputHKLs(ind,1)).LE.RTolerance.AND.&
               ABS(Rhkl(IOutputReflections(knd),2)-RInputHKLs(ind,2)).LE.RTolerance.AND.&
               ABS(Rhkl(IOutputReflections(knd),3)-RInputHKLs(ind,3)).LE.RTolerance) THEN
              IFound = 1
              EXIT
            END IF
          END DO
          IF(IFound.EQ.0) THEN
            IFind = IFind +1
            IOutputReflections(IFind) = jnd
          END IF
          EXIT
        ELSE
          IF((jnd.EQ.SIZE(Rhkl,DIM=1).AND.IWriteFLAG.GE.3.AND.my_rank.EQ.0).or.&
                   (jnd.EQ.SIZE(Rhkl,DIM=1).AND.IWriteFLAG.GE.10)) THEN
            PRINT*,"DiffractionPatternDefinitions(",my_rank,&
                      ") Could Not Find Requested HKL ",&
                      NINT(RInputHKLs(ind,:))," Will Ignore and Continue"
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
        END IF
     END IF
     IF(INoOfLacbedPatterns.NE.IFind) THEN
        INoOfLacbedPatterns = IFind
     END IF
  END IF
  
END SUBROUTINE SpecificReflectionDetermination

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DiffractionPatternCalculation (IErr)
  !RB This is now redundant
  USE MyNumbers
  USE WriteToScreen
  USE IConst
  
  USE IPara; USE RPara;
  
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,IErr
  CHARACTER*20 :: Sind
  
  DO ind =1,SIZE(Rhkl,DIM=1)
     RgDotNorm(ind) = DOT_PRODUCT(RgPool(ind,:),RNormDirM)
  END DO   
  
  ! smallest g is gmag(2) IF 000 beam is included !!!add error catch here
  RMinimumGMag = RgPoolMag(2)
  
  IF (nReflections.LT.INoOfLacbedPatterns) THEN
     nReflections = INoOfLacbedPatterns
  END IF
  
  ! resolution in k space
  RDeltaK = RMinimumGMag*RConvergenceAngle/REAL(IPixelCount,RKIND)

  RETURN

END SUBROUTINE DiffractionPatternCalculation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE HKLCount(Ihklmax,Rhkl0Vec,INhkl,RHOLZAcceptanceAngle,IErr)
!Counts the number of reflections limited by Ihklmax and the acceptance angle, returns it as INhkl

  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  USE MPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,Ihklmax,ind,jnd,knd,INhkl
  REAL(RKIND) :: RHOLZAcceptanceAngle
  REAL(RKIND), DIMENSION(ITHREE) :: Rhkl0Vec,RhklDummyUnitVec,RhklDummyVec,Rhkl0UnitVec

  INhkl = 0
  Rhkl0UnitVec= Rhkl0Vec/SQRT(DOT_PRODUCT(REAL(Rhkl0Vec,RKIND),REAL(Rhkl0Vec,RKIND)))

  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1          
           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)          
           IF( (ind.NE.0).OR.(jnd.NE.0).OR.(knd.NE.0) ) THEN
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
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1
                 END IF
              END IF
			  
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("P")! Primitive
              IF(IHolzFLAG.EQ.0) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                   .LE.SIN(RHOLZAcceptanceAngle)) THEN
                 INhkl = INhkl +1       
              END IF
			  
           CASE DEFAULT
              PRINT*,"HKLMake: unknown space group", SSpaceGroupName, ", aborting"
              IErr=1
              RETURN
           END SELECT
		   
        END DO
     END DO
  END DO
END SUBROUTINE HKLCount

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE HKLMake(Ihklmax,Rhkl0Vec,RHOLZAcceptanceAngle,IErr)
!Fills the list of reciprocal space vectors Rhkl 
!N.B. Rhkl is simply real versions of h,k,l 

  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  USE MPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,Ihklmax,ind,jnd,knd,lnd
  REAL(RKIND) :: RHOLZAcceptanceAngle
  REAL(RKIND),DIMENSION(ITHREE) :: Rhkl0Vec,RhklDummyUnitVec,RhklDummyVec,Rhkl0UnitVec

  lnd = 0
  Rhkl0UnitVec= Rhkl0Vec/SQRT(DOT_PRODUCT(REAL(Rhkl0Vec,RKIND),REAL(Rhkl0Vec,RKIND)))!?isn't RhklDummyVec already real??

  DO ind=-Ihklmax,Ihklmax,1
    DO jnd=-Ihklmax,Ihklmax,1
      DO knd=-Ihklmax,Ihklmax,1
        RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)
        IF(ind.NE.0.AND.jnd.NE.0.AND.knd.NE.0) THEN
          RhklDummyUnitVec= RhklDummyVec / &
            SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))!?isn't RhklDummyVec already real??
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
              IF(IHolzFLAG.EQ.0) THEN
                IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                  lnd=lnd+1
                  Rhkl(lnd,:)=RhklDummyVec
                END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                lnd=lnd+1
                Rhkl(lnd,:)=RhklDummyVec                 
              END IF
            END IF
			  
          CASE("I")! Body Centred
            IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
              IF(IHolzFLAG.EQ.0) THEN
                IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                  lnd=lnd+1
                  Rhkl(lnd,:)= RhklDummyVec
                END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                lnd =  lnd + 1
                Rhkl(lnd,:) = RhklDummyVec                 
              END IF
            END IF
			  
          CASE("A")! A-Face Centred
            IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
              IF(IHolzFLAG.EQ.0) THEN
                IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                  lnd=lnd+1
                  Rhkl(lnd,:)= RhklDummyVec
                END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                lnd =  lnd + 1
                Rhkl(lnd,:) = RhklDummyVec                 
              END IF
            END IF
			  
           CASE("B")! B-Face Centred
             IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
               IF(IHolzFLAG.EQ.0) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                   lnd=lnd+1
                   Rhkl(lnd,:)= RhklDummyVec
                 END IF
               ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                 lnd =  lnd + 1
                 Rhkl(lnd,:) = RhklDummyVec                 
               END IF
             END IF
			  
           CASE("C")! C-Face Centred
             IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY) THEN
               IF(IHolzFLAG.EQ.0) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                   lnd=lnd+1
                   Rhkl(lnd,:)= RhklDummyVec
                 END IF
               ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                 lnd =  lnd + 1
                 Rhkl(lnd,:) = RhklDummyVec                 
               END IF
             END IF
			  
           CASE("R")! Rhombohedral Reverse
             IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
               IF(IHolzFLAG.EQ.0) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                   lnd=lnd+1
                   Rhkl(lnd,:)= RhklDummyVec
                 END IF
               ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                 lnd =  lnd + 1
                 Rhkl(lnd,:) = RhklDummyVec                 
               END IF
             END IF
			  
           CASE("V")! Rhombohedral Obverse
             IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
               IF(IHolzFLAG.EQ.0) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                   lnd=lnd+1
                   Rhkl(lnd,:)= RhklDummyVec
                 END IF
               ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                 lnd =  lnd + 1
                 Rhkl(lnd,:) = RhklDummyVec                 
               END IF
             END IF
			  
           CASE("P")! Primitive
             IF(IHolzFLAG.EQ.0) THEN
               IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                 lnd=lnd+1
                 Rhkl(lnd,:)= RhklDummyVec
               END IF
             ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
               lnd =  lnd + 1
               Rhkl(lnd,:) = RhklDummyVec                 
             END IF
			  
           CASE DEFAULT!RB should never get here since already been through HKLcount
             PRINT*,"HKLMake(): unknown space group", SSpaceGroupName, "--- aborting"
             IErr=1
             RETURN
         END SELECT

        END DO
     END DO
  END DO

END SUBROUTINE HKLmake

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE NewHKLMake(Ihklmax,Rhkl0Vec,RHOLZAcceptanceAngle,IErr)
  
  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  USE MPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr, Ihklmax,ind,jnd,knd,INhkl
  REAL(RKIND) :: RHOLZAcceptanceAngle
  REAL(RKIND), DIMENSION(ITHREE) :: Rhkl0Vec,RhklDummyUnitVec,RhklDummyVec,Rhkl0UnitVec

  INhkl = 0
  
  Rhkl0UnitVec= Rhkl0Vec/SQRT(DOT_PRODUCT(REAL(Rhkl0Vec,RKIND),REAL(Rhkl0Vec,RKIND)))

!First count the number of reflections in the acceptance angle
!??? doesn't this only work for cubic systems?  Where are the magnitudes of the reciprocal lattice vectors?
  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1          
           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)          
           IF( (ind.NE.0).OR.(jnd.NE.0).OR.(knd.NE.0) ) THEN
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
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1
                 END IF
              END IF
			  
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 END IF
              END IF
			  
           CASE("P")! Primitive
              IF(IHolzFLAG.EQ.0) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                   .LE.SIN(RHOLZAcceptanceAngle)) THEN
                 INhkl = INhkl +1       
              END IF
           CASE DEFAULT
              PRINT*,"HKLMake: unknown space group", SSpaceGroupName, ", aborting"
              IErr=1
              RETURN
           END SELECT

        END DO
     END DO
  END DO

!RB now allocate the hkl list...  
  ALLOCATE(Rhkl(INhkl,ITHREE),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"hklMake(",my_rank,")error allocating Rhkl"
     RETURN
  END IF

!RB ...and calculate it all again, filling Rhkl  
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
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                IF(IHolzFLAG.EQ.0) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("P")! Primitive

		   IF(IHolzFLAG.EQ.0) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                    Rhkl(INhkl,:)= RhklDummyVec
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                 INhkl =  INhkl + 1
                 Rhkl(INhkl,:) = RhklDummyVec                 
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
