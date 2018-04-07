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

!>
!! Module-description: 
!!
MODULE setup_reflections_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: HKLList, HKLCount, HKLMake, HKLSort

  CONTAINS

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Assign numbers to the different reflections in IOutputReflections
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!  
  SUBROUTINE HKLList( IErr )

    ! This procedure is called once in felixrefine setup
    USE MyNumbers
    USE message_mod

    ! global outputs (or inout)
    USE IPARA, ONLY : IOutputReflections, INoOfLacbedPatterns

    ! global inputs
    USE IPARA, ONLY : IHKLSelectFLAG
    USE RPARA, ONLY : RInputHKLs, RTolerance, Rhkl

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
            IF (IFound.EQ.0) THEN
              IFind = IFind +1
              IOutputReflections(IFind) = jnd
            END IF
            EXIT
          ELSE
            IF ( jnd.EQ.SIZE(Rhkl,DIM=1) ) THEN
              CALL message(LM,"No requested HKLs found",NINT(RInputHKLs(ind,:)) )
              CALL message(LM,"Will Ignore and Continue")
            END IF
            CYCLE
          END IF

        END DO
      END DO
       
      IF (IFind.LE.0) THEN
        IErr=1; IF(l_alert(IErr,"HKLList","No requested HKLs found")) RETURN
      END IF
      
      INoOfLacbedPatterns = IFind

    END IF
    
  END SUBROUTINE HKLList
  
  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Counts the number of reflections limited by Ihklmax
  !! and the acceptance angle, returns it as INhkl
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!  
  SUBROUTINE HKLCount(Ihklmax, Rhkl0Vec, INhkl, RHOLZAcceptanceAngle, IErr)

    ! This procedure is called once in felixrefine setup
    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE IPARA, ONLY : IHolzFLAG
    USE RPARA, ONLY : Rhkl
    USE SPARA, ONLY : SSpaceGroupName
    
    IMPLICIT NONE
    INTEGER(IKIND),INTENT(OUT) :: INhkl, IErr
    INTEGER(IKIND),INTENT(IN) :: Ihklmax
    REAL(RKIND),INTENT(IN) :: Rhkl0Vec(ITHREE), RHOLZAcceptanceAngle
    INTEGER(IKIND) :: ind, jnd, knd
    REAL(RKIND),DIMENSION(ITHREE) :: RhklDummyUnitVec, RhklDummyVec, Rhkl0UnitVec

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
                IF( ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),THREE)) &
                .LE.TINY) THEN
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
                IErr=1; IF(l_alert(IErr,"HKLCount",&
                      "SSpaceGroupName unrecognised")) RETURN
             END SELECT
		     
          END DO
       END DO
    END DO

  END SUBROUTINE HKLCount

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Fills the list of reciprocal space vectors Rhkl
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!  
  SUBROUTINE HKLMake(Ihklmax, Rhkl0Vec, RHOLZAcceptanceAngle, IErr)   
    ! This procedure is called once in felixrefine setup

    ! N.B. Rhkl is simply real versions of h,k,l 

    USE MyNumbers
    USE message_mod
    
    ! global outputs (or inout)
    USE RPARA, ONLY : Rhkl
      
    ! global inputs
    USE SPARA, ONLY : SSpaceGroupName
    USE IPARA, ONLY : IHolzFLAG
    USE Iconst
    
    IMPLICIT NONE
    
    INTEGER(IKIND) :: IErr, Ihklmax, ind, jnd, knd, lnd
    REAL(RKIND) :: RHOLZAcceptanceAngle
    REAL(RKIND),DIMENSION(ITHREE) :: Rhkl0Vec, RhklDummyUnitVec, RhklDummyVec, Rhkl0UnitVec

    lnd = 0
    Rhkl0UnitVec = Rhkl0Vec / SQRT(DOT_PRODUCT(Rhkl0Vec,Rhkl0Vec))

    DO ind=-Ihklmax,Ihklmax,1
      DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1
          RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)
          IF(ind.NE.0.AND.jnd.NE.0.AND.knd.NE.0) THEN
            RhklDummyUnitVec = RhklDummyVec / &
                  SQRT(DOT_PRODUCT(RhklDummyVec,RhklDummyVec))
          ELSE
            RhklDummyUnitVec = RhklDummyVec
          END IF
             
          SELECT CASE(SSpaceGroupName)
            CASE("F") ! Face Centred
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
                ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.sin(RHOLZAcceptanceAngle)) THEN
                  lnd=lnd+1
                  Rhkl(lnd,:)=RhklDummyVec                 
                END IF
              END IF
			    
            CASE("I") ! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                IF(IHolzFLAG.EQ.0) THEN
                  IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    lnd=lnd+1
                    Rhkl(lnd,:)= RhklDummyVec
                  END IF
                ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.sin(RHOLZAcceptanceAngle)) THEN
                  lnd =  lnd + 1
                  Rhkl(lnd,:) = RhklDummyVec                 
                END IF
              END IF
			    
            CASE("A") ! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                IF(IHolzFLAG.EQ.0) THEN
                  IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    lnd=lnd+1
                    Rhkl(lnd,:)= RhklDummyVec
                  END IF
                ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.sin(RHOLZAcceptanceAngle)) THEN
                  lnd =  lnd + 1
                  Rhkl(lnd,:) = RhklDummyVec                 
                END IF
              END IF
			    
            CASE("B") ! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                IF(IHolzFLAG.EQ.0) THEN
                  IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    lnd=lnd+1
                    Rhkl(lnd,:)= RhklDummyVec
                  END IF
                ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.sin(RHOLZAcceptanceAngle)) THEN
                  lnd =  lnd + 1
                  Rhkl(lnd,:) = RhklDummyVec                 
                END IF
              END IF
			    
            CASE("C") ! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY) THEN
                IF(IHolzFLAG.EQ.0) THEN
                  IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    lnd=lnd+1
                    Rhkl(lnd,:)= RhklDummyVec
                  END IF
                ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.sin(RHOLZAcceptanceAngle)) THEN
                  lnd =  lnd + 1
                  Rhkl(lnd,:) = RhklDummyVec                 
                END IF
              END IF
			    
             CASE("R") ! Rhombohedral Reverse
               IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),THREE)) &
                    .LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                   IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                     lnd=lnd+1
                     Rhkl(lnd,:)= RhklDummyVec
                   END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.sin(RHOLZAcceptanceAngle)) THEN
                   lnd =  lnd + 1
                   Rhkl(lnd,:) = RhklDummyVec                 
                 END IF
               END IF
			    
             CASE("V") ! Rhombohedral Obverse
               IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),THREE)) &
                    .LE.TINY) THEN
                 IF(IHolzFLAG.EQ.0) THEN
                   IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                     lnd=lnd+1
                     Rhkl(lnd,:)= RhklDummyVec
                   END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.sin(RHOLZAcceptanceAngle)) THEN
                   lnd =  lnd + 1
                   Rhkl(lnd,:) = RhklDummyVec                 
                 END IF
               END IF
			    
             CASE("P") ! Primitive
               IF(IHolzFLAG.EQ.0) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                   lnd=lnd+1
                   Rhkl(lnd,:)= RhklDummyVec
                 END IF
               ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                    .LE.sin(RHOLZAcceptanceAngle)) THEN
                 lnd =  lnd + 1
                 Rhkl(lnd,:) = RhklDummyVec                 
               END IF
			    
             CASE DEFAULT ! RB should never get here since already been through HKLcount
               IErr=1; IF(l_alert(IErr,"HKLmake",&
                     "SSpaceGroupName unrecognised")) RETURN

           END SELECT

          END DO
       END DO
    END DO

  END SUBROUTINE HKLmake

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Sorts array into descending order
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE HKLSort(LocalRhkl,N,IErr )

    !--------------------------------------------------------------------
    !	Sort: is based on ShellSort from "Numerical Recipes", routine SHELL().
    !---------------------------------------------------------------------  

    USE MyNumbers
      
    USE SConst; USE IConst
    USE IPara; USE RPara

    USE IChannels

    USE MPI
    USE MyMPI

    IMPLICIT NONE

    INTEGER (IKIND) :: IErr,NN,M,L,K,J,I,LOGNB2,ind
    INTEGER (IKIND),INTENT(IN) :: N
    REAL(RKIND),INTENT(INOUT) :: LocalRhkl(N,ITHREE)
    REAL(RKIND) :: RhklSearch(ITHREE), RhklCompare(ITHREE)
    REAL(RKIND) :: ALN2I,LocalTINY,dummy
    PARAMETER (ALN2I=1.4426950D0, LocalTINY=1.D-5)

    LOGNB2=INT(LOG(REAL(N))*ALN2I+LocalTINY)
    M=N
    DO NN=1,LOGNB2
      M=M/2
      K=N-M
      DO J=1,K
        I=J
3       CONTINUE
        L=I+M
        RhklSearch = LocalRhkl(L,1)*RarVecO + &
           LocalRhkl(L,2)*RbrVecO+LocalRhkl(L,3)*RcrVecO    
        RhklCompare = LocalRhkl(I,1)*RarVecO + &
           LocalRhkl(I,2)*RbrVecO+LocalRhkl(I,3)*RcrVecO
        IF( DOT_PRODUCT(RhklSearch(:),RhklSearch(:)) .LT. &
              DOT_PRODUCT(RhklCompare(:),RhklCompare(:))) THEN
          DO ind=1,ITHREE
            dummy=LocalRhkl(I,ind)
            LocalRhkl(I,ind)=LocalRhkl(L,ind)
            LocalRhkl(L,ind)=dummy
          ENDDO
          I=I-M
          IF(I.GE.1) GOTO 3
        ENDIF
      ENDDO
    ENDDO
    
    RETURN

  END SUBROUTINE HKLSort
  
END MODULE setup_reflections_mod
