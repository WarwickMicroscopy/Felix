!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-17, all rights reserved
!
! Version: 2.0
! Date: 19-12-2022
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build: cRED
! Author:  r.beanland@warwick.ac.uk
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

!>
!! Module-description: 
!!
MODULE setup_reflections_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: HKLList, SelectionRules, HKLMake, HKLSort

  CONTAINS
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Fills the beam pool list RgPoolList for the current frame
  !! and the list of output reflections IgOutList (global variables)
  !!
  !! Major-Authors: Richard Beanland (2021)
  !!  
  SUBROUTINE HKLMake(IFrame, RGPoolLimit, RGOutLimit, IErr)   
    ! This setup procedure fills up RgPoolList and IgOutList for each frame

    USE MyNumbers
    USE message_mod
    
    ! global parameters Rhkl, input reciprocal lattice vectors & wave vector
    USE RPARA, ONLY : RzDirC, RgPoolList, RarVecM, RbrVecM, RcrVecM, RInputHKLs,&
        RElectronWaveVectorMagnitude, IgOutList
      
    ! global inputs
    USE SPARA, ONLY : SPrintString
    USE IPARA, ONLY : INhkl, INoOfHKLsAll
    USE Iconst
    
    IMPLICIT NONE

    INTEGER(IKIND),INTENT(IN) :: IFrame
    REAL(RKIND),INTENT(IN) :: RGPoolLimit, RGOutLimit
    INTEGER(IKIND) :: IErr, ISel, Ih, Ik, Il, inda,indb,indc, jnd, knd,lnd
    REAL(RKIND) :: RarMag, RbrMag, RcrMag, RShell, RGtestMag, RDev
    REAL(RKIND),DIMENSION(ITHREE) :: Rk, RGtest, RGtestM, RGplusk 
   
    ! RGPoolLimit is the upper limit for g-vector magnitudes
    ! If the g-vectors we are counting are bigger than this there is something wrong
    ! probably the tolerance for proximity to the Ewald sphere needs increasing
    ! could be an input in felix.inp

    ! Rk is the k-vector for the incident beam
    ! we are working in the microscope reference frame so k is along z
    Rk=(/ 0.0,0.0,RElectronWaveVectorMagnitude /)
    ! we work our way out from the origin in shells
    ! the shell increment is the smallest basis vector
    RShell=MINVAL( (/ RarMag,RbrMag,RcrMag /) )

    !first g is always 000
    RgPoolList(IFrame,1,:)=(/ 0.0,0.0,0.0 /)
    knd=1!number of reflections in the pool 
    lnd=0!number of the shell 
    !maximum a*,b*,c* limit is determined by the G magnitude limit
    inda=NINT(RGPoolLimit/RarMag)
    indb=NINT(RGPoolLimit/RbrMag)
    indc=NINT(RGPoolLimit/RcrMag)

    !fill the RgPoolList with beams near the Bragg condition
    DO WHILE (REAL(lnd)*RShell.LT.RGPoolLimit)
      !increment the shell
      lnd = lnd+1
!DBG  IF(my_rank.EQ.0)PRINT*,REAL(lnd-1)*RShell,"to",REAL(lnd)*RShell
      !Make a hkl
      DO Ih = -inda,inda
         DO Ik = -indb,indb
            DO Il = -indc,indc
              !check that it's allowed by selection rules
              ISel=0
              CALL SelectionRules(Ih, Ik, Il, ISel, IErr)
              !and check that we have space 
              IF (ISel.EQ.1 .AND. knd.LT.INhkl) THEN
                !Make a g-vector
                RGtest = REAL( (/ Ih,Ik,Il /),RKIND )!Miller indices
                RGtestM = REAL(Ih)*RarVecM+REAL(Ik)*RbrVecM+REAL(Il)*RcrVecM!in microscope frame
                RGtestMag = SQRT(DOT_PRODUCT(RGtestM,RGtestM))
                !is it in the shell
                IF (RGtestMag.GT.REAL(lnd-1)*RShell.AND.RGtestMag.LE.REAL(lnd)*RShell) THEN
                  !is it near a Laue condition |k+g|=|k|
                  RGplusk=RGtestM+Rk
                  !divide by |g| to get a measure independent of incident beam tilt
                  RDev=ABS(RElectronWaveVectorMagnitude-SQRT(DOT_PRODUCT(RGplusk,RGplusk)))/RGtestMag
                  ! Using global variable RDevLimit to include or not
!DBG              IF(my_rank.EQ.0)PRINT*,(/ Ih,Ik,Il /),RDev
                  IF ((RDev-RDevLimit).LT.TINY) THEN !it's near the ZOLZ 
                    !add it to the pool and increment the counter
                    knd=knd+1
                    RgPoolList(IFrame,knd,:)=RGtest
                    IF (RGtestMag.LT.RGOutLimit) THEN
                      IgOutList(IFrame,knd)=1
                    END IF
                  END IF
                END IF
              END IF
            END DO
         END DO
      END DO
    END DO

    ! it is possible that we reach RGPoolLimit before filling up the pool, in which case 
    ! fill up any remaining beam pool places with an enormous g-vector, diabolically
    ! the idea being that this g-vector will never be near any possible Ewald sphere
    IF (knd.LT.INhkl) THEN
      RGtest = REAL( (/ 666,666,666 /),RKIND )
      DO jnd = knd+1, INhkl
        Rhkl(jnd,:)=RGtest
      END DO
    END IF
    
    !output which required output hkl's are in this frame
    CALL message(LM, "Reflection list:")
    DO knd = 1, INhkl
      IF (lnd.EQ.1) THEN
        WRITE(SPrintString,'(I3,1X,I3,1X,I3)') NINT(RInputHKLs(jnd,:))
        CALL message(LM, SPrintString)
      END IF
    END DO
  END SUBROUTINE HKLmake

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Sorts Rhkl array into descending order
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
  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Assign numbers to the different reflections
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!  
  SUBROUTINE HKLList(IFrame, IErr )

    ! This procedure is called once in felixrefine setup
    USE MyNumbers
    USE message_mod

    ! global parameters
    USE IPARA, ONLY : IhklsFrame, INoOfHKLsAll, IhklsAll, INoOfHKLsFrame,&
            ILiveList
    USE SPARA, ONLY : SPrintString

    ! global inputs
    USE RPARA, ONLY : RInputHKLs, Rhkl

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(IN) :: IFrame
    INTEGER(IKIND) :: IFind,IDuplicate,ind,jnd,knd,IErr

    !--------------------------------------------------------------------
    ! IhklsFrame links the list in felix.hkl with the beam pool for the current frame
    ! It has length INoOfHKLsAll but only has entries up to the number of reflections
    ! found, INoOfHKLsFrame.  We ignore any duplicates in felix.hkl 
    IhklsFrame = 0! reset flag for the frame
    IhklsAll = 0!reset global flag
    IFind = 0
    DO ind = 1,INoOfHKLsAll! the reflections in felix.hkl
      DO jnd = 1,SIZE(Rhkl,DIM=1)! the beam pool
        IF(ABS(Rhkl(jnd,1)-RInputHKLs(ind,1)).LE.TINY.AND.&
           ABS(Rhkl(jnd,2)-RInputHKLs(ind,2)).LE.TINY.AND.&
           ABS(Rhkl(jnd,3)-RInputHKLs(ind,3)).LE.TINY) THEN
          ! this reflection is in both lists
          IDuplicate = 0! start from the assumption that it is not a duplicate
          DO knd = 1,IFind!check the list to see if we already have it
            IF (ABS(Rhkl(IhklsFrame(knd),1)-RInputHKLs(ind,1)).LE.TINY.AND.&
                ABS(Rhkl(IhklsFrame(knd),2)-RInputHKLs(ind,2)).LE.TINY.AND.&
                ABS(Rhkl(IhklsFrame(knd),3)-RInputHKLs(ind,3)).LE.TINY) THEN
              IDuplicate = 1!yes we do
              IF (IFrame.EQ.1) CALL message(LS,"Duplicate HKL found, ignoring: ",NINT(RInputHKLs(ind,:)) )
              EXIT
            END IF
          END DO
          IF (IDuplicate.EQ.0) THEN! not a duplicate, we append it to the list
            IFind = IFind +1
            IhklsFrame(IFind) = jnd! for this frame: the index of the reflection in the beam pool
            IhklsAll(IFind) = ind! so we can find it in the felix.hkl list
          END IF
          EXIT
        END IF
      END DO
    END DO

    IF (IFind.LE.0) THEN
      IErr=1
      IF(l_alert(IErr,"HKLList","No requested HKLs found")) RETURN
    END IF
      
    INoOfHKLsFrame = IFind
    IF(IFrame.GT.9999)THEN
      WRITE(SPrintString,'(I3,A22,I5.5)') INoOfHKLsFrame," reflections in Frame ",IFrame
    ELSEIF(IFrame.GT.999)THEN
      WRITE(SPrintString,'(I3,A22,I4.4)') INoOfHKLsFrame," reflections in Frame ",IFrame
    ELSEIF(IFrame.GT.99)THEN
      WRITE(SPrintString,'(I3,A22,I3.3)') INoOfHKLsFrame," reflections in Frame ",IFrame
    ELSEIF(IFrame.GT.9)THEN
      WRITE(SPrintString,'(I3,A22,I2.2)') INoOfHKLsFrame," reflections in Frame ",IFrame
    ELSE
      WRITE(SPrintString,'(I3,A22,I1.1)') INoOfHKLsFrame," reflections in Frame ",IFrame
    END IF
    CALL message(LM,TRIM(ADJUSTL(SPrintString)))
 
  END SUBROUTINE HKLList
  
  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Checks a g-vector Ih,Ik,Il
  !! against the selection rules for the global variable SSpaceGroupName
  !! IFlag comes in as zero and goes out as 1 if it is an allowed reflection
  !!
  !! Major-Authors: Richard Beanland (2021)
  !!  
  SUBROUTINE SelectionRules(Ih, Ik, Il, ISel, IErr)

    ! This procedure is called from HKLMake
    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE SPARA, ONLY : SSpaceGroupName
    
    IMPLICIT NONE

    INTEGER (IKIND),INTENT(IN) :: Ih, Ik, Il
    INTEGER (IKIND),INTENT(INOUT) :: ISel, IErr

    SELECT CASE(SSpaceGroupName)

      CASE("F") !Face Centred, all odd or all even
        IF( (MOD(Ih+Ik,2).EQ.0).AND.&
          (MOD(Ik+Il,2).EQ.0).AND.&
          (MOD(Il+Ih,2).EQ.0) ) ISel=1

      CASE("I")! Body Centred
        IF(MOD(Ih+Ik+Il,2).EQ.0) ISel=1

      CASE("A")! A-Face Centred
        IF(MOD(Ik+Il,2).EQ.0) ISel=1

      CASE("B")! B-Face Centred
        IF(MOD(Ih+Il,2).EQ.0) ISel=1

      CASE("C")! C-Face Centred
        IF(MOD(Ih+Ik,2).EQ.0) ISel=1

      CASE("R")! Rhombohedral Reverse
        IF(MOD(Ih-Ik+Il,3).EQ.0) ISel=1

      CASE("V")! Rhombohedral Obverse
        IF(MOD(-Ih+Ik+Il,3).EQ.0) ISel=1

      CASE("P")! Primitive
        ISel=1

      CASE DEFAULT
      IErr=1
      IF(l_alert(IErr,"SelectionRules",&
          "Space Group Name unrecognised")) RETURN
          
    END SELECT
     
  END SUBROUTINE SelectionRules
  
END MODULE setup_reflections_mod

