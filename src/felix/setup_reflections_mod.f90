!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-17, all rights reserved
!
! Version: 2.0
! Date: 31-08-2022
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
  PUBLIC :: HKLList, SelectionRules, HKLMake, HKLSort

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
    USE RPARA, ONLY : RInputHKLs, RTolerance, Rhkl

    IMPLICIT NONE

    INTEGER(IKIND) :: IFind,IFound,ind,jnd,knd,IErr
      
    IFind = 0

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
  
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Fills the list of reciprocal space vectors Rhkl
  !!
  !! Major-Authors: Richard Beanland (2021)
  !!  
  SUBROUTINE HKLMake(RGlimit, IErr)   
    ! This procedure is called once in felixrefine setup

    USE MyNumbers
    USE message_mod
    
    ! global output Rhkl, input reciprocal lattice vectors & wave vector
    USE RPARA, ONLY : RzDirC, Rhkl, RarVecM, RbrVecM, RcrVecM, RInputHKLs,&
        RElectronWaveVectorMagnitude
      
    ! global inputs
    USE SPARA, ONLY : SSpaceGroupName
    USE IPARA, ONLY : INhkl,INoOfLacbedPatterns
    USE Iconst
    
    IMPLICIT NONE

    REAL(RKIND),INTENT(IN) :: RGlimit   
    INTEGER(IKIND) :: IErr, ISel, Ih, Ik, Il, inda,indb,indc, jnd, knd,lnd
    REAL(RKIND) :: RarMag, RbrMag, RcrMag, RShell, RGtestMag, RDev
    REAL(RKIND),DIMENSION(ITHREE) :: Rk, RGtest, RGtestM, RGplusk 
   
    !The upper limit for g-vector magnitudes
    !If the g-vectors we are counting are bigger than this there is something wrong
    !probably the tolerance for proximity to the Ewald sphere needs increasing
    !could be an input in felix.inp
!    RGlimit = 20.0*TWOPI  ! reciprocal Angstroms * 2pi
    
    !the k-vector for the incident beam
    !we are working in the microscope reference frame so k is along z
    Rk=(/ 0.0,0.0,RElectronWaveVectorMagnitude /)
    !get the size of the reciprocal lattice basis vectors
    RarMag=SQRT(DOT_PRODUCT(RarVecM,RarVecM))!magnitude of a*
    RbrMag=SQRT(DOT_PRODUCT(RbrVecM,RbrVecM))!magnitude of b*
    RcrMag=SQRT(DOT_PRODUCT(RcrVecM,RcrVecM))!magnitude of c*
    !we work our way out from the origin in shells
    !the shell increment is the smallest basis vector
    RShell=MINVAL( (/ RarMag,RbrMag,RcrMag /) )

    !first g is always 000
    Rhkl(1,:)=(/ 0.0,0.0,0.0 /)
    knd=1!number of reflections in the pool 
    lnd=0!number of the shell 
    !maximum a*,b*,c* limit is determined by the G magnitude limit
    inda=NINT(RGlimit/RarMag)
    indb=NINT(RGlimit/RbrMag)
    indc=NINT(RGlimit/RcrMag)
    !fill the Rhkl with beams near the Bragg condition
    DO WHILE (knd.LT.INhkl .AND. REAL(lnd)*RShell.LT.RGlimit)
      !increment the shell
      lnd = lnd+1
!DBG    IF(my_rank.EQ.0)PRINT*,REAL(lnd-1)*RShell,"to",REAL(lnd)*RShell
      !Make a hkl
      DO Ih = -inda,inda
         DO Ik = -indb,indb
            DO Il = -indc,indc
              !check that it's allowed by selection rules
              ISel=0
              CALL SelectionRules(Ih, Ik, Il, ISel, IErr)
              !NB still need to check that we have space in Rhkl because the
              !while doesn't kick in until we finish the do loops
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
                  ! Tolerance of 0.08 here is rather arbitrary, might need
                  ! revisiting
!DBG    IF(my_rank.EQ.0)PRINT*,(/ Ih,Ik,Il /),RDev
                  IF ((RDev-0.08).LT.TINY) THEN !it's near the ZOLZ 
                    !add it to the pool and increment the counter
                    !need to check that we have space in Rhkl because the while doesn't
                    !kick in until we finish the do loops
                    knd=knd+1
                    Rhkl(knd,:)=RGtest
!DBG    IF(my_rank.EQ.0)PRINT*,NINT(Rhkl(knd,:)),knd
                  END IF
                END IF
              END IF
            END DO
         END DO
      END DO
    END DO
!DBG    IF(my_rank.EQ.0)PRINT*,"total ",knd,"reflections in the pool"
    ! fill up any remaining beam pool places with an enormous g-vector, diabolically
    ! the idea being that this g-vector will never be near any possible Ewald sphere
    RGtest = REAL( (/ 666,666,666 /),RKIND )
    DO jnd = knd+1, INhkl
      Rhkl(jnd,:)=RGtest
    END DO

    !now check that the required output hkl's are in this hkl list
    DO jnd = 1,INoOfLacbedPatterns
      lnd = 0!using this as a flag now
      DO knd = 1, INhkl
        IF (RInputHKLs(jnd,1)-Rhkl(knd,1).LT.TINY .AND. &
            RInputHKLs(jnd,2)-Rhkl(knd,2).LT.TINY .AND. &
            RInputHKLs(jnd,3)-Rhkl(knd,3).LT.TINY) lnd = 1
      END DO
      IF (lnd.EQ.0) THEN
        CALL message( LS, "Input hkl not found: ", NINT(RInputHKLs(jnd,:)))
      END IF
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
