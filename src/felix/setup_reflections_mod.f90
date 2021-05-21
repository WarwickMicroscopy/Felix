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
      IF(l_alert(IErr,"HKLCount",&
          "SSpaceGroupName unrecognised")) RETURN
          
    END SELECT
     
  END SUBROUTINE SelectionRules
  
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Fills the list of reciprocal space vectors Rhkl
  !!
  !! Major-Authors: Richard Beanland (2021)
  !!  
  SUBROUTINE HKLMake(RTol, IErr)   
    ! This procedure is called once in felixrefine setup
    ! working out from the origin in reciprocal space it fills in Rhkl using reflections
    ! that are closer to the Ewald sphere than the tolerance
    ! RTol in reciprocal Angstroms

    USE MyNumbers
    USE message_mod
    
    ! global output Rhkl, input reciprocal lattice vectors & wave vector
    USE RPARA, ONLY : Rhkl,RarVecM,RbrVecM,RcrVecM,RElectronWaveVectorMagnitude
      
    ! global inputs
    USE SPARA, ONLY : SSpaceGroupName
    USE IPARA, ONLY : INhkl
    USE Iconst
    
    IMPLICIT NONE

    REAL(RKIND),INTENT(IN) :: RTol   
    INTEGER(IKIND) :: IErr, ISel, Ih, Ik, Il, inda,indb,indc, jnda,jndb,jndc, knd,lnd
    REAL(RKIND) :: RarMag, RbrMag, RcrMag, RGtestMag, RGpluskMag, RGlimit, RShell
    REAL(RKIND),DIMENSION(ITHREE) :: RGtest, RGplusk, Rk
    
    !The upper limit for g-vector magnitudes
    !If the g-vectors we are counting are bigger than this there is something wrong
    !probably the tolerance for proximity to the Ewald sphere needs increasing
    !could be an input in felix.inp
    RGlimit = 10.0*TWOPI  ! reciprocal Angstroms * 2pi
    
    !the k-vector for the incident beam
    !we are working in the microscope reference frame so k is along z
    Rk=(/ 0.0,0.0,RElectronWaveVectorMagnitude /)
    
    !get the size of the reciprocal lattice basis vectors
    RarMag=SQRT(DOT_PRODUCT(RarVecM,RarVecM))!magnitude of a*
    RbrMag=SQRT(DOT_PRODUCT(RbrVecM,RbrVecM))!magnitude of b*
    RcrMag=SQRT(DOT_PRODUCT(RcrVecM,RcrVecM))!magnitude of c*
    !limit the shell increment by the smallest basis vector
    RShell=MINVAL( (/ RarMag,RbrMag,RcrMag /) )

    !we work our way out from the origin in shells
    Rhkl(1,:)=(/ 0.0,0.0,0.0 /)
    knd=1!current number of reflections in the pool 
    lnd=0!current number of times we've expanded the search 
    ! first shell
    RGtestMag=0.0!we check this against RGlimit
    !current a*,b*,c* limits
    inda=1
    indb=1
    indc=1
    !previous a*,b*,c* limits
    jnda=0
    jndb=0
    jndc=0
    DO WHILE (knd.LT.INhkl.AND.RGtestMag.LE.RGlimit)
      lnd=lnd+1
!dbg IF(my_rank.EQ.0)PRINT*,lnd
      DO Ih=-inda,inda
         DO Ik=-indb,indb
            DO Il=-indc,indc
              !skip if we've already looked at this g
              IF (abs(Ih).GT.jnda .OR. abs(Ik).GT.jndb .OR. abs(Il).GT.jndc) THEN
                !Make a g-vector hkl and add it to k
                RGtest=REAL(Ih)*RarVecM+REAL(Ik)*RbrVecM+REAL(Il)*RcrVecM
                RGtestMag=SQRT(DOT_PRODUCT(RGtest,RGtest))
                RGplusk=RGtest+Rk
                !does it satisfy the Laue condition |k+g|=|k| within tolerance?
                RGpluskMag=SQRT(DOT_PRODUCT(RGplusk,RGplusk))
                IF (ABS(RGpluskMag-RElectronWaveVectorMagnitude).LT.RTol) THEN
                  !check that it's allowed by selection rules
                  ISel=0
                  CALL SelectionRules(Ih, Ik, Il, ISel, IErr)
                  !add it to the pool and increment the counter
                  !need to check that we have space in Rhkl because the while doesn't
                  !kick in until we finish the do loops
                  IF (ISel.EQ.1 .AND. knd.LT.INhkl) THEN
                    knd=knd+1
                    Rhkl(knd,:)=REAL((/ Ih,Ik,Il /),RKIND)
!dbg IF(my_rank.EQ.0)PRINT*,Ih,Ik,Il,RGpluskMag-RElectronWaveVectorMagnitude,RGtestMag
                  END IF
                END IF
              END IF
            END DO
         END DO
      END DO
      !expand the range for the next shell
      !update the limits of the previous shell
      jnda=inda
      jndb=indb
      jndc=indc
      !limits for the next shell
      inda=NINT(REAL(lnd)*RarMag/RShell)
      indb=NINT(REAL(lnd)*RbrMag/RShell)
      indc=NINT(REAL(lnd)*RcrMag/RShell)
    END DO
IF(my_rank.EQ.0)PRINT*,"total ",knd,"reflections in the pool"

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
