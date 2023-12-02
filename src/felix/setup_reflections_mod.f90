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
  PUBLIC :: HKLMake

  CONTAINS
  
  !!$%%HKLMake%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description:
  !! 1) Fills the beam pool list RgPoolList for each frame (global variable)
  !! 2) Fills the list of output reflections IgOutList & writes to hkl_list.txt
  !! 3) Makes a set of simple kinematic frames
  !!
  !! Major-Authors: Richard Beanland (2023)
  !!  
  SUBROUTINE HKLMake(RDevLimit, RGOutLimit, RgPoolLimit, IErr)   

    USE MyNumbers
    USE message_mod
    USE myMPI
    USE crystallography_mod

    ! global inputs/outputs
    USE SPARA, ONLY : SPrintString,SChemicalFormula
    USE IPARA, ONLY : INhkl,ILN,INFrames,ICurrentZ,INAtomsUnitCell,IAtomicNumber  ! inputs
    USE IPARA, ONLY : Ig,IgOutList,IgPoolList  ! outputs
    USE RPARA, ONLY : RXDirO,RYDirO,RZDirO,RarVecO,RbrVecO,RcrVecO,RarMag,RbrMag,RcrMag,RFrameAngle,RBigK,&
        RAtomCoordinate,RIsoDW, RgPoolSg ! only RgPoolSg is an output
    USE Iconst
    USE IChannels, ONLY : IChOutIhkl,IChOutIM
    
    IMPLICIT NONE

    REAL(RKIND),INTENT(IN) :: RDevLimit, RGOutLimit, RgPoolLimit
    INTEGER(IKIND) :: IErr,ind,jnd,knd,lnd,mnd,nnd,ond,ISel,Ifull(INFrames),IMaxNg,inda,indb,indc,Ifound
    REAL(RKIND) :: RAngle,Rk(INFrames,ITHREE),Rk0(ITHREE),Rp(ITHREE),RSg,Rphi,Rg(ITHREE),RIkin,&
                   RKplusg(ITHREE),RgMag,RShell,RgMin,Rfq
    COMPLEX :: CFg
    CHARACTER(200) :: path
    CHARACTER(100) :: fString
   
    !-1------------------------------------------------------------------
    ! calculate reflection list g by g
    !--------------------------------------------------------------------
    ! this produces a list of g-vectors Ig, their deviation parameters RgPoolSg, and fills IgPoolList, IgOutList
    ! Initialise variables
    Rk = ZERO  ! list of incident beam k-vectors
    Ig = 0  ! list of reflections (covers all beam pools)
    IgPoolList = 0  ! index giving Ig for each beam pool (1 to INhkl,1 to INFrames)
    IgOutList = 0  ! index giving Ig for each output (1 to INhkl,1 to INFrames)
    IFull = 0  ! is the beam pool full
    RgPoolSg = TEN  !Sg corresponding to Ig (1 to INhkl,1 to INFrames)
    ! we make a list of the incident k-vectors for each frame
    DO ind = 1,INFrames
      RAngle = REAL(ind-1)*DEG2RADIAN*RFrameAngle
      ! Rkk is the k-vector for the incident beam, which we write here in the orthogonal frame O
      Rk(ind,:) = RBigK*(RZDirO*COS(RAngle)+RXDirO*SIN(RAngle))
      ! the 000 beam is the first g-vector in every frame
    END DO
    ! the 000 beam is the first g-vector in every frame
    Ig(1,:) = (/0,0,0/)
    IgPoolList(1,:) = 1
    IgOutList(1,:) = 1
    lnd = 1  ! index counting reflexions as they are added to the list Ig   
    ! we work our way out in shells of the smallest reciprocal lattice vector
    RShell = MINVAL( (/RarMag,RbrMag,RcrMag/) )
    mnd = 1 ! shell count

    ! start of the loop of incrementing shells
1   inda = NINT(REAL(mnd)*RShell/RarMag)
    indb = NINT(REAL(mnd)*RShell/RbrMag)
    indc = NINT(REAL(mnd)*RShell/RcrMag)
    RgMin = 666.666  ! smallest gMag checker set to initial large value 
    DO ind = -inda,inda
      DO jnd = -indb,indb
        DO knd = -indc,indc
          ! take out systematic absences from the lattice
          ! but keep forbidden reflections because they may appear through multiple scattering
          ISel = 0
          CALL SelectionRules(ind, jnd, knd, ISel, IErr)  ! in this module
          IF (ISel.EQ.0) CYCLE
          Ifound = 1  ! flag to indicate this reflexion is active
          Rg = ind*RarVecO + jnd*RbrVecO + knd*RcrVecO
          RgMag = SQRT(DOT_PRODUCT(Rg,Rg))
          ! Is this reflexion in the current shell
          IF (RgMag.GT.REAL(mnd-1)*RShell .AND. RgMag.LE.REAL(mnd)*RShell .AND. &
                  RgMag.LE.RgPoolLimit) THEN
            ! is it the smallest g in this shell
            IF (RgMag.LT.RgMin) THEN
              RgMin = RgMag
            END IF

            ! go through the frames and see if it appears
            ! probably a more elegant/faster way of doing this using matrices rather than a loop
            DO nnd = 1,INFrames
              ! Is the beam pool already full for this frame
              IF (Ifull(nnd).EQ.1) CYCLE
              ! Calculate Sg by getting the vector k0, which is coplanar with k and g and
              ! corresponds to an incident beam at the Bragg condition
              ! First we need the vector component of k perpendicular to g, which we call p 
              Rp = Rk(nnd,:) - DOT_PRODUCT(Rk(nnd,:),Rg)*Rg/RgMag
              ! and now make k0 by adding vectors parallel to g and p
              ! i.e. k0 = (p/|p|)*(k^2-g^2/4)^0.5 - g/2
              Rk0 = SQRT(RBigK**2-QUARTER*RgMag**2)*Rp/SQRT(DOT_PRODUCT(Rp,Rp)) - HALF*Rg
              ! The angle phi between k and k0 is how far we are from the Bragg condition
              Rphi = ACOS(DOT_PRODUCT(Rk(nnd,:),Rk0)/(RBigK**2))
              ! and now Sg is 2g sin(phi/2), with the sign of K-|K+g|
              RKplusg = Rk(nnd,:) + Rg
              RSg = TWO*RgMag*SIN(HALF*Rphi)*SIGN(ONE,RBigK-SQRT(DOT_PRODUCT(RKplusg,RKplusg)))

              ! now check to see if it qualifies for the beam pool and output
              IF (ABS(RSg).LT.RDevLimit) THEN  ! it's in the beam pool for this frame
                lnd = lnd+Ifound  ! increment the reflexion counter on the first occurrance only
                Ifound = 0
                Ig(lnd,:) = (/ind,jnd,knd/)  ! add it to the list of reflexions
                ond = 2 !counter to find the next slot for the g pool
2               IF (IgPoolList(ond,nnd).NE.0) THEN
                  ond = ond + 1
                  IF (ond.EQ.INhkl) THEN  ! we have filled the beam pool
                    Ifull(nnd) = 1
                    IF (my_rank.EQ.0)PRINT*,"Beam pool full for frame",nnd
                    CYCLE
                  ELSE
                    GOTO 2
                  END IF
                ELSE
                  IgPoolList(ond,nnd) = lnd
                  RgPoolSg(ond,nnd) = RSg
                  IF (ABS(RSg).LT.RGOutLimit) THEN
                    IgOutList(ond,nnd) = lnd
                  END IF
                END IF
              END IF
            END DO
          END IF
        END DO
      END DO
    END DO
    IF (RgMin.GT.RgPoolLimit.AND.my_rank.EQ.0)PRINT*,"G-pool limit reached at",RgMin
    IF (lnd.LT.INhkl.AND.RgMin.LE.RgPoolLimit) THEN
      mnd = mnd + 1
      GOTO 1
    END IF

  END SUBROUTINE HKLmake

  !!$%%HKLList%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: List the frames for each reflection
  !!
  !! Major-Authors: Richard Beanland (2023)
  !!  
  SUBROUTINE HKLList( IErr )

  USE MyNumbers

  IMPLICIT NONE

  INTEGER (IKIND) :: IErr
 
  END SUBROUTINE HKLList
  
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
  
END MODULE setup_reflections_mod


