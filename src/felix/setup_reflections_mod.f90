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
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Fills the beam pool list RgPoolList for each frame
  !! and the list of output reflections IgOutList (global variables)
  !!
  !! Major-Authors: Richard Beanland (2023)
  !!  
  SUBROUTINE HKLMake(RDevLimit, RGOutLimit, IErr)   

    USE MyNumbers
    USE message_mod

    ! global inputs/outputs
    USE SPARA, ONLY : SPrintString,SChemicalFormula
    USE IPARA, ONLY : INhkl,IgOutList,IgPoolList,IhklLattice,INFrames,InLattice,ILN,IByteSize
    USE RPARA, ONLY : RXDirO,RYDirO,RZDirO,RcrVecM,RLatMag,RFrameAngle,&
        RBigK,RgLatticeO,RgPoolSg
    USE Iconst
    USE IChannels, ONLY : IChOutIhkl,IChOutIM
    
    IMPLICIT NONE

    REAL(RKIND),INTENT(IN) :: RDevLimit, RGOutLimit
    INTEGER(IKIND) :: IErr,ind,jnd,knd,lnd,ISim,Ix,Iy
    REAL(RKIND) :: RAngle,Rk(ITHREE),Rk0(ITHREE),Rp(ITHREE),RSg,Rphi,Rg(ITHREE),Rmos
    REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RSim
    CHARACTER(200) :: path
    CHARACTER(40) :: fString
   
    !--------------------------------------------------------------------
    ! calculate reflection list frame by frame
    !--------------------------------------------------------------------
    ! In the subroutine ReciprocalLattice we generated all reciprocal lattice vectors
    ! and put them in ascending order of magnitude RLatMag
    ! with indices of IhklLattice and vector RgLatticeO in the orthogonal ref frame
    ! IgPoolList says which reflections are close to the Ewald sphere
    IgPoolList = 0  ! Initialise lists to zero
    IgOutList = 0
    RgPoolSg = ZERO
    Rk0 = ZERO
    DO ind = 1,INFrames
      WRITE(SPrintString, FMT='(A30,I3,A3)') "Counting reflections in frame ",ind,"..."
      CALL message(LS,dbg3,SPrintString)
      RAngle = REAL(ind-1)*DEG2RADIAN*RFrameAngle
      ! Rk is the k-vector for the incident beam, which we write here in the orthogonal frame O
      Rk = RBigK*(RZDirO*COS(RAngle)+RXDirO*SIN(RAngle))
      ! Fill the list of reflections IgPoolList until we have filled the beam pool
      knd = 1
      DO jnd = 1,InLattice  ! work through reflections in ascending order
        ! Calculate Sg by getting the vector k0, which is coplanar with k and g and
        ! corresponds to an incident beam at the Bragg condition
        ! First we need the vector component of k perpendicular to g, which we call p 
        Rp = Rk - DOT_PRODUCT(Rk,RgLatticeO(jnd,:))*RgLatticeO(jnd,:)/(RLatMag(jnd)**2)
        ! and now make k0 by adding vectors parallel to g and p
        ! i.e. k0 = (p/|p|)*(k^2-g^2/4)^0.5 - g/2
        Rk0 = SQRT(RBigK**2-QUARTER*RLatMag(jnd)**2)*Rp/SQRT(DOT_PRODUCT(Rp,Rp)) - &
              HALF*RgLatticeO(jnd,:)
        ! The angle phi between k and k0 is how far we are from the Bragg condition
        Rphi = ACOS(DOT_PRODUCT(Rk,Rk0)/(RBigK**2))
        ! and now Sg is 2g sin(phi/2)
        RSg = TWO*RLatMag(jnd)*SIN(HALF*Rphi)
        IF (ABS(RSg).LT.RDevLimit) THEN
          IF (knd.LE.INhkl) THEN ! while the beam pool isn't full
            IgPoolList(ind,knd) = jnd  ! add it to the list
            RgPoolSg(ind,knd) = RSg  ! put Sg in its list also
            ! Is this reflection small enough to be in the output list
            IF (RLatMag(jnd).LT.RGOutLimit) THEN
              IgOutList(ind,knd) = jnd
            END IF
            knd = knd + 1
          END IF
        END IF
      END DO
      IF(my_rank.EQ.0)PRINT*,"Found",knd-1,"reflections"

      CALL message(LM, "Reflection list:")
      DO knd = 1, INhkl
        IF (IgPoolList(ind,knd).NE.0) THEN
          WRITE(SPrintString,'(I3,1X,I3,1X,I3)') IhklLattice(IgPoolList(ind,knd),:)
          CALL message(LM, SPrintString)
        END IF
      END DO
    END DO

    !output the hkl lists for the frames as a text file and an image
    IF(my_rank.EQ.0) THEN
    
      ! text file
      path = SChemicalFormula(1:ILN) // "/hkl_list.txt"
      OPEN(UNIT=IChOutIhkl, ACTION='WRITE', POSITION='APPEND', STATUS= 'UNKNOWN', &
          FILE=TRIM(ADJUSTL(path)),IOSTAT=IErr)
      WRITE(IChOutIhkl,*) "List of hkl in each frame"
      WRITE(IChOutIhkl,*) "No: h k l |g| Sg"
      DO ind = 1,INFrames
        WRITE(IChOutIhkl,"(A6,I4)") "Frame ",ind
        DO knd = 1, INhkl
          ! version writing output g's only
          IF (IgOutList(ind,knd).NE.0) THEN
            lnd = IgPoolList(ind,knd)
            RSg = RgPoolSg(ind,knd)
            WRITE(fString,"(I4,A2,I3,1X,I3,1X,I3,2X,F8.4,2X,F8.4,2X,F8.2)") &
                    lnd,": ",IhklLattice(lnd,:),RLatMag(lnd)/TWOPI,RSg
          END IF
          ! version writing the full beam pool
          !IF (IgPoolList(ind,knd).NE.0) THEN
          !  lnd = IgPoolList(ind,knd)
          !  RSg = RgPoolSg(ind,knd)
          !  IF (IgOutList(ind,knd).NE.0) THEN  ! output g's are indicated by a *
          !    WRITE(fString,"(I4,A2,I3,1X,I3,1X,I3,2X,F8.4,2X,F8.4,2X,F8.2)") &
          !          lnd,"* ",IhklLattice(lnd,:),RLatMag(lnd)/TWOPI,RSg
          !  ELSE
          !    WRITE(fString,"(I4,A2,I3,1X,I3,1X,I3,2X,F8.4,2X,F8.4)") lnd, ": ", IhklLattice(lnd,:),RSg
          !  END IF
          !  WRITE(IChOutIhkl,*) TRIM(ADJUSTL(fString))
          !END IF
        END DO
      END DO
      CLOSE(IChOutIhkl,IOSTAT=IErr)
      
      ! image
      ISim = 256_IKIND  ! NB HALF the output image size = output |g| limit
      ALLOCATE(RSim(2*ISim,2*ISim),STAT=IErr)
      IF(l_alert(IErr,"HKLmake","allocate RSim")) RETURN
      ! Mosaicity - sets the FWHM  of a kinematic rocking curve
      Rmos = 800.0
      DO ind = 1,INFrames
        RSim = ZERO
        ! Direct beam
        RSim(ISim-1:ISim+1,ISim-1:ISim+1) = 1
        ! output g's
        DO knd = 1, INhkl
          IF (IgOutList(ind,knd).NE.0) THEN
            Rg = RgLatticeO(IgPoolList(ind,knd),:)
            RSg = RgPoolSg(ind,knd)
            ! x- and y-coords in the image are swapped 
            Iy = NINT(DOT_PRODUCT(Rg,RXDirO)*REAL(ISim)/RGOutLimit)  
            Ix = -NINT(DOT_PRODUCT(Rg,RYDirO)*REAL(ISim)/RGOutLimit)
            RSim(ISim+Ix-1:ISim+Ix+1,ISim+Iy-1:ISim+Iy+1) = EXP(-RMos*RSg*RSg)
          END IF
        END DO
        ! write to disk
        IF (ind.LT.10) THEN
          WRITE(path, FMT="(A,A15,I1,A4)") TRIM(ADJUSTL(SChemicalFormula(1:ILN))),&
                  "/Simulations/S_",ind,".bin"
        ELSE IF (ind.LT.100) THEN
          WRITE(path, FMT="(A,A15,I2,A4)") TRIM(ADJUSTL(SChemicalFormula(1:ILN))),&
                  "/Simulations/S_",ind,".bin"
        ELSE IF (ind.LT.1000) THEN
          WRITE(path, FMT="(A,A15,I3,A4)") TRIM(ADJUSTL(SChemicalFormula(1:ILN))),&
                  "/Simulations/S_",ind,".bin"
        ELSE
          WRITE(path, FMT="(A,A15,I4,A4)") TRIM(ADJUSTL(SChemicalFormula(1:ILN))),&
                  "/Simulations/S_",ind,".bin"
        END IF
        OPEN(UNIT=IChOutIM, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(path)),&
          FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=2*ISim*IByteSize)
        IF(l_alert(IErr,"WriteIterationOutput","OPEN() output .bin file")) RETURN      
        DO jnd = 1,2*ISim
          WRITE(IChOutIM,rec=jnd) RSim(jnd,:)
        END DO
        CLOSE(IChOutIM,IOSTAT=IErr) 
        IF(l_alert(IErr,"WriteIterationOutput","CLOSE() output .bin file")) RETURN
      END DO      
    END IF

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
  !! Procedure-description: List the reflections in each frame that form the beam pool and output
  !!
  !! Major-Authors: Richard Beanland (2023)
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    ! it is possible that we reach RGPoolLimit before filling up the pool, in which case 
    ! fill up any remaining beam pool places with an enormous g-vector, diabolically
    ! the idea being that this g-vector will never be near any possible Ewald sphere
!    IF (knd.LE.INhkl) THEN
!      DO jnd = knd+1, INhkl
!        IhklLattice(jnd,:) = (/ 666,666,666 /)
!        RgLatticeO(jnd,:) = REAL( (/ 666,666,666 /),RKIND )
!        RLatMag(jnd) = 666666.666
!      END DO
!    END IF    
    
    
    
    
    
    
    
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
  
  
END MODULE setup_reflections_mod

