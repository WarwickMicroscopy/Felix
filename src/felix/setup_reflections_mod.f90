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
    USE myMPI

    ! global inputs/outputs
    USE SPARA, ONLY : SPrintString,SChemicalFormula
    USE IPARA, ONLY : INhkl,IgOutList,IgPoolList,IhklLattice,INFrames,InLattice,ILN,IByteSize,ISort
    USE RPARA, ONLY : RXDirO,RYDirO,RZDirO,RcrVecM,RLatMag,RFrameAngle,&
        RBigK,RgLatticeO,RgPoolSg
    USE CPARA, ONLY : CFg
    USE Iconst
    USE IChannels, ONLY : IChOutIhkl,IChOutIM
    
    IMPLICIT NONE

    REAL(RKIND),INTENT(IN) :: RDevLimit, RGOutLimit
    INTEGER(IKIND) :: IErr,ind,jnd,knd,lnd,mnd,ISim,Ix,Iy,ILocalFrameMin,ILocalFrameMax,&
                      ILocalNFrames
    INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: Inum,Ipos,ILocalgPool,ITotalgPool
    REAL(RKIND) :: RAngle,Rk(ITHREE),Rk0(ITHREE),Rp(ITHREE),RSg,Rphi,Rg(ITHREE),Rmos,RIkin,&
                   RKplusg(ITHREE)
    REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RSim
    REAL(RKIND), DIMENSION(:), ALLOCATABLE :: RLocalSgPool,RTotalSgPool
    CHARACTER(200) :: path
    CHARACTER(100) :: fString
   
    !--------------------------------------------------------------------
    ! calculate reflection list frame by frame
    !--------------------------------------------------------------------
    ! In the subroutine ReciprocalLattice we generated all reciprocal lattice vectors
    ! and put them in ascending order of magnitude RLatMag
    ! with indices of IhklLattice and vector RgLatticeO in the orthogonal ref frame
    ! IgPoolList says which reflections are close to the Ewald sphere
    ! IgOutList says which reflections are to be saved (|g|<RGOutLimit)
    IgPoolList = 0  ! Initialise lists to zero
    RgPoolSg = ZERO
    Rk0 = ZERO
    !--------------------------------------------------------------------
    ! set up frame-parallel calculations
    ! The frames to be calculated by this core
    ILocalFrameMin = (INFrames*(my_rank)/p)+1
    ILocalFrameMax = (INFrames*(my_rank+1)/p)
    ILocalNFrames = ILocalFrameMax-ILocalFrameMin+1
    ! Calculations are done in 1D arrays that are reshaped later
    ALLOCATE(ILocalgPool(INhkl*ILocalNFrames),STAT=IErr)
    IF(l_alert(IErr,"HKLmake","allocate ILocalgPool")) RETURN
    ALLOCATE(ITotalgPool(INhkl*INFrames),STAT=IErr)
    IF(l_alert(IErr,"HKLmake","allocate ITotalgPool")) RETURN
    ALLOCATE(RLocalSgPool(INhkl*ILocalNFrames),STAT=IErr)
    IF(l_alert(IErr,"HKLmake","allocate ILocalgPool")) RETURN
    ALLOCATE(RTotalSgPool(INhkl*INFrames),STAT=IErr)
    IF(l_alert(IErr,"HKLmake","allocate ITotalgPool")) RETURN
    ALLOCATE(Inum(p),Ipos(p),STAT=IErr)
    IF(l_alert(IErr,"HKLMake","allocate ILocalNhkl")) RETURN
    DO ind = 1,p
      Ipos(ind) = INhkl*INFrames*(ind-1)/p
      Inum(ind) = INhkl*(INFrames*ind/p - INFrames*(ind-1)/p)
    END DO
    DO ind = 1,ILocalNFrames
      RAngle = REAL(ILocalFrameMin+ind-2)*DEG2RADIAN*RFrameAngle
      ! Rk is the k-vector for the incident beam, which we write here in the orthogonal frame O
      Rk = RBigK*(RZDirO*COS(RAngle)+RXDirO*SIN(RAngle))
      ! Fill the list of reflections ILocalgPoolList until we have filled the beam pool
      knd = 1  ! size of beam pool for this frame
      DO mnd = 1,InLattice
        jnd = ISort(mnd)  ! work through reflections in ascending order
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
        ! and now Sg is 2g sin(phi/2), with the sign of K-|K+g|
        RKplusg = Rk + RgLatticeO(jnd,:)
        RSg = TWO*RLatMag(jnd)*SIN(HALF*Rphi)*SIGN(ONE,RBigK-SQRT(DOT_PRODUCT(RKplusg,RKplusg)))
        IF (ABS(RSg).LT.RDevLimit) THEN
          IF (knd.LE.INhkl) THEN ! while the beam pool isn't full
            ILocalgPool((ind-1)*INhkl+knd) = jnd  ! add the reflection to the list
            RLocalSgPool((ind-1)*INhkl+knd) = RSg  ! put Sg in its list also
            knd = knd + 1
          END IF
        END IF
      END DO
    END DO

    !==================== ! MPI gatherv into 1D arrays ========================
    CALL MPI_GATHERV(ILocalgPool,SIZE(ILocalgPool),MPI_INTEGER,ITotalgPool,&
            Inum,Ipos,MPI_INTEGER,root,MPI_COMM_WORLD,IErr)
    CALL MPI_GATHERV(RLocalSgPool,SIZE(RLocalSgPool),MPI_DOUBLE_PRECISION,RTotalSgPool,&
            Inum,Ipos,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,IErr)
    IgPoolList = RESHAPE(ITotalgPool, (/INFrames,INhkl/) )
    RGPoolSg = RESHAPE(RTotalSgPool, (/INFrames,INhkl/) )

    ! fill IgOutList & output as a text file and a frame series
    IF(my_rank.EQ.0) THEN
      CALL message(LS,dbg3,"Writing hkl list and images")
      path = SChemicalFormula(1:ILN) // "/hkl_list.txt"
      OPEN(UNIT=IChOutIhkl, ACTION='WRITE', POSITION='APPEND', STATUS= 'UNKNOWN', &
          FILE=TRIM(ADJUSTL(path)),IOSTAT=IErr)
      WRITE(IChOutIhkl,*) "List of hkl in each frame"
      WRITE(IChOutIhkl,*) "No: h k l  Fg  |g|  Sg"
    END IF
    
    IgOutList = 0
    DO ind = 1,INFrames
      ! output to slurm if requested
      CALL message(LM, "Reflection list:")
      DO knd = 1, INhkl
        IF (IgPoolList(ind,knd).NE.0) THEN
          WRITE(SPrintString,'(I3,1X,I3,1X,I3)') IhklLattice(IgPoolList(ind,knd),:)
          CALL message(LM, SPrintString)
        END IF
      END DO
      IF (my_rank.EQ.0) WRITE(IChOutIhkl,"(A6,I4)") "Frame ",ind
      DO knd = 1, INhkl
        IF (IgPoolList(ind,knd).NE.0) THEN
          ! Is this reflection small enough to be in the output list
          IF (RLatMag(IgPoolList(ind,knd)).LT.RGOutLimit) THEN
            IgOutList(ind,knd) = IgPoolList(ind,knd)
          END IF
          lnd = IgPoolList(ind,knd)
          RSg = RgPoolSg(ind,knd)
          WRITE(fString,"(I6,A2, 3(I3,1X),2X, F8.4,A1,F8.4,A3, F6.2,2X, F8.4)") &
                  lnd,": ",IhklLattice(lnd,:), REAL(CFg(lnd)),"+",AIMAG(CFg(lnd)),"i  ",&
                  RLatMag(lnd)/TWOPI, RSg
          IF (my_rank.EQ.0) WRITE(IChOutIhkl,*) TRIM(ADJUSTL(fString))
        END IF
      END DO
    END DO
    IF (my_rank.EQ.0) CLOSE(IChOutIhkl,IOSTAT=IErr)

    IF(my_rank.EQ.0) THEN
      ! image
      ISim = 256_IKIND  ! NB HALF the output image size = output |g| limit/0.98
      ALLOCATE(RSim(2*ISim,2*ISim),STAT=IErr)
      IF(l_alert(IErr,"HKLmake","allocate RSim")) RETURN
      ! Mosaicity - sets the FWHM  of a kinematic rocking curve
      Rmos = 3000.0
      DO ind = 1,INFrames
        RAngle = REAL(ind-1)*DEG2RADIAN*RFrameAngle
        RSim = ZERO
        ! Direct beam
        RSim(ISim-1:ISim+1,ISim-1:ISim+1) = 1
        ! output g's
        DO knd = 1, INhkl
          IF (IgOutList(ind,knd).NE.0) THEN
            lnd = IgPoolList(ind,knd)  ! index of reflection in the reciprocal lattice
            Rg = RgLatticeO(lnd,:)  ! g-vector
            RIkin = CFg(lnd)*CONJG(CFg(lnd))  ! simple kinematic intensity
            RSg = RgPoolSg(ind,knd)  !Sg
            ! x- and y-coords (NB swapped in the image!)
            Rp = RXDirO*COS(RAngle)-RZDirO*SIN(RAngle)  ! unit vector horizontal in the image
            ! position of the spot, 2% leeway to avoid going over the edge of the image
            Ix = ISim-0.98*NINT(DOT_PRODUCT(Rg,Rp)*REAL(ISim)/RGOutLimit)  
            Iy = ISim+0.98*NINT(DOT_PRODUCT(Rg,RYDirO)*REAL(ISim)/RGOutLimit)
            RSim(Iy-1:Iy+1,Ix-1:Ix+1) = EXP(-RMos*RSg*RSg)*RIkin
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

    ! Clean up
    DEALLOCATE(ILocalgPool,ITotalgPool,Rsim,RLocalSgPool,RTotalSgPool)

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

