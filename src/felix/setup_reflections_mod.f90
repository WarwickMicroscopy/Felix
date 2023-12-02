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
  PUBLIC :: HKLMake,HKLList

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
  SUBROUTINE HKLMake(RDevLimit, RGOutLimit, IErr)   

    USE MyNumbers
    USE message_mod
    USE myMPI

    ! global inputs/outputs
    USE SPARA, ONLY : SPrintString,SChemicalFormula
    USE IPARA, ONLY : INhkl,IgOutList,IgPoolList,IhklLattice,INFrames,InLattice,ILN,IByteSize,ISort
    USE RPARA, ONLY : RXDirO,RYDirO,RZDirO,RcrVecM,RgMagLattice,RFrameAngle,&
        RBigK,RgLatticeO,RgPoolSg
    USE CPARA, ONLY : CFgLattice
    USE Iconst
    USE IChannels, ONLY : IChOutIhkl,IChOutIM
    
    IMPLICIT NONE

    REAL(RKIND),INTENT(IN) :: RDevLimit, RGOutLimit
    INTEGER(IKIND) :: IErr,ind,jnd,knd,lnd,mnd,ISim,Ix,Iy,ILocalFrameMin,ILocalFrameMax,&
                      ILocalNFrames,IMaxNg
    INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: Inum,Ipos,ILocalgPool,ITotalgPool
    REAL(RKIND) :: RAngle,Rk(ITHREE),Rk0(ITHREE),Rp(ITHREE),RSg,Rphi,Rg(ITHREE),RInst,RIkin,&
                   RKplusg(ITHREE)
    REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RSim
    REAL(RKIND), DIMENSION(:), ALLOCATABLE :: RLocalSgPool,RTotalSgPool
    CHARACTER(200) :: path
    CHARACTER(100) :: fString
   
    !-1------------------------------------------------------------------
    ! calculate reflection list frame by frame
    !--------------------------------------------------------------------
    ! In the subroutine ReciprocalLattice we generated all reciprocal lattice vectors
    ! and put them in ascending order of magnitude RgMagLattice
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
    IF(l_alert(IErr,"HKLMake","allocate MPI variables")) RETURN
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
        Rp = Rk - DOT_PRODUCT(Rk,RgLatticeO(jnd,:))*RgLatticeO(jnd,:)/(RgMagLattice(jnd)**2)
        ! and now make k0 by adding vectors parallel to g and p
        ! i.e. k0 = (p/|p|)*(k^2-g^2/4)^0.5 - g/2
        Rk0 = SQRT(RBigK**2-QUARTER*RgMagLattice(jnd)**2)*Rp/SQRT(DOT_PRODUCT(Rp,Rp)) - &
              HALF*RgLatticeO(jnd,:)
        ! The angle phi between k and k0 is how far we are from the Bragg condition
        Rphi = ACOS(DOT_PRODUCT(Rk,Rk0)/(RBigK**2))
        ! and now Sg is 2g sin(phi/2), with the sign of K-|K+g|
        RKplusg = Rk + RgLatticeO(jnd,:)
        RSg = TWO*RgMagLattice(jnd)*SIN(HALF*Rphi)*SIGN(ONE,RBigK-SQRT(DOT_PRODUCT(RKplusg,RKplusg)))
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
    !This broadcast is not strictly necessary but keeps all cores synchronised
    CALL MPI_BCAST(ITotalgPool,SIZE(ITotalgPool),MPI_INTEGER,root,MPI_COMM_WORLD,IErr)
    CALL MPI_BCAST(RTotalSgPool,SIZE(RTotalSgPool),MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,IErr)
    !IgPoolList = RESHAPE(ITotalgPool, (/INhkl,INFrames/) )
    !RGPoolSg = RESHAPE(RTotalSgPool, (/INhkl,INFrames/) )
    knd = 1
    DO ind = 1,INhkl
      DO jnd = 1,INFrames
        IgPoolList(ind,jnd) = ITotalgPool(knd)
        RGPoolSg(ind,jnd) = RTotalSgPool(knd)
      END DO
    END DO
    ! Clean up
    DEALLOCATE(ILocalgPool,ITotalgPool,RLocalSgPool,RTotalSgPool,Inum,Ipos)

    !-2------------------------------------------------------------------
    ! Fill IgOutList & output as a text file
    IF(my_rank.EQ.0) THEN
      CALL message(LS,dbg3,"Writing hkl list and images")
      path = SChemicalFormula(1:ILN) // "/hkl_list.txt"
      OPEN(UNIT=IChOutIhkl, ACTION='WRITE', POSITION='APPEND', STATUS= 'UNKNOWN', &
          FILE=TRIM(ADJUSTL(path)),IOSTAT=IErr)
      WRITE(IChOutIhkl,*) "List of hkl in each frame"
      WRITE(IChOutIhkl,*) "No: h k l  Fg  |g|  Sg"
    END IF
    IgOutList = 0
    IMaxNg = 0  ! maximum number of outputs in a frame
    DO ind = 1,INFrames
      ! output to slurm if requested
      CALL message(LM, "Reflection list:")
      DO knd = 1, INhkl
        IF (IgPoolList(knd,ind).NE.0) THEN
          WRITE(SPrintString,'(I3,1X,I3,1X,I3)') IhklLattice(IgPoolList(knd,ind),:)
          CALL message(LM, SPrintString)
        END IF
      END DO
      ! write reflections in each frame
      ! h k l Fg(Re Im) |g| Sg
      IF (my_rank.EQ.0) WRITE(IChOutIhkl,"(A6,I4)") "Frame ",ind
      jnd = 0  ! to count reflections per frame
      DO knd = 1, INhkl
        IF (IgPoolList(knd,ind).NE.0) THEN
          lnd = IgPoolList(knd,ind)
          RSg = RgPoolSg(knd,ind)
          ! Is this reflection small enough to be in the output list
          IF (RgMagLattice(IgPoolList(knd,ind)).LT.RGOutLimit) THEN
          jnd = jnd + 1
            IgOutList(knd,ind) = IgPoolList(knd,ind)
          !END IF  !** remove comment to give all reflections, complements ##
            WRITE(fString,"(3(I3,1X),2X, F8.4,A1,F8.4,A3, F6.2,2X, F8.4)") &
                  IhklLattice(lnd,:), REAL(CFgLattice(lnd)),"+",AIMAG(CFgLattice(lnd)),"i  ",&
                  RgMagLattice(lnd)/TWOPI, RSg
            IF (my_rank.EQ.0) WRITE(IChOutIhkl,*) TRIM(ADJUSTL(fString))
          END IF  !## remove comment to give output reflections, complements ** 
        END IF
      END DO
      WRITE(fString,"(A6,I3,A19)") "Total ",jnd," output reflections"
      IF (my_rank.EQ.0) WRITE(IChOutIhkl,*) TRIM(ADJUSTL(fString))
      IF (jnd.GT.IMaxNg) IMaxNg = jnd  ! update max number of outputs if necessary
    END DO
    IF (my_rank.EQ.0) CLOSE(IChOutIhkl,IOSTAT=IErr)

    !-3------------------------------------------------------------------
    ! Write a set of kinematic simulation frames
    ! image
    ISim = 256_IKIND  ! NB HALF the output image size = output |g| limit/0.98
    ALLOCATE(RSim(2*ISim,2*ISim),STAT=IErr)
    IF(l_alert(IErr,"HKLmake","allocate RSim")) RETURN
    ! Instrument broadening term - sets the FWHM  of a kinematic rocking curve
    RInst = 3000.0
    DO ind = 1,INFrames
      RAngle = REAL(ind-1)*DEG2RADIAN*RFrameAngle
      RSim = ZERO
      ! Direct beam
      RSim(ISim-1:ISim+1,ISim-1:ISim+1) = 1
      ! output g's
      DO knd = 1, INhkl
        IF (IgOutList(knd,ind).NE.0) THEN
          lnd = IgPoolList(knd,ind)  ! index of reflection in the reciprocal lattice
          Rg = RgLatticeO(lnd,:)  ! g-vector
          RIkin = CFgLattice(lnd)*CONJG(CFgLattice(lnd))  ! simple kinematic intensity
          RSg = RgPoolSg(knd,ind)  !Sg
          ! x- and y-coords (NB swapped in the image!)
          Rp = RXDirO*COS(RAngle)-RZDirO*SIN(RAngle)  ! unit vector horizontal in the image
          ! position of the spot, 2% leeway to avoid going over the edge of the image
          Ix = ISim-0.98*NINT(DOT_PRODUCT(Rg,Rp)*REAL(ISim)/RGOutLimit)  
          Iy = ISim+0.98*NINT(DOT_PRODUCT(Rg,RYDirO)*REAL(ISim)/RGOutLimit)
          RSim(Iy-1:Iy+1,Ix-1:Ix+1) = EXP(-RInst*RSg*RSg)*RIkin
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
      IF(my_rank.EQ.0) THEN
        OPEN(UNIT=IChOutIM, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(path)),&
          FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=2*ISim*IByteSize)
        IF(l_alert(IErr,"WriteIterationOutput","OPEN() output .bin file")) RETURN      
        DO jnd = 1,2*ISim
          WRITE(IChOutIM,rec=jnd) RSim(jnd,:)
        END DO
        CLOSE(IChOutIM,IOSTAT=IErr) 
        IF(l_alert(IErr,"WriteIterationOutput","CLOSE() output .bin file")) RETURN
      END IF
    END DO
    DEALLOCATE(RSim)
    

  END SUBROUTINE HKLmake

  !!$%%HKLList%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: List the frames for each reflection
  !!
  !! Major-Authors: Richard Beanland (2023)
  !!  
  SUBROUTINE HKLList( IErr )

    ! This procedure is called once in felixrefine setup
    ! 1) get unique g's in all pools
    ! 2) make a new list of unique g's and associated parameters and delete reciprocal lattice
    ! 3) save a set of kinematic rocking curves
    
    USE MyNumbers
    USE message_mod

    ! global parameters
    USE IPARA, ONLY : ILN,INFrames,ISort,INhkl,IgPoolList,IgOutList,Ihkl,IhklLattice,INoOfHKLsAll
    USE SPARA, ONLY : SPrintString,SChemicalFormula
    USE RPARA, ONLY : RgO,RgMag,RgLatticeO,RgMagLattice,RgPoolSg
    USE CPARA, ONLY : CFg,CFgLattice
    USE IChannels, ONLY : IChOutIhkl

    IMPLICIT NONE

    REAL(RKIND) :: RInst,RIkin
    INTEGER(IKIND) :: ind,jnd,knd,lnd,Iy,IErr,Imin,Imax
    INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: IFullgList,IReducedgList,IUniquegList
    CHARACTER(200) :: path
    CHARACTER(100) :: fString

    !-1------------------------------------------------------------------
    ! first get a list of pool reflections without duplicates, IUniquegList
    ind = INFrames*INhkl
    ALLOCATE(IFullgList(ind), STAT=IErr)  ! everything
    ALLOCATE(IReducedgList(ind), STAT=IErr)  ! unique reflections in an oversize matrix
    IF(l_alert(IErr,"HKLlist","allocations 1")) RETURN
    IFullgList = RESHAPE(IgPoolList,[ind])
    !now the list of unique reflections
    Imin = 0
    Imax = MAXVAL(IFullgList)
    ind = 0
    DO WHILE (Imin.LT.Imax)
        ind = ind+1
        Imin = MINVAL(IFullgList, MASK=IFullgList.GT.Imin)
        IReducedgList(ind) = Imin
    END DO
    ALLOCATE(IUniquegList(ind), STAT=IErr)
    IF(l_alert(IErr,"HKLlist","allocations 2")) RETURN
    IUniquegList = IReducedgList(1:ind)
    IF(l_alert(IErr,"HKLlist","allocate IUniquegList")) RETURN
    ! Tidy up
    DEALLOCATE(IFullgList,IReducedgList)
    INoOfHKLsAll = ind
    WRITE(SPrintString, FMT='(I5,A19)') ind, " pool reflections"
    CALL message(LS,SPrintString)

    !-2------------------------------------------------------------------
    ! Make reduced lists of hkl, g-vector, |g| and Fg so we can deallocate the reciprocal lattice
    ALLOCATE(Ihkl(INoOfHKLsAll,ITHREE), STAT=IErr)  ! Miller indices
    ALLOCATE(RgO(INoOfHKLsAll,ITHREE), STAT=IErr)  ! g-vector, orthogonal frame
    ALLOCATE(RgMag(INoOfHKLsAll), STAT=IErr)  ! |g|
    ALLOCATE(CFg(INoOfHKLsAll), STAT=IErr)  ! Fg
    IF(l_alert(IErr,"HKLlist","allocations 3")) RETURN
    DO jnd = 1,INoOfHKLsAll
      Ihkl(jnd,:) = IhklLattice(IUniquegList(jnd),:)
      RgO(jnd,:) = RgLatticeO(IUniquegList(jnd),:)
      RgMag(jnd) = RgMagLattice(IUniquegList(jnd))
      CFg(jnd) = CFgLattice(IUniquegList(jnd))
    END DO
    ! Change the indices for IgPoolList and IgOutList
    ! For a frame [j] and a given reflection in the beam pool [i,j],
    ! we find hkl, g, |g| and Fg at the index given in IgPoolList[i,j].
    ! IgOutList[i,j] gives the number of the output reflection.
    lnd = 0  ! counter for output reflections
    DO ind = 1,INoOfHKLsAll
      Iy = 1  ! flag for counting
      DO jnd = 1, INhkl
        DO knd = 1, INFrames
          IF (IgPoolList(jnd,knd).EQ.IUniquegList(ind)) IgPoolList(jnd,knd) = ind
          IF (IgOutList(jnd,knd).EQ.IUniquegList(ind)) THEN
            IF (Iy.EQ.1) THEN  ! we only count the first appearance
              lnd = lnd + 1
              Iy = 0
            END IF
            IgOutList(jnd,knd) = lnd
          END IF
        END DO
      END DO
    END DO
    WRITE(SPrintString, FMT='(I5,A19)') lnd, " output reflections"
    CALL message(LS,SPrintString)

    !-3------------------------------------------------------------------
    ! kinematic rocking curves  
    RInst = 3000.0  ! instrument broadening term
    IF(my_rank.EQ.0) THEN
      CALL message(LS,dbg3,"Writing kinematic rocking curves")
      path = SChemicalFormula(1:ILN) // "/hkl_K-rocks.txt"
      OPEN(UNIT=IChOutIhkl, ACTION='WRITE', POSITION='APPEND', STATUS= 'UNKNOWN', &
          FILE=TRIM(ADJUSTL(path)),IOSTAT=IErr)
      WRITE(IChOutIhkl,*) "List of kinematic rocking curves"
      DO ind = 1,lnd
        Iy = 1
        DO knd = 1,INFrames
          DO jnd = 1,INhkl
            IF (IgOutList(jnd,knd).EQ.ind) THEN
              IF (Iy.EQ.1) THEN
                WRITE(fString,"(3(I3,1X))") Ihkl(IgPoolList(jnd,knd),:)
                WRITE(IChOutIhkl,*) TRIM(ADJUSTL(fString))
                Iy = 0
              END IF
              RIkin = CFg(IgPoolList(jnd,knd))*CONJG(CFg(IgPoolList(jnd,knd))) * &
                  EXP(-RInst*RgPoolSg(jnd,knd)*RgPoolSg(jnd,knd)) ! Gaussian shape of reflection with Sg
              WRITE(fString,"(I4,A3,F7.3)") knd," : ",RIkin
              WRITE(IChOutIhkl,*) TRIM(ADJUSTL(fString))
            END IF
          END DO
        END DO    
      END DO
      CLOSE(IChOutIhkl,IOSTAT=IErr)
    END IF
    
 
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

