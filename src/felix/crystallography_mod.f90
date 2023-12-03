!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>
!! Module-description: This defines lattice vectors as well as the fractional atomic coordinates
!!
MODULE crystallography_mod
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: ReciprocalVectors, HKLSave, CrystalOrientation, UniqueAtomPositions, gVectors, HKLMake, HKLPlot

  CONTAINS

  !>
  !! Procedure-description: Calculates g-vector matrices, global variables
  !!
  !! Author:  r.beanland@warwick.ac.uk
  !!
  SUBROUTINE gVectors(IErr)

    USE MyNumbers
    USE message_mod
    USE MyMPI 
 
    ! global inputs
    USE RPARA, ONLY : RarVecM,RbrVecM,RcrVecM,RNormDirM,Rhkl,RConvergenceAngle
    USE IPARA, ONLY : INhkl
    USE SPARA, ONLY : SPrintString    
    ! global outputs
    USE RPARA, ONLY : RgPool,RgPoolMag,RgDotNorm,RDeltaK,RgMatrix
    
    IMPLICIT NONE
    
    INTEGER(IKIND) :: IErr,ind,jnd

    IErr=0!No route to throw an error here in fact
    
    !calculate g-vector pool, the magnitudes and component parallel to specimen surface
    DO ind=1,INhkl
      DO jnd=1,ITHREE
        RgPool(ind,jnd) = Rhkl(ind,1)*RarVecM(jnd) + &
            Rhkl(ind,2)*RbrVecM(jnd) + Rhkl(ind,3)*RcrVecM(jnd)
      END DO
      RgPoolMag(ind)= SQRT(DOT_PRODUCT(RgPool(ind,:),RgPool(ind,:)))
      RgDotNorm(ind) = DOT_PRODUCT(RgPool(ind,:),RNormDirM)
    END DO
    
    ! Calculate matrix  of g-vectors that corresponds to the Ug matrix
    DO ind=1,INhkl
      DO jnd=1,INhkl
        RgMatrix(ind,jnd,:)= RgPool(ind,:)-RgPool(jnd,:)
      END DO
    END DO

    !outputs if requested    
    CALL message(LL,dbg3,"first 16 g-vectors", RgMatrix(1:16,1,:)) 
    CALL message(LL,dbg7,"g-vectors and magnitude (1/A), in the microscope reference frame" )
    DO ind = 1,INhkl
      CALL message(LL,dbg7,"hkl  :",NINT(Rhkl(ind,:)))
      CALL message(LL,dbg7,"g mag:",RgPoolMag(ind))
    END DO
    CALL message(LL,dbg7,"g.n list")
    DO ind = 1,INhkl
      CALL message(LL,dbg7,"hkl :",NINT(Rhkl(ind,:)))
      CALL message(LL,dbg7,"g.n :",RgDotNorm(ind))
    END DO
    
    END SUBROUTINE gVectors


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Creates basis reciprocal lattice vectors and initial
  !! microscope reference frame, gives the number of reciprocal lattice points InLattice
  !! that will be calculated later in the ReciprocalLattice subroutine
  !!
  !! Major-Authors: Richard Beanland (2023)
  !!
  SUBROUTINE ReciprocalVectors(IErr)

    USE MyNumbers
    USE MyMPI
    USE message_mod

    ! global inputs
    USE IPARA, ONLY : IVolumeFLAG,INAtomsUnitCell
    USE RPARA, ONLY : RAlpha,RBeta,RGamma,RCellA,RCellB,RCellC,RXDirC_0,RZDirC_0,&
            RAtomCoordinate,RAtomXYZ
    USE SPARA, ONLY : SSpaceGroupName,SPrintString

    ! global outputs
    USE RPARA, ONLY : RaVecO,RbVecO,RcVecO,RVolume,RarVecO,RbrVecO,RcrVecO,&
            RXDirO,RYDirO,RZDirO,RarMag,RbrMag,RcrMag
!    USE IPARA, ONLY : inda,indb,indc

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr,ind
    REAL(RKIND) :: Rt,RxAngle
    REAL(RKIND), DIMENSION(ITHREE,ITHREE) :: RTMatC2O
    
    !direct lattice vectors in an orthogonal reference frame, Angstrom units 
    ! a is parallel to [100]
    RaVecO(1)= RCellA
    RaVecO(2)= ZERO
    RaVecO(3)= ZERO
    ! b lies in the x-y plane 
    RbVecO(1)= RCellB*COS(RGamma)
    RbVecO(2)= RCellB*SIN(RGamma)
    RbVecO(3)= ZERO
    ! c lies... wherever
    RcVecO(1)= RCellC*COS(RBeta)
    RcVecO(2)= RCellC*(COS(RAlpha)-COS(RBeta)*COS(RGamma))/SIN(RGamma)
    RcVecO(3)= RCellC*(SQRT(1.D0-COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)&
      -COS(RGamma)*COS(RGamma)+TWO*COS(RAlpha)*COS(RBeta)*COS(RGamma)) / SIN(RGamma))

    ! RTmatC2O transforms from crystal (implicit units) 
    ! to orthogonal O reference frame (Angstrom units)
    RTMatC2O(:,1) = RaVecO(:)
    RTMatC2O(:,2) = RbVecO(:)
    RTMatC2O(:,3) = RcVecO(:)

    ! Atom coordinates in the orthogonal reference frame
    DO ind=1,INAtomsUnitCell
      RAtomCoordinate(ind,:) = MATMUL(RTMatC2O,RAtomXYZ(ind,:))
    END DO

    !calculate cell volume if required
    IF(IVolumeFLAG .EQ. 0) THEN
       RVolume= RCellA*RCellB*RCellC* &
            SQRT(1.0D0 - &
            COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)-COS(RGamma)*COS(RGamma) + &
            TWO*COS(RAlpha)*COS(RBeta)*COS(RGamma))
    END IF

    !Some checks for rhombohedral cells?
    Rt = DOT_PRODUCT(RaVecO/DOT_PRODUCT(RaVecO,RaVecO),RbVecO/DOT_PRODUCT(RbVecO,RbVecO))*&
         DOT_PRODUCT(RbVecO/DOT_PRODUCT(RbVecO,RbVecO),RcVecO/DOT_PRODUCT(RcVecO,RcVecO))*&
         DOT_PRODUCT(RcVecO/DOT_PRODUCT(RcVecO,RcVecO),RaVecO/DOT_PRODUCT(RaVecO,RaVecO))
       
    IF(SCAN(SSpaceGroupName,'rR').NE.0) THEN
       IF(ABS(Rt).LT.TINY) THEN
        SSpaceGroupName = TRIM(ADJUSTL("V"))
        ! Crystal is either Obverse or Reverse
        ! Selection Rules are not in place to determine the difference, 
        ! assume the crystal is Obverse")
      ELSE
        SSpaceGroupName=TRIM(ADJUSTL('P'))
        ! Primitive setting (Rhombohedral axes)
      END IF
    END IF

    ! Set up Reciprocal Lattice Vectors: orthogonal reference frame in 1/Angstrom units
    ! RarDirO,RbrDirO,RcrDirO vectors are reciprocal lattice vectors 
    ! 2pi/a, 2pi/b, 2pi/c in an orthogonal frame
    ! Note that reciprocal lattice vectors have two pi included,
    ! we are using the physics convention exp(i*g.r)
    RarVecO= TWOPI*CROSS(RbVecO,RcVecO)/DOT_PRODUCT(RbVecO,CROSS(RcVecO,RaVecO))
    RbrVecO= TWOPI*CROSS(RcVecO,RaVecO)/DOT_PRODUCT(RcVecO,CROSS(RaVecO,RbVecO))
    RcrVecO= TWOPI*CROSS(RaVecO,RbVecO)/DOT_PRODUCT(RaVecO,CROSS(RbVecO,RcVecO))
    DO ind=1,ITHREE
       IF (abs(RarVecO(ind)).LT.TINY) THEN
          RarVecO(ind) = ZERO
       END IF
       IF (abs(RbrVecO(ind)).LT.TINY) THEN
          RbrVecO(ind) = ZERO
       END IF
       IF (abs(RcrVecO(ind)).LT.TINY) THEN
          RcrVecO(ind) = ZERO
       END IF
    ENDDO

    ! Set up initial microscope reference frame
    ! X, Y and Z are orthogonal vectors that defines the simulation
    ! Also referred to as the microscope reference frame M.
    ! The electron beam propagates along +Zm.
    ! The alpha rotation axis is along Ym.  Positive alpha rotation moves the field
    ! of view of the simulation along +Xm.
    ! In the crystal reference frame we read in reciprocal vectors RXDirC_0 and RZDirC_0
    ! These define Xm & Zm in the inital reference frame
    ! RXDirO,RYDirO,RZDirO are UNIT reciprocal lattice vectors parallel to X,Y,Z
    RXDirO = RXDirC_0(1)*RarVecO + RXDirC_0(2)*RbrVecO + RXDirC_0(3)*RcrVecO
    RXDirO = RXDirO/SQRT(DOT_PRODUCT(RXDirO,RXDirO))
    RZDirO = RZDirC_0(1)*RaVecO + RZDirC_0(2)*RbVecO + RZDirC_0(3)*RcVecO
    RZDirO = RZDirO/SQRT(DOT_PRODUCT(RZDirO,RZDirO))
    ! Check the input is sensible, i.e. Xm is perpendicular to Zm
    RxAngle = ABS(180.0D0*ACOS(DOT_PRODUCT(RXDirO,RZDirO))/PI)
    IF(ABS(RxAngle-90.0D0).GT.0.1)THEN! with a tolerance of 0.1 degrees
      WRITE(SPrintString,"(A15,F5.1,A27)") "Error: X is at ",RxAngle," degrees to Z, should be 90"
      CALL message(LS,SPrintString)
      IErr = 1
    ELSE!fine correction of x
      ! take off any component parallel to z & renormalise
      RXDirO = RXDirO - DOT_PRODUCT(RXDirO,RZDirO)*RZDirO
      RXDirO = RXDirO/SQRT(DOT_PRODUCT(RXDirO,RXDirO))
    END IF
    RYDirO = CROSS(RZDirO,RXDirO)  ! the rotation axis
    
    ! Reciprocal vector magnitudes and size of the lattice
    RarMag=SQRT(DOT_PRODUCT(RarVecO,RarVecO))!magnitude of a*
    RbrMag=SQRT(DOT_PRODUCT(RbrVecO,RbrVecO))!magnitude of b*
    RcrMag=SQRT(DOT_PRODUCT(RcrVecO,RcrVecO))!magnitude of c*

  END SUBROUTINE ReciprocalVectors


  !!$%%HKLMake%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description:
  !! Fills the beam pool list RgPoolList for each frame (global variable)
  !!
  !! Major-Authors: Richard Beanland (2023)
  !!  
  SUBROUTINE HKLMake(RDevLimit, RGOutLimit, RgPoolLimit, IErr)   

    USE MyNumbers
    USE message_mod
    USE myMPI

    ! global inputs/outputs
    USE SPARA, ONLY : SPrintString
    USE IPARA, ONLY : INhkl,ILN,INFrames  ! inputs
    USE IPARA, ONLY : Ig,IgOutList,IgPoolList  ! outputs
    USE RPARA, ONLY : RXDirO,RYDirO,RZDirO,RarVecO,RbrVecO,RcrVecO,RarMag,RbrMag,RcrMag,RFrameAngle,RBigK,&
          RgPoolSg ! only RgPoolSg is an output
    USE Iconst
    
    IMPLICIT NONE

    REAL(RKIND),INTENT(IN) :: RDevLimit, RGOutLimit, RgPoolLimit
    INTEGER(IKIND) :: IErr,ind,jnd,knd,lnd,mnd,nnd,ond,ISel,Ifull(INFrames),IMaxNg,inda,indb,indc,Ifound
    REAL(RKIND) :: RAngle,Rk(ITHREE),Rk0(ITHREE),Rp(ITHREE),RSg,Rphi,Rg(ITHREE),RIkin,&
                   RKplusg(ITHREE),RgMag,RShell,RgMin
   
    !-1------------------------------------------------------------------
    ! calculate reflection list g by g
    !--------------------------------------------------------------------
    ! this produces a list of g-vectors Ig, their deviation parameters RgPoolSg, and fills IgPoolList, IgOutList
    ! Initialise variables
    Ig = 0  ! list of reflections (covers all beam pools)
    IgPoolList = 0  ! index giving Ig for each beam pool (1 to INhkl,1 to INFrames)
    IgOutList = 0  ! index giving Ig for each output (1 to INhkl,1 to INFrames)
    IFull = 0  ! is the beam pool full
    RgPoolSg = TEN  !Sg corresponding to Ig (1 to INhkl,1 to INFrames)
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
              RAngle = REAL(nnd-1)*DEG2RADIAN*RFrameAngle
              Rk = RBigK*(RZDirO*COS(RAngle)+RXDirO*SIN(RAngle))
              ! Calculate Sg by getting the vector k0, which is coplanar with k and g and
              ! corresponds to an incident beam at the Bragg condition
              ! First we need the vector component of k perpendicular to g, which we call p 
              Rp = Rk - DOT_PRODUCT(Rk,Rg)*Rg/(RgMag**2)
              ! and now make k0 by adding vectors parallel to g and p
              ! i.e. k0 = (p/|p|)*(k^2-g^2/4)^0.5 - g/2
              Rk0 = SQRT(RBigK**2-QUARTER*RgMag**2)*Rp/SQRT(DOT_PRODUCT(Rp,Rp)) - HALF*Rg
              ! The angle phi between k and k0 is how far we are from the Bragg condition
              Rphi = ACOS(DOT_PRODUCT(Rk,Rk0)/(RBigK**2))
              ! and now Sg is 2g sin(phi/2), with the sign of K-|K+g|
              RKplusg = Rk + Rg
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
    IF (lnd.LT.INhkl.AND.RgMin.LE.RgPoolLimit) THEN
      mnd = mnd + 1
      GOTO 1
    END IF

  END SUBROUTINE HKLmake


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: writes hkl_list.txt, the list of reflections in each frame
  !!
  !! Major-Authors: Richard Beanland (2023)
  !!
  SUBROUTINE HKLSave(RgOutLimit,IErr)

    USE MyNumbers
    USE MyMPI
    USE message_mod
    USE ug_matrix_mod

    ! global inputs
    USE IPARA, ONLY : ILN,INFrames,INhkl,Ig,IGPoolList,IgOutList,ICurrentZ,INAtomsUnitCell,&
            IAtomicNumber
    USE RPARA, ONLY : RarVecO,RbrVecO,RcrVecO,RgPoolSg,RAtomCoordinate,RIsoDW
    USE SPARA, ONLY : SPrintString,SChemicalFormula
    USE IChannels, ONLY : IChOutIhkl

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr,ind,jnd,knd,lnd,mnd,ISim,Ix,Iy
    REAL(RKIND) :: Rg(ITHREE),RgMag,RSg,Rfq,RSim(512,512),RAngle,RInst,Rp(ITHREE),RIkin
    REAL(RKIND), INTENT(IN) :: RgOutLimit 
    COMPLEX :: CFg
    CHARACTER(200) :: path
    CHARACTER(100) :: fString

    !-1------------------------------------------------------------------
    ! write reflection lists and rocking curves to hkl_list
    IF(my_rank.EQ.0) THEN
      CALL message(LS,dbg3,"Writing hkl list and images")
      path = SChemicalFormula(1:ILN) // "/hkl_list.txt"
      OPEN(UNIT=IChOutIhkl, ACTION='WRITE', POSITION='APPEND', STATUS= 'UNKNOWN', &
          FILE=TRIM(ADJUSTL(path)),IOSTAT=IErr)
      WRITE(IChOutIhkl,*) "List of hkl in each frame"
      WRITE(IChOutIhkl,*) "No: h k l  Fg  |g|  Sg"
    END IF
    DO ind = 1,INFrames
      ! output to slurm if requested
      CALL message(LM, "Reflection list:")
      DO knd = 1, INhkl
        IF (IgPoolList(knd,ind).NE.0) THEN
          WRITE(SPrintString,'(I3,1X,I3,1X,I3)') Ig(IgPoolList(knd,ind),:)
          CALL message(LM, SPrintString)
        END IF
      END DO
      ! write reflections in each frame
      ! h k l Fg(Re Im) |g| Sg
      IF (my_rank.EQ.0) WRITE(IChOutIhkl,"(A6,I4)") "Frame ",ind
      DO knd = 1, INhkl
        IF (IgOutList(knd,ind).NE.0) THEN
          lnd = IgOutList(knd,ind)
          RSg = RgPoolSg(knd,ind)
          Rg = Ig(lnd,1)*RarVecO + Ig(lnd,2)*RbrVecO + Ig(lnd,3)*RcrVecO
          RgMag = SQRT(DOT_PRODUCT(Rg,Rg))
          ! Calculate structure factor
          CFg = CZERO
          DO mnd=1,INAtomsUnitCell
            ICurrentZ = IAtomicNumber(mnd)
            CALL AtomicScatteringFactor(Rfq,IErr)  ! in ug_matrix_mod
            CFg = CFg+Rfq*EXP(-CIMAGONE*DOT_PRODUCT(Rg,RAtomCoordinate(mnd,:)) ) * &
            ! Isotropic D-W factor exp(-B sin(theta)^2/lamda^2) = exp(-Bg^2/16pi^2)
            EXP(-RIsoDW(mnd)*RgMag**2/(FOURPI**2))
          END DO

          WRITE(fString,"(3(I3,1X),2X, F8.4,A1,F8.4,A3, F6.2,2X, F8.4)") &
                  Ig(lnd,:), REAL(CFg),"+",AIMAG(CFg),"i  ",&
                  RgMag/TWOPI, RSg
          IF (my_rank.EQ.0) WRITE(IChOutIhkl,*) TRIM(ADJUSTL(fString))
        END IF
      END DO
!      IF (jnd.GT.IMaxNg) IMaxNg = jnd  ! update max number of outputs if necessary
    END DO

    ! write frames for each reflexion
    DO ind = 1,INhkl*INFrames
      Rg = Ig(ind,:)
      IF (DOT_PRODUCT(Rg,Rg).GT.TINY) THEN  ! this reflexion is not zero
        IF (my_rank.EQ.0) WRITE(IChOutIhkl,"(A10,I4,2X,3(I3,1X))") "Reflexion ",ind,Ig(ind,:)
        DO jnd = 1,INFrames
          DO knd = 1,INhkl
            IF (IgOutList(knd,jnd).EQ.ind) THEN
              RSg = RgPoolSg(knd,jnd)
              IF (my_rank.EQ.0) WRITE(IChOutIhkl,"(I4 ,2X, F8.4)") jnd,RSg
            END IF
          END DO
        END DO
      END IF
    END DO 
    IF (my_rank.EQ.0) CLOSE(IChOutIhkl,IOSTAT=IErr)

END SUBROUTINE HKLSave


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: writes hkl_list.txt, the list of reflections in each frame
  !!
  !! Major-Authors: Richard Beanland (2023)
  !!
  SUBROUTINE HKLPlot(ISim, RgOutLimit,IErr)

    USE MyNumbers
    USE MyMPI
    USE message_mod
    USE ug_matrix_mod

    ! global inputs
    USE IPARA, ONLY : ILN,INFrames,INhkl,Ig,IgPoolList,IgOutList,ICurrentZ,INAtomsUnitCell,&
            IAtomicNumber,IByteSize
    USE RPARA, ONLY : RarVecO,RbrVecO,RcrVecO,RgPoolSg,RAtomCoordinate,RIsoDW,RxDirO,RyDirO,RzDirO,RFrameAngle
    USE SPARA, ONLY : SPrintString,SChemicalFormula
    USE IChannels, ONLY : IChOutIM

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr,ind,jnd,knd,lnd,mnd,Ix,Iy
    INTEGER(IKIND), INTENT(IN) :: ISim
    REAL(RKIND) :: Rg(ITHREE),RgMag,RSg,Rfq,RSim(2*ISim,2*ISim),RAngle,RInst,Rp(ITHREE),RIkin,RmaxI
    REAL(RKIND), INTENT(IN) :: RgOutLimit 
    COMPLEX :: CFg
    CHARACTER(200) :: path
    CHARACTER(100) :: fString

    !--------------------------------------------------------------------
    ! Write a set of kinematic simulation frames
    ! Instrument broadening term - sets the FWHM  of a kinematic rocking curve
    RInst = 3000.0
    DO ind = 1,INFrames
      RAngle = REAL(ind-1)*DEG2RADIAN*RFrameAngle
      RSim = ZERO
      RmaxI = 0.0  ! max g in the image
      ! output g's
      DO knd = 2, INhkl
        IF (IgOutList(knd,ind).NE.0) THEN
          lnd = IgPoolList(knd,ind)  ! index of reflection in the reciprocal lattice
!DBG      IF(my_rank.EQ.0)PRINT*,ind,":",Ig(lnd,:)
          Rg = Ig(lnd,1)*RarVecO + Ig(lnd,2)*RbrVecO + Ig(lnd,3)*RcrVecO  ! g-vector
          RgMag = SQRT(DOT_PRODUCT(Rg,Rg))
          ! Calculate structure factor
          CFg = CZERO
          DO mnd=1,INAtomsUnitCell
            ICurrentZ = IAtomicNumber(mnd)
            CALL AtomicScatteringFactor(Rfq,IErr)  ! in ug_matrix_mod
            CFg = CFg+Rfq*EXP(-CIMAGONE*DOT_PRODUCT(Rg,RAtomCoordinate(mnd,:)) ) * &
            ! Isotropic D-W factor exp(-B sin(theta)^2/lamda^2) = exp(-Bg^2/16pi^2)
            EXP(-RIsoDW(mnd)*RgMag**2/(FOURPI**2))
          END DO
          RIkin = CFg*CONJG(CFg)  ! simple kinematic intensity
          !IF(RIkin.GT.RmaxI) RmaxI = RIkin
          RSg = RgPoolSg(knd,ind)  !Sg
          ! x- and y-coords (NB swapped in the image!)
          Rp = RXDirO*COS(RAngle)-RZDirO*SIN(RAngle)  ! unit vector horizontal in the image
          ! position of the spot, 2% leeway to avoid going over the edge of the image
          Ix = ISim-0.98*NINT(DOT_PRODUCT(Rg,Rp)*REAL(ISim)/RGOutLimit)  
          Iy = ISim+0.98*NINT(DOT_PRODUCT(Rg,RYDirO)*REAL(ISim)/RGOutLimit)
          RSim(Iy-1:Iy+1,Ix-1:Ix+1) = RIkin*EXP(-RInst*RSg*RSg)
        END IF
      END DO
      ! direct beam
      RSim(ISim-1:ISim+1,ISim-1:ISim+1) = 100.0!RmaxI
      ! write to disk - set up file name
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
      IF(my_rank.EQ.0) THEN  ! open file to write
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
  
  END SUBROUTINE HKLPlot
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE CrystalOrientation(IErr)
    USE MyNumbers
    USE MyMPI
    USE message_mod

    USE RPARA, ONLY : RarVecM,RbrVecM,RcrVecM,RaVecM,RbVecM,RcVecM,RNormDirM,RaVecO,RbVecO,&
          RcVecO,RarVecO,RbrVecO,RcrVecO,RXDirO,RYDirO,RZDirO,RarMag,RbrMag,RcrMag
    USE SPARA, ONLY : SSpaceGroupName

    ! global inputs
    USE IPARA, ONLY : IVolumeFLAG,INAtomsUnitCell
    USE RPARA, ONLY : RAlpha,RBeta,RGamma,RCellA,RCellB,RCellC,RNormDirC,RXDirC,&
          RZDirC,RXDirM,RAtomCoordinate,RAtomXYZ
    USE SPARA, ONLY : SPrintString

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr,ind,jnd
    REAL(RKIND), DIMENSION(ITHREE,ITHREE) :: RTMatC2O,RTMatO2M
    REAL(RKIND) :: RNormAngle

    ! RTmatC2O transforms from crystal (implicit units) 
    ! to orthogonal reference frame (Angstrom units)
    RTMatC2O(:,1) = RaVecO(:)
    RTMatC2O(:,2) = RbVecO(:)
    RTMatC2O(:,3) = RcVecO(:)

    ! RTmatO2M transforms from orthogonal to microscope reference frame
    RTMatO2M(1,:) = RXDirO(:)
    RTMatO2M(2,:) = RYDirO(:)
    RTMatO2M(3,:) = RZDirO(:)

    ! Unit normal to the specimen in REAL space
    ! This is used for all g-vectors as a boundary condition
    RNormDirM = MATMUL(RTMatO2M,MATMUL(RTMatC2O,RNormDirC))
    RNormDirM = RNormDirM/SQRT(DOT_PRODUCT(RNormDirM,RNormDirM)) 

    ! Check to see if the normal is close to perpendicular to the incident beam
    ! remember in the microscope frame z = [001] so n.z is just n(3)
    RNormAngle = ABS(180.0D0*ACOS(RNormDirM(3))/PI)
    IF(RNormAngle.GT.75.0D0 .AND. my_rank.EQ.0) THEN
      WRITE(SPrintString,"(A30,F5.1,A8)") "Warning: surface normal is at ",RNormAngle," degrees"
      CALL message(LS,SPrintString) 
    END IF

    ! Transform atomic coordinates to microscope frame to give RAtomCoordinate
    ! RaVecM, RbVecM, RbVecM unit cell vectors in Angstrom units in the microscope frame
    RaVecM= MATMUL(RTMatO2M,RaVecO)
    RbVecM= MATMUL(RTMatO2M,RbVecO)
    RcVecM= MATMUL(RTMatO2M,RcVecO)
    ! Calculate atomic position vectors RAtomCoordinate
    ! In microscope reference frame, in Angstrom units (NB RAtomXYZ=crystal frame, in .cif)
    DO ind=1,INAtomsUnitCell
      DO jnd=1,ITHREE
        RAtomCoordinate(ind,jnd)= RAtomXYZ(ind,1)*RaVecM(jnd) + &
              RAtomXYZ(ind,2)*RbVecM(jnd)+RAtomXYZ(ind,3)*RcVecM(jnd)
      END DO
    END DO
    
    ! create reciprocal lattice vectors in Microscope reference frame
    ! Note that reciprocal lattice vectors have two pi included,
    ! we are using the optical convention exp(i*g.r)
    RarVecM= TWOPI*CROSS(RbVecM,RcVecM)/DOT_PRODUCT(RbVecM,CROSS(RcVecM,RaVecM))  ! a*
    RbrVecM= TWOPI*CROSS(RcVecM,RaVecM)/DOT_PRODUCT(RcVecM,CROSS(RaVecM,RbVecM))  ! b*
    RcrVecM= TWOPI*CROSS(RaVecM,RbVecM)/DOT_PRODUCT(RaVecM,CROSS(RbVecM,RcVecM))  ! c*
    ! their magnitudes
    RarMag=SQRT(DOT_PRODUCT(RarVecM,RarVecM))!magnitude of a*
    RbrMag=SQRT(DOT_PRODUCT(RbrVecM,RbrVecM))!magnitude of b*
    RcrMag=SQRT(DOT_PRODUCT(RcrVecM,RcrVecM))!magnitude of c*

  END SUBROUTINE CrystalOrientation

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Calculates the full set of possible fractional atomic positions,
  !! the mean inner potential and wavevector in the material K
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016, 2017)
  !!
  SUBROUTINE UniqueAtomPositions(IErr)

    ! updates full crystal arrays from basis atom refinement
    ! NB Anisotropic DWFs yet to be properly implemented

    USE MyNumbers
    USE message_mod
    USE ug_matrix_mod

    ! global outputs
    USE RPARA, ONLY : RAtomCoordinate,ROccupancy,RIsoDW,RAtomXYZ,RMeanInnerPotential,RBigK
    USE IPARA, ONLY : IAtomicNumber,IAnisoDW
    USE SPARA, ONLY : SAtomLabel, SAtomName

    ! global inputs
    USE RPARA, ONLY : RBasisOccupancy,RBasisIsoDW,RSymVec,RBasisAtomPosition,RSymMat, &
          RcVecM,RbVecM,RaVecM,RCurrentGMagnitude,RScattFacToVolts,RElectronWaveVectorMagnitude
    USE SPARA, ONLY : SBasisAtomLabel, SBasisAtomName
    USE IPARA, ONLY : IBasisAtomicNumber, IBasisAnisoDW, IMaxPossibleNAtomsUnitCell, &
         INAtomsUnitCell,ICurrentZ
    USE SPARA, ONLY : SPrintString
    
    IMPLICIT NONE
    
    INTEGER(IKIND) :: IErr,ind,jnd,knd,Imax
    INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: IAllAtomicNumber, IUniqAtomicNumber
    REAL(RKIND), ALLOCATABLE :: RAllAtomPosition(:,:), RAllOccupancy(:), RAllIsoDW(:)!, RAllAnisoDW
    REAL(RKIND), ALLOCATABLE :: RUniqAtomPosition(:,:), RUniqOccupancy(:), RUniqIsoDW(:)!, RUniqAnisoDW(:)
    REAL(RKIND) :: RScatteringFactor
    LOGICAL :: Lunique
    CHARACTER(2), DIMENSION(:), ALLOCATABLE :: SAllAtomName, SUniqAtomName
    CHARACTER(5), DIMENSION(:), ALLOCATABLE :: SAllAtomLabel, SUniqAtomLabel

    Imax=SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)
    ! All atom positions generated by symmetry, including duplicates: local variables
    ALLOCATE( RAllAtomPosition(Imax,ITHREE),SAllAtomName(Imax),SAllAtomLabel(Imax),&
        RAllOccupancy(Imax),RAllIsoDW(Imax),IAllAtomicNumber(Imax), STAT=IErr )!,RAllAnisoDW(Imax)
    ! Unique atom positions reduced from above
    ALLOCATE( RUniqAtomPosition(Imax,ITHREE),SUniqAtomName(Imax),SUniqAtomLabel(Imax),&
        RUniqOccupancy(Imax),RUniqIsoDW(Imax),IUniqAtomicNumber(Imax), STAT=IErr )!,RUniqAnisoDW(Imax)
    IF(l_alert(IErr,"UniqueAtomPositions","allocations")) RETURN
   
    !--------------------------------------------------------------------  
    ! apply symmetry elements to generate all equivalent positions 
    !--------------------------------------------------------------------  
    knd=1
    DO ind=1, SIZE(RSymVec,1)  
      DO jnd=1, SIZE(RBasisAtomPosition,1)
        RAllAtomPosition(knd,:)= MATMUL(RSymMat(ind,:,:),RBasisAtomPosition(jnd,:)) &
              + RSymVec(ind,:)
        SAllAtomLabel(knd) = SBasisAtomLabel(jnd)
        SAllAtomName(knd) = SBasisAtomName(jnd)
        RAllOccupancy(knd) = RBasisOccupancy(jnd)
        RAllIsoDW(knd) = RBasisIsoDW(jnd)
        IAllAtomicNumber(knd) = IBasisAtomicNumber(jnd)
        !RAllAnisoDW(knd) = IBasisAnisoDW(jnd)
        knd=knd+1
      END DO
    END DO
    RAllAtomPosition = MODULO(RAllAtomPosition,ONE)
    WHERE(ABS(RAllAtomPosition).LT.TINY) RAllAtomPosition = ZERO

    !--------------------------------------------------------------------  
    ! Reduce to the set of unique fractional atomic positions
    RUniqAtomPosition(1,:)= RAllAtomPosition(1,:)
    SUniqAtomLabel(1)= SAllAtomLabel(1)
    SUniqAtomName(1)= SAllAtomName(1)
    RUniqIsoDW(1) = RAllIsoDW(1)
    RUniqOccupancy(1) = RAllOccupancy(1)
    IUniqAtomicNumber(1) = IAllAtomicNumber(1)
    !RUniqAnisoDW(1) = RAllAnisoDW(1)
    jnd=1
    ! work through all possible atom coords and check for duplicates
    DO ind=1,IMax
      Lunique=.TRUE.
      DO knd=1,jnd ! check against the unique ones found so far
        IF (SUM(ABS(RAllAtomPosition(ind,:)-RUniqAtomPosition(knd,:))).LE.TINY) THEN ! position same
          IF (SAllAtomLabel(ind).EQ.SUniqAtomLabel(knd)) THEN ! Label is the same too, so not unique
            Lunique=.FALSE.
            EXIT
          END IF
        END IF
      END DO
      IF (Lunique .EQV. .TRUE.) THEN
        jnd=jnd+1
        RUniqAtomPosition(jnd,:)= RAllAtomPosition(ind,:)
        SUniqAtomLabel(jnd)= SAllAtomLabel(ind)
        SUniqAtomName(jnd)= SAllAtomName(ind)
        RUniqIsoDW(jnd) = RAllIsoDW(ind)
        RUniqOccupancy(jnd) = RAllOccupancy(ind)
        IUniqAtomicNumber(jnd) = IAllAtomicNumber(ind)!
        !RUniqAnisoDW(jnd) = RAllAnisoDW(ind)
      END IF
    END DO
    INAtomsUnitCell = jnd ! this is how many unique atoms there are in the unit cell
    WRITE(SPrintString,"(A23,I4,A6)") "The unit cell contains ",INAtomsUnitCell," atoms"
    CALL message(LS,SPrintString)

    ! make the unit cell in arrays of the correct size
    ALLOCATE(RAtomXYZ(jnd,ITHREE),STAT=IErr)
    ! atoms,in microscope reference frame, in Angstrom units
    ALLOCATE(RAtomCoordinate(jnd,ITHREE),STAT=IErr)
    ALLOCATE(SAtomLabel(jnd),STAT=IErr) ! atom label
    ALLOCATE(SAtomName(jnd),STAT=IErr) ! atom name
    ALLOCATE(RIsoDW(jnd),STAT=IErr) ! Isotropic Debye-Waller factor
    ALLOCATE(ROccupancy(jnd),STAT=IErr)
    ALLOCATE(IAtomicNumber(jnd),STAT=IErr)
    ! Anisotropic Debye-Waller factor, not yet functioning
    !ALLOCATE(IAnisoDW(jnd),STAT=IErr)
    IF(l_alert(IErr,"crystallography","allocations")) RETURN
    RAtomXYZ(:,:) = RUniqAtomPosition(1:jnd,:)
    SAtomLabel = SUniqAtomLabel(1:jnd)
    SAtomName = SUniqAtomName(1:jnd)
    RIsoDW = RUniqIsoDW(1:jnd)
    ROccupancy = RUniqOccupancy(1:jnd)
    IAtomicNumber = IUniqAtomicNumber(1:jnd)
    !IAnisoDW = RUniqAnisoDW(1:jnd)
        
    ! Finished with these variables now
    DEALLOCATE(RAllAtomPosition, SAllAtomName, RAllOccupancy, RAllIsoDW, &
         IAllAtomicNumber, RUniqAtomPosition, SUniqAtomName, &
         RUniqOccupancy, RUniqIsoDW, IUniqAtomicNumber, STAT=IErr)!, RAllAnisoDW, RUniqAnisoDW
    IF(l_alert(IErr,"UniqueAtomPositions","deallocations")) RETURN

    DO ind=1,INAtomsUnitCell    
      CALL message( LL, dbg7, "Atom ",ind)
      WRITE(SPrintString,"(A18,F8.4,F8.4,F8.4)") ": Atom position = ", RAtomXYZ(ind,:)
      CALL message( LL, dbg7, SAtomName(ind)//SPrintString )
      CALL message( LL, dbg7, "(DWF, occupancy) = ",(/ RIsoDW(ind), ROccupancy(ind) /) )
    END DO

    !--------------------------------------------------------------------
    ! calculate mean inner potential and wave vector magnitude
    !--------------------------------------------------------------------
    ! calculate the mean inner potential as the sum of scattering factors
    ! at g=0 multiplied by h^2/(2pi*m0*e*CellVolume)
    RMeanInnerPotential=ZERO
    RCurrentGMagnitude=ZERO  ! this is a global variable, sets g=0
    DO ind=1,INAtomsUnitCell
      ICurrentZ = IAtomicNumber(ind)
      CALL AtomicScatteringFactor(RScatteringFactor,IErr)
      CALL message( LL, dbg3, "Atom ",ind)
      CALL message( LL, dbg3, "f(theta) at g=0 ",RScatteringFactor)
      RMeanInnerPotential = RMeanInnerPotential+RScatteringFactor
    END DO
    RMeanInnerPotential = RMeanInnerPotential*RScattFacToVolts
    WRITE(SPrintString,FMT='(A21,F6.2,A6)') "Mean inner potential ",RMeanInnerPotential," Volts"
    SPrintString=TRIM(ADJUSTL(SPrintString))
    CALL message(LS,SPrintString)

    ! Wave vector magnitude in crystal
    ! high-energy approximation (not HOLZ compatible)
    ! K^2=k^2+U0
    RBigK= SQRT(RElectronWaveVectorMagnitude**2)!-RMeanInnerPotential)
    CALL message ( LM, dbg3, "K (Angstroms) = ",RBigK )


  END SUBROUTINE UniqueAtomPositions

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Checks a g-vector Ih,Ik,Il
  !! against the selection rules for the global variable SSpaceGroupName
  !! IFlag comes in as zero and goes out as 1 if it is an allowed reflection
  !!
  !! Major-Authors: Richard Beanland (2021)
  !!  
  SUBROUTINE SelectionRules(Ih, Ik, Il, ISel, IErr)

    ! Systematic absences from the lattice type
    ! returns 1 if the reflection is allowed
    ! returns 0 if it is forbidden
    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE SPARA, ONLY : SSpaceGroupName
    
    IMPLICIT NONE

    INTEGER (IKIND),INTENT(IN) :: Ih, Ik, Il
    INTEGER (IKIND),INTENT(INOUT) :: ISel, IErr

    ISel = 0

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
  
END MODULE crystallography_mod
