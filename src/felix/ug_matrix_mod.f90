!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: 1.3
! Date: 13-05-2024
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build: g-vector limit 
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
!! Module-description: 
!!
MODULE ug_matrix_mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UgMatrix, Absorption, GetVgContributionij, StructureFactorInitialisation
  CONTAINS

  !>
  !! Procedure-description: Calculate complex Ug matrix with absorption
  !!
  !! Major-Authors: Richard Beanland (2019)
  !!
  SUBROUTINE UgMatrix(IErr)

  USE MyNumbers
  USE message_mod

  USE MyMPI

  ! global inputs
  USE IPARA, ONLY : INhkl,IAtomicNumber,INAtomsUnitCell,ICurrentZ,IAnisoDW,IWriteFLAG
  USE RPARA, ONLY : RCurrentGMagnitude,RIsoDW,RDebyeWallerConstant,RVolume,&
    RRelativisticCorrection,Rhkl,ROccupancy,RAtomCoordinate,RgMatrix,RAnisotropicDebyeWallerFactorTensor
    !,RElectronMass,RAngstromConversion,RScattFacToVolts,RElectronCharge,RPlanckConstant,RgMatrixMagnitude
  USE CPARA, ONLY : CUgMatNoAbs,CUgMatPrime
  USE SPARA, ONLY : SPrintString
  ! global outputs
  USE CPARA, ONLY : CUgMat

  IMPLICIT NONE
    
  INTEGER(IKIND) :: ind,jnd,knd,IErr
  REAL(RKIND) :: RScatteringFactor,RPreFactor
  COMPLEX(CKIND),DIMENSION(:,:),ALLOCATABLE :: CTempMat!to avoid problems with transpose

    
  IF (my_rank.EQ.0) THEN!There may be a bug when individual cores calculate UgMat, make it the responsibility of core 0 and broadcast it
    !conversion factor from f to Ug  
    RPreFactor=RRelativisticCorrection/(PI*RVolume)
    CUgMatNoAbs = CZERO
    ! fill lower diagonal of Ug matrix(excluding absorption) with Fourier components of the potential Vg
    DO ind=2,INhkl
      DO jnd=1,ind-1
        RCurrentGMagnitude = SQRT(DOT_PRODUCT(RgMatrix(ind,jnd,:),RgMatrix(ind,jnd,:)))!RgMatrixMagnitude(ind,jnd)
        ! Old version, CVgij contribution from each atom and pseudoatom in Volts
        !CALL GetVgContributionij(RScatteringFactor,ind,jnd,CVgij,IErr) version with pseudoatoms
        DO knd=1,INAtomsUnitCell
          ICurrentZ = IAtomicNumber(knd) ! atomic number, Z, NB passed as a global variable for absorption
          ! Get scattering factor
          CALL AtomicScatteringFactor(RScatteringFactor,IErr)        
          ! Occupancy
          RScatteringFactor = RScatteringFactor*ROccupancy(knd)
          ! Debye-Waller factor
          IF(RIsoDW(knd).GT.10.OR.RIsoDW(knd).LT.0) RIsoDW(knd) = RDebyeWallerConstant!use default in felix.inp for unrealistic values in the cif
          ! Isotropic D-W factor
          ! exp(-B sin(theta)^2/lamda^2) = exp(-Bs^2) = exp(-Bg^2/16pi^2), see e.g. Bird&King
          RScatteringFactor = RScatteringFactor*EXP(-RIsoDW(knd)*(RCurrentGMagnitude**2)/(FOUR*TWOPI**2) )
          ! Here we go directly to Ug's, missing out the Fourier components of the potential Vg
          ! (formerly calculated as CVgij).  If the Vg's are desired they can be obtained from 
          ! multiplying RScatteringFactor by RScattFacToVolts,
          ! or by multiplying Ug's by (RScattFacToVolts/RPreFactor)
          ! The structure factor equation, complex Ug(ind,jnd)=sum(f*exp(-ig.r)) in Volts
          CUgMatNoAbs(ind,jnd)=CUgMatNoAbs(ind,jnd)+RPreFactor*RScatteringFactor*&
              EXP(-CIMAGONE*DOT_PRODUCT(RgMatrix(ind,jnd,:),RAtomCoordinate(knd,:)) )
        END DO
      END DO
    END DO
    ! Only the lower half of the Ug matrix was calculated, this completes the upper half
    ALLOCATE (CTempMat(INhkl,INhkl),STAT=IErr)
    IF(l_alert(IErr,"UgMatrix","allocate CTempMat")) RETURN
    CTempMat = TRANSPOSE(CUgMatNoAbs)! To avoid the bug when conj(transpose) is used
    CUgMatNoAbs = CUgMatNoAbs + CONJG(CTempMat)
    DEALLOCATE(CTempMat)
  END IF
  ind=INhkl*INhkl
  !===================================== ! Send UgMat to all cores
  CALL MPI_BCAST(CUgMatNoAbs,ind,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
  !=====================================
  
  CALL message( LM,dbg3, "Ug matrix, without absorption (nm^-2)" )!LM, dbg3
  DO ind = 1,40
    IF(IWriteFLAG.GE.2) WRITE(SPrintString,FMT='(3(I3,1X),A2,1X,6(F7.4,1X,F7.4,2X))') NINT(Rhkl(ind,:)),": ",100*CUgMatNoAbs(ind,1:6)
    CALL message( LM,dbg3, SPrintString)
  END DO

  END SUBROUTINE UgMatrix

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Select case using IAbsorbFLAG and calculate U'g prime in parallel
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE Absorption(IErr)

    ! calculate CUgMatPrime then ----> CUgMat = CUgMatNoAbs + CUgMatPrime

    USE MyNumbers
    USE message_mod
    USE MyMPI   !?? used for parallel
   
    ! global outputs
    USE CPARA, ONLY : CUgMat, CUgMatPrime

    ! global outputs - seemingly local or minor
    USE RPARA, ONLY : RCurrentGMagnitude, RCurrentB
    USE IPARA, ONLY : ICurrentZ

    ! global inputs
    USE CPARA, ONLY : CUgMatNoAbs
    USE SPARA, ONLY : SPrintString
    USE IPARA, ONLY : IAbsorbFLAG, INAtomsUnitCell, ISymmetryRelations, IEquivalentUgKey, &
          IAtomicNumber
    USE RPARA, ONLY : RAbsorptionPercentage, RAngstromConversion, RElectronCharge, &
          RElectronMass, RElectronVelocity, RPlanckConstant, RRelativisticCorrection, &
          RVolume,RIsoDW,ROccupancy,RgMatrix,Rhkl,RAtomCoordinate,RScattFacToVolts !&
          !, RgMatrixMagnitude

    IMPLICIT NONE

    ! local variables
    INTEGER(IKIND) :: ind,jnd,knd,lnd,IErr,IUniqueUgs,ILocalMin,ILocalMax
    INTEGER(IKIND),DIMENSION(2) :: ILoc
    REAL(RKIND) :: Rintegral,RfPrime,RPreFactor
    REAL(RKIND),DIMENSION(3) :: RCurrentG
    COMPLEX(CKIND),DIMENSION(:),ALLOCATABLE :: CLocalUgPrime,CUgPrime
    REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RLocalUgReal,RLocalUgImag,RUgReal,RUgImag
    INTEGER(IKIND),DIMENSION(:),ALLOCATABLE :: Ipos,Inum
    
    !--------------------------------------------------------------------  
    !  select absorption model
    !--------------------------------------------------------------------
    CUgMatPrime = CZERO
    SELECT CASE (IAbsorbFLAG)

    CASE(0) ! No absorption
      !nothing to do, CUgMatPrime is already 0

    CASE(1) ! Proportional
      CUgMatPrime = CUgMatNoAbs*EXP(CIMAGONE*PI/2)*(RAbsorptionPercentage/100_RKIND)
    
    CASE(2) 
    !--------------------------------------------------------------------  
    !  Bird & King absorption
    !--------------------------------------------------------------------

      !conversion factor from f to Ug  
      RPreFactor=RRelativisticCorrection/(PI*RVolume)
    
      !--------------------------------------------------------------------  
      ! allocate U'g calculated by this core and setup cores for absorption
      !--------------------------------------------------------------------

      ! work through unique Ug's
      IUniqueUgs = SIZE(IEquivalentUgKey)
      ! allocations for the U'g to be calculated by this core  
      ILocalMin = (IUniqueUgs*(my_rank)/p)+1
      ILocalMax = (IUniqueUgs*(my_rank+1)/p)
      ALLOCATE(Ipos(p),Inum(p),STAT=IErr)
      IF(l_alert(IErr,"Absorption","allocate Ipos")) RETURN
      ! U'g list for this core [Re,Im]
      ALLOCATE(CLocalUgPrime(ILocalMax-ILocalMin+1),STAT=IErr)
      IF(l_alert(IErr,"Absorption","allocate CLocalUgPrime")) RETURN
      ! U'g list for this core [Re]
      ALLOCATE(RLocalUgReal(ILocalMax-ILocalMin+1),STAT=IErr)
      IF(l_alert(IErr,"Absorption","allocate RLocalUgReal")) RETURN
      ! U'g list for this core [Im]
      ALLOCATE(RLocalUgImag(ILocalMax-ILocalMin+1),STAT=IErr)
      IF(l_alert(IErr,"Absorption","allocate RLocalUgImag")) RETURN
      ALLOCATE(CUgPrime(IUniqueUgs),STAT=IErr) ! complete U'g list
      IF(l_alert(IErr,"Absorption","allocate CUgPrime")) RETURN
      ALLOCATE(RUgReal(IUniqueUgs),STAT=IErr) ! complete U'g list [Re]
      IF(l_alert(IErr,"Absorption","allocate RUgReal")) RETURN
      ALLOCATE(RUgImag(IUniqueUgs),STAT=IErr) ! complete U'g list [Im]
      IF(l_alert(IErr,"Absorption","allocate RUgImag")) RETURN

      ! setup position and number for each core
      DO ind = 1,p ! p is the number of cores
        Ipos(ind) = IUniqueUgs*(ind-1)/p ! position in the MPI buffer
        Inum(ind) = IUniqueUgs*(ind)/p - IUniqueUgs*(ind-1)/p ! number of U'g components
      END DO

      !--------------------------------------------------------------------  
      ! fill each core's U'g list 
      !--------------------------------------------------------------------

      DO ind=ILocalMin,ILocalMax ! Different U'g s for each core

        CLocalUgPrime(ind-ILocalMin+1)=CZERO
        ! number of this Ug
        jnd=IEquivalentUgKey(ind)
        ! find the position of this Ug in the matrix
        ILoc = MINLOC(ABS(ISymmetryRelations-jnd))
        RCurrentG = RgMatrix(ILoc(1),ILoc(2),:) ! g-vector, local variable
        RCurrentGMagnitude = SQRT(DOT_PRODUCT(RgMatrix(ILoc(1),ILoc(2),:),RgMatrix(ILoc(1),ILoc(2),:)))! g-vector magnitude

        ! Structure factor calculation for absorptive form factors
        DO knd=1,INAtomsUnitCell
          ICurrentZ = IAtomicNumber(knd) ! Atomic number, global variable
          RCurrentB = RIsoDW(knd) ! Debye-Waller constant, global variable
          ! Get absorptive form factor f'
          CALL AbsorptiveScatteringFactor(RfPrime,IErr) ! NB uses Kirkland scattering factors
          IF(l_alert(IErr,"Absorption","CALL AbsorptiveScatteringFactor")) RETURN
          ! Occupancy
          RfPrime=RfPrime*ROccupancy(knd)
          ! Debye Waller factor, isotropic only 
          RfPrime=RfPrime*EXP(-RIsoDW(knd)*(RCurrentGMagnitude**2)/(4*TWOPI**2) )
          ! Here we go directly to Ug's, missing out the Fourier components of the potential Vg
          ! (formerly calculated as CVgPrime).  If the Vg's are desired they can be obtained from 
          ! multiplying the scattering factor RfPrime by RScattFacToVolts,
          ! or by multiplying U'g's by (RScattFacToVolts/RPreFactor)
          ! Absorptive Structure factor equation giving imaginary potential
          CLocalUgPrime(ind-ILocalMin+1)=CLocalUgPrime(ind-ILocalMin+1)+CIMAGONE*RPreFactor*Rfprime * &
                EXP(-CIMAGONE*DOT_PRODUCT(RCurrentG,RAtomCoordinate(knd,:)) )
        END DO
      END DO

      !--------------------------------------------------------------------  
      ! join the U'g lists of each core to CUgPrime
      !--------------------------------------------------------------------

      !?? RB I give up trying to MPI a complex number, do it with two REAL ones
      RLocalUgReal=REAL(CLocalUgPrime)
      RLocalUgImag=AIMAG(CLocalUgPrime)
      ! MPI gatherv the new U'g s into CUgPrime
      ! NB MPI_GATHERV(BufferToSend,No.of elements,datatype,ReceivingArray,No.of elements,)
      CALL MPI_GATHERV(RLocalUgReal,SIZE(RLocalUgReal),MPI_DOUBLE_PRECISION,&
                     RUgReal,Inum,Ipos,MPI_DOUBLE_PRECISION,&
                     root,MPI_COMM_WORLD,IErr)
      IF(l_alert(IErr,"Absorption","MPI_GATHERV RLocalUgReal")) RETURN
      CALL MPI_GATHERV(RLocalUgImag,SIZE(RLocalUgImag),MPI_DOUBLE_PRECISION,&
                     RUgImag,Inum,Ipos,MPI_DOUBLE_PRECISION,&
                     root,MPI_COMM_WORLD,IErr)
      IF(l_alert(IErr,"Absorption","MPI_GATHERV RLocalUgImag")) RETURN
      !===================================== send out the full list to all cores
      CALL MPI_BCAST(RUgReal,IUniqueUgs,MPI_DOUBLE_PRECISION,&
                     root,MPI_COMM_WORLD,IErr)
      CALL MPI_BCAST(RUgImag,IUniqueUgs,MPI_DOUBLE_PRECISION,&
                     root,MPI_COMM_WORLD,IErr)
      !=====================================

      DO ind=1,IUniqueUgs
        CUgPrime(ind)=CMPLX(RUgReal(ind),RUgImag(ind))
      END DO

      !--------------------------------------------------------------------  
      ! construct CUgMatPrime
      !--------------------------------------------------------------------

      DO ind=1,IUniqueUgs
        ! number of this Ug
        jnd=IEquivalentUgKey(ind)
        ! Fill CUgMatPrime
        WHERE(ISymmetryRelations.EQ.jnd)
          CUgMatPrime = CUgPrime(ind)
        END WHERE
        ! NB for imaginary potential U'(g)=-U'(-g)*
        WHERE(ISymmetryRelations.EQ.-jnd)
          CUgMatPrime = -CONJG(CUgPrime(ind))
        END WHERE
      END DO

    CASE DEFAULT ! Default case is no absorption, do nothing
      CALL message( LS,dbg3, "No absorption correction" )

    END SELECT

    !--------------------------------------------------------------------  
    ! the final Ug matrix with absorption
    !--------------------------------------------------------------------

    CUgMat = CUgMatNoAbs + CUgMatPrime 
    IF(my_rank.EQ.0) THEN
      CALL message( LM, dbg3, "Ug matrix, including absorption (nm^-2)" )
      DO ind = 1,40
        WRITE(SPrintString,FMT='(3(I5,1X),A2,1X,8(F9.4,1X))') NINT(Rhkl(ind,:)),": ",100*CUgMat(ind,1:4)
        CALL message( LM, dbg3, SPrintString )
      END DO
    END IF

  END SUBROUTINE Absorption

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Calculate a Fourier component of the potential Vg
  !! for atoms and pseudoatoms, NOW REDUNDANT
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE GetVgContributionij(RScatteringFactor,ind,jnd,CVgij,IErr) 

    ! ----> CVgij, used to make/update CUgMatNoAbs
    ! used each SimulateAndFit update scattering matrix Ug
    ! used once felixrefine setup via StructureFactorInitialisation 

    USE MyNumbers
    USE message_mod

    ! global outputs
    USE IPARA, ONLY : ICurrentZ

    ! global inputs
    USE IPARA, ONLY : INAtomsUnitCell,IAtomicNumber,IAnisoDW
    USE RPARA, ONLY : RIsoDW,RCurrentGMagnitude,RgMatrix,RVolume, &
          RAnisotropicDebyeWallerFactorTensor,RAtomCoordinate,ROccupancy, &
          RDebyeWallerConstant,RScattFacToVolts
    USE CPARA, ONLY : CIMAGONE
    USE RConst, ONLY : RPlanckConstant,RAngstromConversion,RElectronMass,RElectronCharge

    IMPLICIT NONE

    REAL(RKIND),INTENT(INOUT) :: RScatteringFactor
    INTEGER(IKIND),INTENT(IN) :: ind, jnd
    COMPLEX(CKIND),INTENT(OUT) :: CVgij
    INTEGER(IKIND),INTENT(OUT) :: IErr
    COMPLEX(CKIND) :: CFpseudo
    INTEGER(IKIND) :: knd, INumPseudAtoms=0

    
    CVgij = CZERO!this is in Volts
    ! Sums CVgij contribution from each atom and pseudoatom
    DO knd=1,INAtomsUnitCell
      ICurrentZ = IAtomicNumber(knd) ! atomic number, Z, NB passed as a global variable for absorption
      IF (ICurrentZ.LT.105) THEN ! It's not a pseudoatom
        ! Get scattering factor
        CALL AtomicScatteringFactor(RScatteringFactor,IErr)
        !IF (my_rank.EQ.0) PRINT*, knd,RCurrentGMagnitude,RScatteringFactor
        ! Occupancy
        RScatteringFactor = RScatteringFactor*ROccupancy(knd)
        ! Isotropic Debye-Waller factor
        IF(RIsoDW(knd).GT.10.OR.RIsoDW(knd).LT.0) RIsoDW(knd) = RDebyeWallerConstant!use default in felix.inp for unrealistic values in the cif
        ! Isotropic D-W factor
        ! exp(-B sin(theta)^2/lamda^2) = exp(-Bs^2) = exp(-Bg^2/16pi^2), see e.g. Bird&King
        RScatteringFactor = RScatteringFactor*EXP(-RIsoDW(knd) * &
              (RCurrentGMagnitude**2)/(FOUR*TWOPI**2) )
        ! The structure factor equation, complex Vg(ind,jnd)=sum(f*exp(-ig.r)) in Volts
        CVgij=CVgij+RScatteringFactor*RScattFacToVolts*EXP(-CIMAGONE*DOT_PRODUCT(RgMatrix(ind,jnd,:),&
              RAtomCoordinate(knd,:)) )
      ELSE ! pseudoatom
        INumPseudAtoms=INumPseudAtoms+1
        CALL PseudoAtom(CFpseudo,ind,jnd,INumPseudAtoms,IErr)
        ! Occupancy
        CFpseudo = CFpseudo*ROccupancy(knd)
        ! Error check: only isotropic Debye-Waller for pseudoatoms currently
        !?? DW factor: Need to work out how to get it from the REAL atom at same site
        ! assume it is the next atom in the list, for now
        CFpseudo = CFpseudo*EXP(-RIsoDW(knd+1)*(RCurrentGMagnitude**2)/(FOUR*TWOPI**2) )
        
        CVgij = CVgij + CFpseudo * &
              EXP(-CIMAGONE*DOT_PRODUCT(RgMatrix(ind,jnd,:), RAtomCoordinate(knd,:)) )

      END IF
    ENDDO

  END SUBROUTINE GetVgContributionij

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: 
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE StructureFactorInitialisation(IErr)

    ! sets up pseudoatoms,CUgMatNoAbs, mean inner potential,IEquivalentUgKey  
    USE MyNumbers
    USE message_mod
    USE MyFFTW
    USE utilities_mod, ONLY : Gaussian, Lorentzian, ReSortUgs

    ! global outputs
    USE CPARA, ONLY : CUgMatNoAbs,CUgMatPrime,CUniqueUg,CPseudoAtom,CPseudoScatt
    USE IPARA, ONLY : ICurrentZ, ISymmetryRelations
    USE RPARA, ONLY : RMeanInnerPotential,RgSumMat
    USE SPARA, ONLY : SPrintString
    USE BlochPara, ONLY : RBigK

    ! global inputs
    USE IPARA, ONLY : IInitialSimulationFLAG,INAtomsUnitCell,&
          INhkl,IAtomicNumber,IEquivalentUgKey,IWriteFLAG,IAnisoDW
    USE RPARA, ONLY : RAngstromConversion,RElectronCharge,RElectronMass,&
          RVolume,RIsoDW,ROccupancy,&
          RElectronWaveVectorMagnitude,RgMatrix,RDebyeWallerConstant,RTolerance,&
          RAtomCoordinate,Rhkl,RAnisotropicDebyeWallerFactorTensor,RScattFacToVolts,&
          RLengthX,RLengthY,RLengthZ
    USE IChannels, ONLY : IChOutWIImage

    ! global should be local to Ug.f90
    USE IPARA, ONLY : IPsize
    USE RPARA, ONLY : RCurrentGMagnitude,RPScale
          
    IMPLICIT NONE

    INTEGER(IKIND) :: ind,jnd,knd,lnd,mnd,oddindlorentz,evenindlorentz,&
          oddindgauss,evenindgauss,currentatom,IErr,Iuid,Iplan_forward,IPseudo
    INTEGER(IKIND),DIMENSION(2) :: IPos,ILoc
    COMPLEX(CKIND) :: CVgij,CFpseudo
    REAL(RKIND) :: RMeanInnerPotentialVolts,RScatteringFactor,&
          RPMag,Rx,Ry,Rr,RPalpha,RTheta,Rfold
    REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: RTempMat!to avoid problems with transpose ffs
    
    !--------------------------------------------------------------------
    ! count pseudoatoms & allocate pseudoatom arrays
    !--------------------------------------------------------------------

    ! count Pseudoatoms
    IPseudo=0
    DO jnd=1,INAtomsUnitCell
      IF (IAtomicNumber(jnd).GE.105) IPseudo = IPseudo + 1
    END DO

      !--------------------------------------------------------------------
      ! calculate pseudo potential and pseudo factor for any pseudoatoms
      !--------------------------------------------------------------------
    IF (IPseudo.GT.0) THEN!Calculate pseudoatom potentials
      ! size of the array used to calculate the pseudoatom FFT, global variable  
      IPsize=1024 ! must be an EVEN number (preferably 2^n)!
      ! matrices with Pseudoatom potentials (REAL space)  
      ALLOCATE(CPseudoAtom(IPsize,IPsize,IPseudo),STAT=IErr)
      IF(l_alert(IErr,"StructureFactorInitialisation","allocate CPseudoAtom")) RETURN
      ! matrices with Pseudoatom scattering factor (reciprocal space)
      ALLOCATE(CPseudoScatt(IPsize,IPsize,IPseudo),STAT=IErr)
      IF(l_alert(IErr,"StructureFactorInitialisation","allocate CPseudoScatt")) RETURN
      RPScale = 0.01 ! one picometre per pixel, working in Angstroms
      ! Magnitude of pseudoatom potential, in volts
      RPMag = 0.01815 ! set such that a Ja gives the same Ug matrix as a hydrogen atom
      mnd = 0 ! pseudoatom counter
        DO lnd=1,INAtomsUnitCell 
          IF (IAtomicNumber(lnd).GE.105) THEN ! we have a pseudoatom
          ! intialise pseudoatom variables
          mnd=mnd+1
          ! The Debye-Waller factor is used to determine alpha for pseudoatoms
          RPalpha=10.0*RIsoDW(lnd)
          IF (IAtomicNumber(lnd).EQ.105) Rfold=ZERO   !Ja
          IF (IAtomicNumber(lnd).EQ.106) Rfold=ONE    !Jb
          IF (IAtomicNumber(lnd).EQ.107) Rfold=TWO    !Jc
          IF (IAtomicNumber(lnd).EQ.108) Rfold=THREE  !Jd
          IF (IAtomicNumber(lnd).EQ.109) Rfold=FOUR   !Je
          IF (IAtomicNumber(lnd).EQ.110) Rfold=SIX    !Jf
          ! potential in REAL space
          DO ind=1,IPsize
            DO jnd=1,IPsize ! x&y run from e.g. -511.5 to +511.5 picometres
              Rx=RPScale*(REAL(ind-(IPsize/2))-HALF)
              Ry=RPScale*(REAL(jnd-(IPsize/2))-HALF)
              Rr=SQRT(Rx*Rx+Ry*Ry)
              Rtheta=ACOS(Rx/Rr)
              !?? Easier to make a complex input to fftw rather than fanny around
              !?? with the different format needed for a REAL input. Lazy.
              CPseudoAtom(ind,jnd,mnd)=CMPLX(RPMag*RPalpha*Rr*EXP(-RPalpha*Rr) * &
                  COS(Rfold*Rtheta),ZERO)        
            END DO
          END DO
          ! output each pseudo potential .img to check
          IF (my_rank.EQ.0) THEN 
            WRITE(SPrintString,FMT='(A15,I1,A3)') "PseudoPotential",mnd,".img"
            OPEN(UNIT=IChOutWIImage, STATUS= 'UNKNOWN', FILE=SPrintString,&
                  FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=8192)
            IF(l_alert(IErr,"StructureFactorInitialisation",&
                  "WRITE a pseudo factor.img")) RETURN
            DO jnd = 1,IPsize
              WRITE(IChOutWIImage,rec=jnd) REAL(CPseudoAtom(jnd,:,mnd))
            END DO
            CLOSE(IChOutWIImage,IOSTAT=IErr)
          END IF
          ! CPseudoScatt = a 2d fft of RPseudoAtom
          CALL dfftw_plan_dft_2d_ ( Iplan_forward, IPsize,IPsize,& 
                CPseudoAtom(:,:,mnd),CPseudoScatt(:,:,mnd),&
                FFTW_FORWARD,FFTW_ESTIMATE )
          CALL dfftw_execute_ (Iplan_forward)
          CALL dfftw_destroy_plan_ (Iplan_forward) ! could be moved to clean up?
          ! output each pseudo factor .img to check
          IF (my_rank.EQ.0) THEN
            WRITE(SPrintString,FMT='(A12,I1,A3)') "PseudoFactor",mnd,".img"
            OPEN(UNIT=IChOutWIImage, STATUS= 'UNKNOWN', FILE=SPrintString,&
                FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=8192)
            IF(l_alert(IErr,"StructureFactorInitialisation",&
                "WRITE a pseudo factor.img")) RETURN
            DO jnd = 1,IPsize
              WRITE(IChOutWIImage,rec=jnd) ABS(CPseudoScatt(jnd,:,mnd))
            END DO
            CLOSE(IChOutWIImage,IOSTAT=IErr)
          END IF
        END IF
      END DO
    END IF

    !--------------------------------------------------------------------
    ! calculate Ug matrix (excluding absorption)
    CALL UgMatrix(IErr)
    IF(l_alert(IErr,"Structure factor initialise","UgMatrix")) RETURN
   
    !--------------------------------------------------------------------
    ! calculate mean inner potential and wave vector magnitude
    !--------------------------------------------------------------------
    ! calculate the mean inner potential as the sum of scattering factors
    ! at g=0 multiplied by h^2/(2pi*m0*e*CellVolume)
    RMeanInnerPotential=ZERO
    RCurrentGMagnitude=ZERO
    DO ind=1,INAtomsUnitCell
      ICurrentZ = IAtomicNumber(ind)
      IF(ICurrentZ.LT.105) THEN ! It's not a pseudoatom
        CALL AtomicScatteringFactor(RScatteringFactor,IErr)
        CALL message( LL, dbg3, "Atom ",ind)
        CALL message( LL, dbg3, "f(theta) at g=0 ",RScatteringFactor)
        RMeanInnerPotential = RMeanInnerPotential+RScatteringFactor
      END IF
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

    !--------------------------------------------------------------------
    
    IF (IInitialSimulationFLAG.EQ.1) THEN

      !--------------------------------------------------------------------
      ! count equivalent Ugs
      !--------------------------------------------------------------------
      ! Matrix of sums of indices - for symmetry equivalence in the Ug matrix
      ALLOCATE(RgSumMat(INhkl,INhkl),STAT=IErr) 
      IF(l_alert(IErr,"felixrefine","allocate RgSumMat")) RETURN      
      ! IEquivalentUgKey is used later in absorption case 2 Bird & king
      RgSumMat = ZERO
      ! equivalent Ug's are identified by abs(h)+abs(k)+abs(l)+a*h^2+b*k^2+c*l^2...
      DO ind = 1,INhkl
        DO jnd = 1,ind
          RgSumMat(ind,jnd)=ABS(Rhkl(ind,1)-Rhkl(jnd,1))+ABS(Rhkl(ind,2)-Rhkl(jnd,2))+ABS(Rhkl(ind,3)-Rhkl(jnd,3))+ &
            RLengthX*(Rhkl(ind,1)-Rhkl(jnd,1))**TWO+RLengthY*(Rhkl(ind,2)-Rhkl(jnd,2))**TWO+ &
            RLengthZ*(Rhkl(ind,3)-Rhkl(jnd,3))**TWO
        END DO
      END DO
      ! it's symmetric
      ALLOCATE (RTempMat(INhkl,INhkl),STAT=IErr)
      RTempMat = TRANSPOSE(RgSumMat)
      RgSumMat = RgSumMat+RTempMat
      DEALLOCATE (RTempMat)
      CALL message ( LL, dbg3, "hkl: g Sum matrix" )
      DO ind =1,16
        IF(IWriteFLAG.GE.4) WRITE(SPrintString,FMT='(3(I2,1X),A2,1X,12(F6.1,1X))') NINT(Rhkl(ind,:)),": ",RgSumMat(ind,1:12)
        CALL message ( LL, dbg3, SPrintString )!LM, dbg3
      END DO

      ISymmetryRelations = 0_IKIND 
      Iuid = 0_IKIND 
      DO jnd = 1,INhkl
        DO ind = 1,INhkl
          IF(ISymmetryRelations(ind,jnd).NE.0) THEN
            CYCLE
          ELSE
            Iuid = Iuid + 1_IKIND
            ! fill the symmetry relation matrix with incrementing numbers
            ! that have the sign of the imaginary part
            WHERE (ABS(RgSumMat-RgSumMat(ind,jnd)).LE.RTolerance)
              ISymmetryRelations = Iuid*SIGN(1_IKIND,NINT(AIMAG(CUgMatNoAbs)/(TINY)))
            END WHERE
          END IF
        END DO
      END DO
      DEALLOCATE (RgSumMat)
      WRITE(SPrintString,FMT='(I6,A25)') Iuid," unique structure factors"
      SPrintString=TRIM(ADJUSTL(SPrintString))
      CALL message ( LS, SPrintString )
      CALL message ( LM, dbg3, "hkl: symmetry matrix" )
      DO ind =1,40
        WRITE(SPrintString,FMT='(3(I4,1X),A2,1X,16(I4,1X))') NINT(Rhkl(ind,:)),": ",ISymmetryRelations(ind,1:16)
        CALL message ( LM,dbg3, SPrintString )
      END DO

      ! link each key with its Ug, from 1 to the number of unique Ug's Iuid
      ALLOCATE(IEquivalentUgKey(Iuid),STAT=IErr)
      IF(l_alert(IErr,"StructureFactorInitialisation","allocate IEquivalentUgKey")) RETURN
      ALLOCATE(CUniqueUg(Iuid),STAT=IErr)
      IF(l_alert(IErr,"StructureFactorInitialisation","allocate CUniqueUg")) RETURN

      DO ind = 1,Iuid
        ILoc = MINLOC(ABS(ISymmetryRelations-ind))
        IEquivalentUgKey(ind) = ind
        CUniqueUg(ind) = CUgMatNoAbs(ILoc(1),ILoc(2))
      END DO

      ! put them in descending order of magnitude
      ! IEquivalentUgKey is used later in absorption case 2 Bird & king  
      CALL ReSortUgs(IEquivalentUgKey,CUniqueUg,Iuid) ! modifies those arrays
    
    END IF

    RETURN

  END SUBROUTINE StructureFactorInitialisation

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  !>
  !! Procedure-description: Choose scattering factors using IScatterFactorMethodFLAG
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  SUBROUTINE  AtomicScatteringFactor(RScatteringFactor,IErr)  

    USE MyNumbers
    USE utilities_mod, ONLY : Gaussian, Lorentzian

    ! global inputs
    USE IPARA, ONLY : ICurrentZ, IScatterFactorMethodFLAG
    USE RPARA, ONLY : RCurrentGMagnitude, RScattFactors
    
    IMPLICIT NONE

    REAL(RKIND),INTENT(OUT) :: RScatteringFactor
    INTEGER(IKIND) :: ind,jnd,knd,IErr
    
    ! select scattering factor method
    RScatteringFactor = ZERO
    SELECT CASE (IScatterFactorMethodFLAG)
            
    CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians
      ! NB Kirkland scattering factor is in Angstrom units
      ! NB atomic number and g-vector passed as global variables
      RScatteringFactor = Kirkland(RCurrentGMagnitude)
        
    CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
      DO ind = 1,4
        ! Peng Method uses summation of 4 Gaussians
        RScatteringFactor = RScatteringFactor + &
          GAUSSIAN(RScattFactors(ICurrentZ,ind),RCurrentGMagnitude,ZERO, & 
          SQRT(2/RScattFactors(ICurrentZ,ind+4)),ZERO)
      END DO
    
    CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)
      ! NB DoyleTurner scattering factor is in Angstrom units
      ! NB atomic number and g-vector passed as global variables
      RScatteringFactor = DoyleTurner(RCurrentGMagnitude)
        
    CASE(3) ! 10 Parameter method with Scattering Parameters from Lobato et al. 2014
      !?? update github wiki
      DO ind = 1,5
        jnd=ind+5
        RScatteringFactor = RScatteringFactor + &
          LORENTZIAN(RScattFactors(ICurrentZ,ind)* &
         (TWO+RScattFactors(ICurrentZ,jnd)*(RCurrentGMagnitude**TWO)),ONE, &
          RScattFactors(ICurrentZ,jnd)*(RCurrentGMagnitude**TWO),ZERO)
      END DO

    END SELECT

  END SUBROUTINE AtomicScatteringFactor

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Returns the absorptive form factor f'
  !! evaluated for RCurrentG, RCurrentB and ICurrentZ using the Bird & King method
  !! this is a integration over k-space called 
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  SUBROUTINE AbsorptiveScatteringFactor(RfPrime,IErr) 
    ! used in each (case 2 Bird & King) absorption

    USE MyNumbers
    USE message_mod
    USE RPARA, ONLY : RPlanckConstant,RAngstromConversion,RElectronMass,RElectronVelocity

    IMPLICIT NONE

    INTEGER(IKIND), PARAMETER :: inf=1
    INTEGER(IKIND), PARAMETER :: limit=500
    INTEGER(IKIND), PARAMETER :: lenw= limit*4
    INTEGER(IKIND) :: IErr,Ieval,last, iwork(limit)
    REAL(RKIND) :: RAccuracy,RError,RfPrime,RAbsPreFactor
    REAL(RKIND) :: work(lenw)
    
    RAbsPreFactor = TWO*RPlanckConstant*RAngstromConversion/(RElectronMass*RElectronVelocity)
    RAccuracy=0.00000001D0 ! accuracy of integration
    ! use single integration IntegrateBK as an external function of one variable
    ! Quadpack integration 0 to infinity
    CALL dqagi(IntegrateBK,ZERO,inf,0,RAccuracy,RfPrime,RError,Ieval,IErr,&
         limit, lenw, last, iwork, work )
    ! The integration required is actually -inf to inf in 2 dimensions
    ! We used symmetry to just do 0 to inf, so multiply by 4
    RfPrime=RfPrime*4*RAbsPreFactor
    
  END SUBROUTINE AbsorptiveScatteringFactor

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !>
  !! Procedure-description: Returns a PseudoAtom scattering factor 
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  SUBROUTINE PseudoAtom(CFpseudo,i,j,k,IErr)

    ! Reads a scattering factor from the kth Stewart pseudoatom in CPseudoScatt 
    ! RCurrentGMagnitude is passed as a global variable
    ! IPsize is the size, RPscale the scale factor of:
    ! the matrix holding the scattering factor, global variable
    ! i and j give the location of the g-vector in the appropriate matrices

    USE MyNumbers

    USE RPARA, ONLY : RElectronWaveLength,RCurrentGMagnitude,RPscale,RgMatrix
    USE IPARA, ONLY : IPsize
    USE CPARA, ONLY : CPseudoScatt
    
    IMPLICIT NONE
    
    COMPLEX(CKIND),INTENT(OUT) :: CFpseudo
    INTEGER(IKIND) :: i,j,k,Ix,Iy,IErr
    REAL(RKIND) :: RPMag,Rx,Ry,Rr,Rtheta

    IF (RElectronWaveLength*ABS(RCurrentGMagnitude)*REAL(IPsize)*RPscale/TWOPI.GE.ONE) THEN
      ! g vector is out of range of the fft
      CFpseudo=CZERO
    ELSE ! Find the pixel corresponding to g
      Ix=NINT(RElectronWaveLength*RgMatrix(i,j,1) * &
            REAL(IPsize)*REAL(IPsize)*RPscale/(TWO*TWOPI))
      Iy=NINT(RElectronWaveLength*RgMatrix(i,j,2) * &
            REAL(IPsize)*REAL(IPsize)*RPscale/(TWO*TWOPI))
      ! fft has the origin at [0,0], negative numbers wrap around from edges
      IF (Ix.LE.0) Ix=Ix+IPsize 
      IF (Iy.LE.0) Iy=Iy+IPsize
      CFpseudo=CPseudoScatt(Ix,Iy,k)
    END IF
    
  END SUBROUTINE PseudoAtom

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Returns a Doyle-Turner scattering factor 
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  FUNCTION DoyleTurner(Rg)
    ! used in each (case 2 DoyleTurner) AtomicScatteringFactor
    ! P.A. Doyle and P.S. Turner Acta Cryst A24,390 (1968)
    ! ICurrentZ is atomic number, passed as a global variable
    ! Rg is magnitude of scattering vector in 1/A
    ! (NB exp(-i*g.r), physics negative convention), global variable
    ! DoyleTurner scattering factor is in Angstrom units
    USE MyNumbers

    ! global inputs
    USE IPARA, ONLY : ICurrentZ 
    USE RPARA, ONLY : RScattFactors
    
    IMPLICIT NONE
    
    INTEGER(IKIND) :: ind,IErr
    REAL(RKIND) :: DoyleTurner,Ra,Rb,Rs,Rg

    ! NB DoyleTurner scattering factors are calculated using s = sin(theta)/lambda = (d*)/2 = g/4pi
    Rs = Rg / FOURPI
    DoyleTurner=ZERO;
    
    DO ind = 1,4
      Ra=RScattFactors(ICurrentZ,ind*2-1);
      Rb=RScattFactors(ICurrentZ,ind*2);
      DoyleTurner = DoyleTurner + Ra*EXP(-(Rb*Rs**2))
    END DO
    
  END FUNCTION DoyleTurner
  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Returns a Kirkland scattering factor 
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  FUNCTION Kirkland(Rg)
    ! used in each (case 0 Kirkland) AtomicScatteringFactor
    ! used in each (case 2 Bird & King) absorption via BirdKing

    ! From Appendix C of Kirkland, "Advanced Computing in Electron Microscopy", 2nd ed.
    ! ICurrentZ is atomic number, passed as a global variable
    ! Rg is magnitude of scattering vector in 1/A
    ! (NB exp(-i*g.r), physics negative convention), global variable
    ! Kirkland scattering factor is in Angstrom units
    USE MyNumbers

    ! global inputs
    USE IPARA, ONLY : ICurrentZ 
    use RPARA, ONLY : RScattFactors
    
    IMPLICIT NONE
    
    INTEGER(IKIND) :: ind,IErr
    REAL(RKIND) :: Kirkland,Ra,Rb,Rc,Rd,Rq,Rg

    ! NB Kirkland scattering factors are calculated in the optics convention exp(2*pi*i*q.r)
    Rq = Rg / TWOPI
    Kirkland=ZERO;
    ! Equation C.15
    DO ind = 1,3
      Ra=RScattFactors(ICurrentZ,ind*2-1);
      Rb=RScattFactors(ICurrentZ,ind*2);
      Rc=RScattFactors(ICurrentZ,ind*2+5);
      Rd=RScattFactors(ICurrentZ,ind*2+6);
      Kirkland = Kirkland + Ra/((Rq**2)+Rb) + Rc*EXP(-(Rd*Rq**2))
    END DO
    
  END FUNCTION Kirkland

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Used as part of numerical integration to calculate
  !! absorptive form factor f' for Bird & King absorption method
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  FUNCTION IntegrateBK(Sy) 
    ! used in each (case 2 Bird & King) absorption (indirectly)
    ! 'CALL dqagi(IntegrateBK,---)' each DoubleIntegrateBK

    USE MyNumbers

    USE RPARA, ONLY : RSprimeY

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr,Ieval
    INTEGER(IKIND), PARAMETER :: inf = 1
    INTEGER(IKIND), PARAMETER :: limit = 500
    INTEGER(IKIND), PARAMETER :: lenw = limit*4
    REAL(RKIND) :: RAccuracy, RError, IntegrateBK, Sy
    INTEGER(IKIND) last, iwork(limit)
    REAL(RKIND) work(lenw)

    RSprimeY = Sy
    RAccuracy = 0.00000001D0 ! accuracy of integration
    ! use BirdKing as an external function of one variable
    ! Quadpack integration 0 to infinity
    CALL dqagi(BirdKing,ZERO,inf,0,RAccuracy,IntegrateBK,RError,Ieval,IErr,&
         limit, lenw, last, iwork, work )
    
  END FUNCTION IntegrateBK

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Used as part of numerical integration to calculate
  !! absorptive form factor f' for Bird & King absorption method. Defines a Bird & King 
  !! integrand to calculate an absorptive scattering factor 
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  FUNCTION BirdKing(RSprimeX)
    ! used in each (case 2 Bird & King) absorption (indirectly)
    ! 'CALL dqagi(BirdKing,---)' each IntegrateBK

    ! From Bird and King, Acta Cryst A46, 202 (1990)
    ! ICurrentZ is atomic number, global variable
    ! RCurrentB is Debye-Waller constant b=8*pi*<u^2>, 
    ! where u is mean square thermal vibration amplitude in Angstroms, global variable
    ! RCurrentGMagnitude is magnitude of scattering vector in 1/A 
    ! NB exp(-i*g.r), physics negative convention, global variable
    ! RSprime is dummy parameter for integration [s'x s'y]
    ! NB can't print from here as it is called EXTERNAL in Integrate
    USE MyNumbers

    ! global inputs
    USE RPARA, ONLY : RCurrentB, RCurrentGMagnitude, RSprimeY
    
    IMPLICIT NONE
    
    INTEGER(IKIND) :: ind
    REAL(RKIND):: BirdKing, Rs, Rg1, Rg2, RsEff
    REAL(RKIND), INTENT(IN) :: RSprimeX
    REAL(RKIND),DIMENSION(2) :: RGprime
    
    ! NB Kirkland scattering factors in optics convention
    RGprime = 2 * TWOPI * [RSprimeX,RSprimeY]
    ! Since [s'x s'y]  is a dummy parameter for integration I can assign s'x //g
    Rg1 = SQRT( (RCurrentGMagnitude/2+RGprime(1))**2 + RGprime(2)**2 )
    Rg2 = SQRT( (RCurrentGMagnitude/2-RGprime(1))**2 + RGprime(2)**2 )
    RsEff = RSprimeX**2+RSprimeY**2-RCurrentGMagnitude**2/(16*TWOPI**2)
    BirdKing = Kirkland(Rg1)*Kirkland(Rg2)*(1-EXP(-2*RCurrentB*RsEff ) )
    
  END FUNCTION BirdKing

END MODULE ug_matrix_mod
