!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: 2.0/ 1.15 /
! Date: 19-12-2022 16-01-2019
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build: cRED Mode F: test different lattice types" 
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
!! Module-description: Solely for BlochCoefficientCalculation procedure
!!
MODULE bloch_mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BlochCoefficientCalculation

  CONTAINS

  !>
  !! Procedure-description: Simulates the electron beam and calculates Bloch
  !! coefficients for a specified pixel for each LACBED pattern and for each thickness.
  !! This can be used in parallel such that each core calculates different pixels.
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE BlochCoefficientCalculation(IYPixelIndex,IXPixelIndex,IPixelNumber,&
                    IFirstPixelToCalculate,nBeams, RThickness, RKn, IErr)

    ! accesses procedure ZGEMM
    USE MyNumbers
    USE MyMPI
    USE message_mod
  
    ! globals - output
    USE RPara, ONLY : RIndividualReflections ! RIndividualReflections( LACBED_ID, thickness_ID, local_pixel_ID )
    USE CPara, ONLY : CAmplitudeandPhase

    ! globals - input  
    USE CPara, ONLY : CUgMat
    USE RPara, ONLY : RDeltaK,RDeltaThickness,RInitialThickness,RNormDirM,RgDotNorm,RgPool,&
                      RgPoolMag,Rhkl,RgMatrix,RMeanInnerPotential,RDevPara,RBigK
    USE IPara, ONLY : IHolzFLAG,IMinStrongBeams,IMinWeakBeams,&
                      INoOfHKLsFrame,ISizeX,ISizeY,IThicknessCount,INhkl,&
                      IhklsFrame
    USE SPARA, ONLY : SPrintString
    USE RPARA, ONLY : RgDotNorm

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(IN) :: IYPixelIndex,IXPixelIndex,IPixelNumber,&
          IFirstPixelToCalculate
    INTEGER(IKIND),INTENT(OUT) :: nBeams, IErr
	REAL, INTENT(OUT) :: RThickness, RKn
    
    COMPLEX(CKIND),ALLOCATABLE :: CBeamProjectionMatrix(:,:),&
          CDummyBeamMatrix(:,:),CUgSgMatrix(:,:),CEigenVectors(:,:),CEigenValues(:),&
          CInvertedEigenVectors(:,:),CAlphaWeightingCoefficients(:),&
          CEigenValueDependentTerms(:,:)
    COMPLEX(CKIND) :: CFullWaveFunctions(INhkl)
    REAL(RKIND) :: RFullWaveIntensity(INhkl),RTiltedK(ITHREE)
    INTEGER(IKIND) :: IStrongBeamList(INhkl),IWeakBeamList(INhkl),&
          nWeakBeams
    INTEGER(IKIND) :: ind,jnd,knd,pnd,IThickness,IThicknessIndex,ILowerLimit,&
          IUpperLimit       
    REAL(RKIND) :: Rk0(3),RkPrime(3),RK,RKg,Rd1,Rd2
    COMPLEX(CKIND) sumC,sumD
    COMPLEX(CKIND), DIMENSION(:,:), ALLOCATABLE :: CBeamTranspose,CUgMatPartial,CDummyEigenVectors
    COMPLEX(CKIND), DIMENSION(:,:), ALLOCATABLE :: CStructureMatrix
    CHARACTER*40 surname
    CHARACTER*100 SindString,SjndString,SPixelCount,SnBeams,SWeakBeamIndex
    
    IErr=0

    ! TiltedK is the vector of the incoming tilted beam
    ! in units of (1/A), in the microscope ref frame(NB exp(i*k.r), physics convention)
    ! x-position in k-space
    RTiltedK(1)= (REAL(IXPixelIndex,RKIND)-0.5_RKIND*REAL(ISizeX,RKIND)-0.5_RKIND)*RDeltaK
    ! y-position in k-space
    RTiltedK(2)= (REAL(IYPixelIndex,RKIND)-0.5_RKIND*REAL(ISizeY,RKIND)-0.5_RKIND)*RDeltaK 
    RTiltedK(3)= SQRT(RBigK**2 - RTiltedK(1)**2 - RTiltedK(2)**2) 
    RKn = DOT_PRODUCT(RTiltedK,RNormDirM)
    Rk0 = ZERO
    RkPrime=ZERO
    !IF(my_rank.EQ.0) PRINT*,RTiltedK
    ! Compute the deviation parameter for reflection pool
    ! NB RDevPara is in units of (1/A)
    ! in the microscope ref frame(NB exp(i*s.r), physics convention)
    DO knd=1,INhkl
      ! Version without small angle approximation
      ! Sg=(g/k)*[2(k^2-k0.k')]^0.5
      ! k0 is defined by the Bragg condition
      Rk0(1) = -RgPoolMag(knd)/2
      Rk0(3) = SQRT(RBigK**2-Rk0(1)**2)
      ! k' is from RTiltedK
      RkPrime(1)=DOT_PRODUCT(RTiltedK,RgPool(knd,:))/RgPoolMag(knd)!Gives NaN for 000
      RkPrime(3) = SQRT(RBigK**2-RkPrime(1)**2)
      RDevPara(knd)=-SIGN(ONE,(2*DOT_PRODUCT(RgPool(knd,:),RTiltedK)+RgPoolMag(knd)**2))*&
                    RgPoolMag(knd)*SQRT(2*(RBigK**2-DOT_PRODUCT(Rk0,RkPrime)))/RBigK
      IF (RgPoolMag(knd).EQ.ZERO) RDevPara(knd)=ZERO!Avoid NaN for 000
      !IF(my_rank.EQ.0) PRINT*, knd,RgPool(knd,1),RgPool(knd,2)
      !IF(my_rank.EQ.0) PRINT*, "new",RDevPara(knd),&
      !      SIGN(ONE,(2*DOT_PRODUCT(RgPool(knd,:),RTiltedK)-RgPoolMag(knd)**2))
      ! Old version, Sg parallel to z: Sg=-[k'z+gz-sqrt( (k'z+gz)^2-2k'.g-g^2)]
      !RDevPara(knd)= -RTiltedK(3)-RgPool(knd,3)+&
      !  SQRT( (RTiltedK(3)+RgPool(knd,3))**2 - &
      !  2*DOT_PRODUCT(RgPool(knd,:),RTiltedK) - RgPoolMag(knd)**2 )
      !IF(my_rank.EQ.0) PRINT*, "old", RDevPara(knd)
      ! Debugging output
      IF(knd.EQ.2.AND.IYPixelIndex.EQ.10.AND.IXPixelIndex.EQ.10) THEN
        CALL message(LM,"RBigK ",RBigK)!LM,dbg7
        CALL message(LM,"Rhkl(knd) ",Rhkl(knd:knd,:))
        CALL message(LM,"RgPool(knd) ",RgPool(knd:knd,:))
        CALL message(LM,"RTiltedK ",RTiltedK)
        CALL message(LM,"RDevPara ",RDevPara(knd))
      END IF
    END DO

    ! select only those beams where the Ewald sphere is close to the
    ! reciprocal lattice, i.e. within RBSMaxDeviationPara
    CALL StrongAndWeakBeamsDetermination(INhkl,IMinWeakBeams,&
                    IMinStrongBeams,CUgMat,&
                    IStrongBeamList,IWeakBeamList,nBeams,nWeakBeams,IErr)
    IF(l_alert(IErr,"BlochCoefficientCalculation",&
          "StrongAndWeakBeamsDetermination()")) RETURN
    CALL message(LXL,dbg7,"strong beams",nBeams)
    CALL message(LXL,dbg7,"weak beams",nWeakBeams)
    CALL message(LXL,dbg7,"INhkl",INhkl)

    !--------------------------------------------------------------------
    ! ALLOCATE memory for eigen problem
    !--------------------------------------------------------------------

    ! now nBeams determined, allocate complex arrays
    ALLOCATE( CBeamProjectionMatrix(nBeams,INhkl), STAT=IErr )
    ALLOCATE( CDummyBeamMatrix(nBeams,INhkl), STAT=IErr )
    ALLOCATE( CUgSgMatrix(nBeams,nBeams), STAT=IErr )
    ALLOCATE( CEigenValues(nBeams), STAT=IErr )
    ALLOCATE( CEigenVectors(nBeams,nBeams), STAT=IErr )
    ALLOCATE( CDummyEigenVectors(nBeams,nBeams), STAT=IErr )
    ALLOCATE( CInvertedEigenVectors(nBeams,nBeams), STAT=IErr )
    ALLOCATE( CBeamTranspose(INhkl,nBeams), STAT=IErr )
    ALLOCATE( CUgMatPartial(INhkl,nBeams), STAT=IErr )
    ALLOCATE( CAlphaWeightingCoefficients(nBeams), STAT=IErr )
    ALLOCATE( CEigenValueDependentTerms(nBeams,nBeams), STAT=IErr )
    ALLOCATE( CStructureMatrix(nBeams, nBeams), STAT = IErr)
    IF(l_alert(IErr,"BlochCoefficientCalculation","allocations")) RETURN

    ! compute the effective Ug matrix by selecting only those beams
    ! for which IStrongBeamList has an entry
    CBeamProjectionMatrix= CZERO
    DO knd=1,nBeams
      CBeamProjectionMatrix(knd,IStrongBeamList(knd))=CONE
    ENDDO

    CUgSgMatrix = CZERO
    CBeamTranspose=TRANSPOSE(CBeamProjectionMatrix)
    ! reduce the matrix to just include strong beams using some nifty matrix multiplication
    ! CUgMatPartial = CUgMat * CBeamTranspose
    CALL ZGEMM('N','N',INhkl,nBeams,INhkl,CONE,CUgMat, &
              INhkl,CBeamTranspose,INhkl,CZERO,CUgMatPartial,INhkl)
    ! CUgSgMatrix = CBeamProjectionMatrix * CUgMatPartial
    CALL ZGEMM('N','N',nBeams,nBeams,INhkl,CONE,CBeamProjectionMatrix, &
              nBeams,CUgMatPartial,INhkl,CZERO,CUgSgMatrix,nBeams)

    !--------------------------------------------------------------------
    ! Constructing the UgSg (aka Structure) matrix
    !--------------------------------------------------------------------
    ! replace the diagonal parts with strong beam deviation parameters
    !The 4pi^2 is a result of using h, not hbar, in the conversion from VG(ij) to Ug(ij)
    DO ind=1,nBeams
      CUgSgMatrix(ind,ind) = TWO*RBigK*RDevPara(IStrongBeamList(ind))/(TWOPI*TWOPI)
    ENDDO

    ! Weak beams: add perturbatively for the 1st column (sumC) and
    ! the diagonal elements (sumD)
    DO knd=2,nBeams
      sumC=CZERO
      sumD=CZERO
      DO ind=1,nWeakBeams
        ! Zuo&Weickenmeier Ultramicroscopy 57 (1995) 375-383 eq.4
        sumC=sumC + &
        CUgMat(IStrongBeamList(knd),IWeakBeamList(ind))*&
        CUgMat(IWeakBeamList(ind),1)/(TWO*RBigK*RDevPara(IWeakBeamList(ind)))
        ! Zuo&Weickenmeier Ultramicroscopy 57 (1995) 375-383 eq.5
        sumD = sumD + &
        CUgMat(IStrongBeamList(knd),IWeakBeamList(ind))*&
        CUgMat(IWeakBeamList(ind),IStrongBeamList(knd))/&
        (TWO*RBigK*RDevPara(IWeakBeamList(ind)))
      ENDDO
      ! Replace the Ug's
      WHERE (CUgSgMatrix.EQ.CUgSgMatrix(knd,1))
        CUgSgMatrix = CUgSgMatrix(knd,1) - sumC
      END WHERE
      ! Replace the Sg's
      CUgSgMatrix(knd,knd)= CUgSgMatrix(knd,knd) - TWO*RBigK*sumD/(TWOPI*TWOPI)
    ENDDO
    !Divide by 2K so off-diagonal elementa are Ug/2K, diagonal elements are Sg, Spence's (1993) 'Structure matrix'
    CUgSgMatrix = TWOPI*TWOPI*CUgSgMatrix/(TWO*RBigK)

    ! Recalculation of UgSg matrix for tilted foil
    DO ind = 1, nBeams
      Rd1 = SQRT(1+RgDotNorm(IStrongBeamList(ind))/RKn)
      DO jnd = 1, nBeams
        Rd2 = SQRT(1+RgDotNorm(IStrongBeamList(jnd))/RKn)
IF(CStructureMatrix(ind,jnd).NE.CStructureMatrix(ind,jnd))PRINT*,my_rank,"NaN!!! i,j= ",ind,jnd
        CStructureMatrix(ind,jnd) = CUgSgMatrix(ind,jnd)/(Rd1*Rd2)
      END DO
    END DO

    !--------------------------------------------------------------------
    ! diagonalize the UgMatEffective
    !--------------------------------------------------------------------
    CALL EigenSpectrum(nBeams,CStructureMatrix,CEigenValues(:),CEigenVectors(:,:),IErr)
    IF(l_alert(IErr,"BlochCoefficientCalculation","EigenSpectrum()")) RETURN
    ! NB destroys CUgSgMatrix

    !THIS MAY BE REDUNDANT?
    IF (IHolzFLAG.EQ.1) THEN ! higher order laue zone included so adjust Eigen values/vectors
      CEigenValues = CEigenValues * RKn/RBigK
      DO knd = 1,nBeams
        CEigenVectors(knd,:) = CEigenVectors(knd,:) / &
              SQRT(1+RgDotNorm(IStrongBeamList(knd))/RKn)
      END DO
    END IF

    ! Invert the EigenVector matrix
    CDummyEigenVectors = CEigenVectors
    CALL INVERT(nBeams,CDummyEigenVectors(:,:),CInvertedEigenVectors,IErr)

    !--------------------------------------------------------------------
    ! fill RIndividualReflections( LACBED_ID , thickness_ID, local_pixel_ID ) 
    !--------------------------------------------------------------------
   
    ! Calculate intensities for different specimen thicknesses
    DO IThicknessIndex=1,IThicknessCount,1

      RThickness = RInitialThickness + REAL((IThicknessIndex-1),RKIND)*RDeltaThickness
      IThickness = NINT(RThickness,IKIND)

      CALL CreateWaveFunctions(RThickness,RKn,RFullWaveIntensity,CFullWaveFunctions,&
                    INhkl,nBeams,IStrongBeamList,CEigenVectors,CInvertedEigenVectors,CEigenValues,IErr)
      IF(l_alert(IErr,"BlochCoefficientCalculation","CreateWaveFunctions")) RETURN
      ! Collect Intensities from all thickness for later writing
      DO ind = 1,INoOfHKLsFrame
        RIndividualReflections(ind,IThicknessIndex,&
          (IPixelNumber-IFirstPixelToCalculate)+1) = RFullWaveIntensity(IhklsFrame(ind))
        END DO
    END DO

    ! DEALLOCATE eigen problem memory
    DEALLOCATE(CUgSgMatrix,CBeamTranspose, CUgMatPartial, &
         CInvertedEigenVectors, CAlphaWeightingCoefficients, &
         CEigenValues,CEigenVectors,CEigenValueDependentTerms, &
         CBeamProjectionMatrix, CDummyBeamMatrix, CStructureMatrix, STAT=IErr)
    IF(l_alert(IErr,"BlochCoefficientCalculation","deallocating arrays")) RETURN
    
  END SUBROUTINE BlochCoefficientCalculation

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Calculates intensity for a specific incident beam 
  !! orientation and a single thickness
  !! Now accounts for non-parallel incident electron beam
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016), Chris Cumming (2022)
  !!
  SUBROUTINE CreateWaveFunctions(RThickness,RKn,RFullWaveIntensity,CFullWaveFunctions,&
                    INhkl,nBeams,IStrongBeamList,CEigenVectors,CInvertedEigenVectors,CEigenValues,IErr)

    USE MyNumbers
    USE MyMPI
    USE message_mod

    USE RPARA, ONLY : RgDotNorm

    IMPLICIT NONE
    
    REAL(RKIND),INTENT(IN) :: RThickness, RKn
    REAL(RKIND),INTENT(OUT) :: RFullWaveIntensity(INhkl)
    COMPLEX(CKIND),INTENT(OUT) :: CFullWaveFunctions(INhkl) 
    INTEGER(IKIND),INTENT(IN) :: INhkl,nBeams,IStrongBeamList(INhkl)
    COMPLEX(CKIND),INTENT(IN) :: CEigenVectors(nBeams,nBeams),CInvertedEigenVectors(nBeams,nBeams),CEigenValues(nBeams)
    INTEGER(IKIND),INTENT(OUT) :: IErr 
    COMPLEX(CKIND) :: CPsi0(nBeams),CWaveFunctions(nBeams),CEigenValueDependentTerms(nBeams,nBeams), &
          CMmatrix(nBeams,nBeams), CInvertedM(nBeams,nBeams)
    INTEGER(IKIND) :: ind
    
    IErr=0
    ! The top surface boundary conditions, only valid for singular incident electron beam
    CPsi0 = CZERO ! All diffracted beams are zero
    CPsi0(1) = CONE ! The 000 beam has unit amplitude

    ! Form eigenvalue diagonal matrix
    CEigenValueDependentTerms= CZERO
    DO ind=1,nBeams
      CEigenValueDependentTerms(ind,ind) = &
           EXP((CIMAGONE*CMPLX(RThickness,ZERO,CKIND)*CEigenValues(ind)))
    END DO

    ! The M matrix and its inverse account for surface normal, see Zuo & Spence, or Palatinus 2015(?)
    CMmatrix = CZERO
    CInvertedM = CZERO
    DO ind = 1, nBeams
      CMmatrix(ind,ind) = 1/SQRT(1+RgDotNorm(IStrongBeamList(ind))/RKn)
      CInvertedM(ind,ind) = SQRT(1+RgDotNorm(IStrongBeamList(ind))/RKn)
    END DO

    !--------------------------------------------------------------------
    ! Wave functions
    !--------------------------------------------------------------------
    ! we get all nBeams at once from this matrix equation
    CWaveFunctions(:) = MATMUL(CMmatrix, MATMUL(CEigenVectors, MATMUL(CEigenValueDependentTerms, &
    MATMUL(CInvertedEigenVectors, MATMUL(CInvertedM, CPsi0)))))

    !?? possible time saving here by only calculating the (tens of) output
    !?? reflections rather than all strong beams (hundreds)
    !--------------------------------------------------------------------
    ! we return the wave functions/intensities in CFullWaveFunctions and RFullWaveIntensity
    ! which have size INhkl, the beam pool
    CFullWaveFunctions = CZERO
    RFullWaveIntensity = ZERO
    DO ind=1,nBeams
       CFullWaveFunctions(IStrongBeamList(ind)) = CWaveFunctions(ind)
       RFullWaveIntensity(IStrongBeamList(ind)) = CWaveFunctions(ind)*CONJG(CWaveFunctions(ind))
    ENDDO

  END SUBROUTINE CreateWavefunctions

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !>
  !! Procedure-description: Determines number of weak and strong beams. Uses Sg and
  !! perturbation strengths and iterates over the number of weak and strong until
  !! there are enough.
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE StrongAndWeakBeamsDetermination(INhkl,IMinWeakBeams,&
                    IMinStrongBeams,CUgMat,&
                    IStrongBeamList,IWeakBeamList,nBeams,nWeakBeams,IErr)
    
    ! select only those beams where the Ewald sphere is close to the
    ! reciprocal lattice, i.e. within RBSMaxDeviationPara

    USE MyNumbers
    USE MyMPI
    USE message_mod

    USE RPara, ONLY : RDevPara

    INTEGER(IKIND),INTENT(IN) :: INhkl
    COMPLEX(CKIND),DIMENSION(INhkl,INhkl),INTENT(IN) :: CUgMat
    INTEGER(IKIND),INTENT(IN) :: IMinWeakBeams, IMinStrongBeams
    INTEGER(IKIND),DIMENSION(INhkl),INTENT(OUT) :: IStrongBeamList,IWeakBeamList
    INTEGER(IKIND),INTENT(OUT) :: nBeams,nWeakBeams,IErr
    INTEGER(IKIND) :: ind,jnd
    INTEGER(IKIND),DIMENSION(:) :: IStrong(INhkl),IWeak(INhkl)
    REAL(RKIND) :: RMaxSg,RMinPertStrong,RMinPertWeak
    REAL(RKIND),DIMENSION(:) :: RPertStrength0(INhkl)

    !----------------------------------------------------------------------------
    ! strong beams
    !----------------------------------------------------------------------------

    ! Use Sg and perturbation strength to define strong beams
    ! PerturbationStrength Eq. 8 Zuo Ultramicroscopy 57 (1995) 375, |Ug/2KSg|
    ! Here use |Ug/Sg| since 2K is a constant
    ! NB RPertStrength0 is an array of perturbation strengths for all reflections
    RPertStrength0 = ABS(CUgMat(:,1)/(RDevPara))
    ! 000 beam is NaN otherwise, always included by making it a large number
    RPertStrength0(1) = 1000.0

    ! NB IStrong is an array listing the strong beams (1=Strong, 0=Not strong)
    IStrong=0_IKIND
    ! start with a small deviation parameter limit
    RMaxSg = 0.005
    RMinPertStrong=0.0025/RMaxSg ! Gives additional beams based on perturbation strength

    ! main calculation
    ! now increase RMaxSg until we have enough strong beams
    DO WHILE (SUM(IStrong).LT.IMinStrongBeams)
      WHERE (ABS(RDevPara).LT.RMaxSg.OR.RPertStrength0.GE.RMinPertStrong)
        IStrong=1_IKIND
      END WHERE
      RMaxSg=RMaxSg+0.005
      ! RMinPertStrong=0.0025/RMaxSg
    END DO

    ! give the strong beams a number in IStrongBeamList
    IStrongBeamList=0_IKIND
    ind=1_IKIND
    DO jnd=1,INhkl
      IF (IStrong(jnd).EQ.1) THEN
        IStrongBeamList(ind)=jnd
        ind=ind+1
      END IF
    END DO
    !The no. of strong beams gives the dimension of the Bloch wave problem
    nBeams=ind-1  
    
    CALL message(LXL,dbg7,"Strong Beam List",IStrongBeamList)
    CALL message(LXL,dbg7,"Sg limit for strong beams = ",RMaxSg)
    CALL message(LXL,dbg7,"Smallest strong perturbation strength = ",RMinPertStrong)
    IF(SUM(IStrong)+IMinWeakBeams.GT.INhkl) IErr = 1
    IF(l_alert(IErr,"StrongAndWeakBeamsDetermination",&
          "Insufficient reflections to accommodate all Strong and Weak Beams")) RETURN
    
    !----------------------------------------------------------------------------
    ! weak beams
    !----------------------------------------------------------------------------

    ! Decrease perturbation strength until we have enough weak beams
    ! NB IWeak is an array listing the weak beams (1=Weak, 0=Not weak)
    IWeak=0_IKIND
    RMinPertWeak=0.9*RMinPertStrong
    DO WHILE (SUM(IWeak).LT.IMinWeakBeams)
      WHERE (RPertStrength0.GE.RMinPertWeak.AND.IStrong.NE.1_IKIND)
        IWeak=1
      END WHERE
      RMinPertWeak=0.9*RMinPertWeak
    END DO

    CALL message(LXL,dbg7,"weak beams",SUM(IWeak))
    CALL message(LXL,dbg7,"Smallest weak perturbation strength = ",RMinPertWeak)

    ! give the weak beams a number in IWeakBeamList
    IWeakBeamList=0_IKIND
    ind=1_IKIND
    DO jnd=1,INhkl
      IF (IWeak(jnd).EQ.1) THEN
        IWeakBeamList(ind)=jnd
        ind=ind+1
      END IF
    END DO
    nWeakBeams=ind-1

    CALL message(LXL,dbg7,"Weak Beam List",IWeakBeamList)

  END SUBROUTINE StrongAndWeakBeamsDetermination

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Returns eigenvalues and eigenvectors of matrix.
  !!
  !! Major-Authors: Rudo Roemer (2014), Richard Beanland (2016)
  !!
  SUBROUTINE EigenSpectrum(IMatrixDimension, MatrixToBeDiagonalised, EigenValues,&
                    EigenVectors, IErr)

    ! accesses procedure ZGEEV
    USE MyNumbers
    USE message_mod
    USE MyMPI

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(IN) :: IMatrixDimension
    COMPLEX(RKIND),INTENT(IN) :: MatrixToBeDiagonalised(IMatrixDimension,IMatrixDimension)
    COMPLEX(RKIND),INTENT(OUT) :: EigenValues(IMatrixDimension),&
          EigenVectors(IMatrixDimension,IMatrixDimension)
    INTEGER(IKIND),INTENT(OUT) :: IErr
    INTEGER(IKIND) :: WorkSpaceDimension
    ! dummy vector outputs used while finding respective eigenvectors/values
    COMPLEX(CKIND),DIMENSION(:), ALLOCATABLE :: CWorkSpace 
    REAL(RKIND), DIMENSION(:), ALLOCATABLE :: WorkSpace
    EXTERNAL ZGEEV

    ! find optimum size of arrays
    WorkSpaceDimension=1
    ALLOCATE(CWorkSpace(WorkSpaceDimension),STAT = IErr)
    IF(l_alert(IErr,"EigenSpectrum","allocate CWorkSpace")) RETURN
    ALLOCATE(WorkSpace(2*IMatrixDimension),STAT = IErr)
    IF(l_alert(IErr,"EigenSpectrum","allocate WorkSpace")) RETURN

    WorkSpaceDimension=-1

    CALL ZGEEV('N','V', IMatrixDimension, MatrixToBeDiagonalised, IMatrixDimension,&
         EigenValues, 0,1, EigenVectors,IMatrixDimension, &
         CWorkSpace, WorkSpaceDimension, WorkSpace, IErr )
    IF(l_alert(IErr,"EigenSpectrum","ZGEEV()")) RETURN

    WorkSpaceDimension = INT(CWorkSpace(1))

    ! REALLOCATE necessary memory
    DEALLOCATE(CWorkSpace,STAT=IErr)
    IF(l_alert(IErr,"EigenSpectrum","deallocate CWorkSpace")) RETURN
    ALLOCATE(CWorkSpace(WorkSpaceDimension),STAT = IErr)
    IF(l_alert(IErr,"EigenSpectrum","allocate CWorkSpace")) RETURN

    ! do the actual call to get the spectrum
    CALL ZGEEV('N','V', IMatrixDimension, MatrixToBeDiagonalised, IMatrixDimension,&
         EigenValues, 0,1, EigenVectors,IMatrixDimension, &
         CWorkSpace, WorkSpaceDimension, WorkSpace, IErr )
    IF(l_alert(IErr,"EigenSpectrum","ZGEEV()")) RETURN

    DEALLOCATE(CWorkSpace,WorkSpace,STAT = IErr)
    IF(l_alert(IErr,"EigenSpectrum","deallocate")) RETURN

    RETURN

  END SUBROUTINE EigenSpectrum

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Invert an M*M Complex Matrix
  !!
  !! Major-Authors: Rudo Roemer (2014), Richard Beanland (2016)
  !!
  SUBROUTINE INVERT(MatrixSize,Matrix,InvertedMatrix,IErr)  

    ! accesses procedure ZGETRF
    USE MyNumbers
    USE message_mod
    USE MyMPI
    
    IMPLICIT NONE
    
    INTEGER(IKIND),INTENT(IN) :: MatrixSize
    COMPLEX(CKIND),INTENT(INOUT) :: Matrix(MatrixSize,MatrixSize) ! destroyed in process
    COMPLEX(CKIND),INTENT(OUT) :: InvertedMatrix(1:MatrixSize,1:MatrixSize)
    INTEGER(IKIND),INTENT(OUT) :: IErr

    INTEGER :: LWORK, INFO, I
    INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
    COMPLEX(CKIND), DIMENSION(:), ALLOCATABLE :: WORK
    
    ALLOCATE(IPIV(MatrixSize),STAT=IErr)
    IF(l_alert(IErr,"EigenSpectrum","allocate IPIV")) RETURN
    
    CALL ZGETRF(MatrixSize,MatrixSize,Matrix,MatrixSize,IPIV,IErr)
    IF(l_alert(IErr,"EigenSpectrum","ZGETRF()")) RETURN

    LWORK = MatrixSize*MatrixSize
    ALLOCATE(WORK(LWORK),STAT=IErr)   
    IF(l_alert(IErr,"EigenSpectrum","WORK")) RETURN
    
    CALL ZGETRI(MatrixSize,Matrix,MatrixSize,IPIV,WORK,LWORK,IErr)
    IF(l_alert(IErr,"EigenSpectrum","ZGETRI()")) RETURN

    DEALLOCATE(IPIV,WORK,STAT=IErr)
    IF(l_alert(IErr,"EigenSpectrum","deallocate IPIV")) RETURN

    InvertedMatrix = Matrix  
    RETURN

  END SUBROUTINE INVERT

END MODULE bloch_mod
