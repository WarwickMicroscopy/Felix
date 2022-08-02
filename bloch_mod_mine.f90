!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: :VERSION: RB_coord / 1.15 /
! Date:    :DATE: 16-01-2019
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build:   :BUILD: Mode F: test different lattice types" 
! Author:  :AUTHOR: r.beanland
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
!! Module-description: Solely for BlochCoefficientCalculation procedure, calls all other subroutines later defined.
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
                    IFirstPixelToCalculate,IErr)

    ! Imports modules
	! Accesses procedure ZGEMM
    USE MyNumbers
    USE MyMPI
    USE message_mod
    
    USE test_koch_mod
  
    ! Globals - output
    USE RPara, ONLY : RIndividualReflections ! RIndividualReflections( LACBED_ID, thickness_ID, local_pixel_ID )
    USE CPara, ONLY : CAmplitudeandPhase
    USE IPara, ONLY : IPixelComputed

    ! Globals - input  
    USE CPara, ONLY : CUgMat
    USE RPara, ONLY : RDeltaK,RDeltaThickness,RInitialThickness,RNormDirM,RgDotNorm,RgPool,&
                      RgPoolMag,Rhkl
    USE IPara, ONLY : IHolzFLAG,IMinStrongBeams,IMinWeakBeams,&
                      INoOfLacbedPatterns,IPixelCount,IThicknessCount,INhkl,&
                      IOutputReflections,IBlochMethodFLAG
    USE BlochPara, ONLY : RBigK            
    USE SPARA, ONLY : SPrintString
    
    IMPLICIT NONE
	
	! Defining local variables
    
    INTEGER(IKIND),INTENT(IN) :: IYPixelIndex,IXPixelIndex,IPixelNumber,&
          IFirstPixelToCalculate
    INTEGER(IKIND),INTENT(OUT) :: IErr
    
    COMPLEX(CKIND),ALLOCATABLE :: CBeamProjectionMatrix(:,:),&
          CDummyBeamMatrix(:,:),CUgSgMatrix(:,:),CEigenVectors(:,:),CEigenValues(:),&
          CInvertedEigenVectors(:,:),CAlphaWeightingCoefficients(:),&
          CEigenValueDependentTerms(:,:)
    COMPLEX(CKIND) :: CFullWaveFunctions(INhkl)
    REAL(RKIND) :: RFullWaveIntensity(INhkl),RDevPara(INhkl),&
          RTiltedK(ITHREE)
    INTEGER(IKIND) :: IStrongBeamList(INhkl),IWeakBeamList(INhkl),&
          nBeams,nWeakBeams
    INTEGER(IKIND) :: ind,knd,pnd,IThickness,IThicknessIndex,ILowerLimit,&
          IUpperLimit       
    REAL(RKIND) :: RThickness,RKn,Rk0(3),RkPrime(3)
    COMPLEX(CKIND) sumC,sumD
    COMPLEX(CKIND), DIMENSION(:,:), ALLOCATABLE :: CBeamTranspose,CUgMatPartial,CDummyEigenVectors
    CHARACTER*40 surname
    CHARACTER*100 SindString,SjndString,SPixelCount,SnBeams,SWeakBeamIndex

    ! Variables used for koch spence method development
    COMPLEX(CKIND),ALLOCATABLE :: CDiagonalSgMatrix(:,:), COffDiagonalSgMatrix(:,:)
    COMPLEX(CKIND) :: CScatteringElement
    INTEGER(IKIND) :: ScatterMatrixRow
     
    IErr=0
    ! we are inside the mask
    IPixelComputed= IPixelComputed + 1

    ! TiltedK is the vector of the incoming tilted beam
    ! in units of (1/A), in the microscope ref frame(NB exp(i*k.r), physics convention)
	! Are these the wrong way round?
    ! x-position in k-space
    RTiltedK(1)= (REAL(IYPixelIndex,RKIND)-REAL(IPixelCount,RKIND)-0.5_RKIND)*RDeltaK
    ! y-position in k-space
    RTiltedK(2)= (REAL(IXPixelIndex,RKIND)-REAL(IPixelCount,RKIND)-0.5_RKIND)*RDeltaK 
	! Finding z component in k space(?)
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

    ! Select only those beams where the Ewald sphere is close to the
    ! reciprocal lattice, i.e. within RBSMaxDeviationPara
    CALL StrongAndWeakBeamsDetermination(INhkl,IMinWeakBeams,&
                    IMinStrongBeams,RDevPara,CUgMat,&
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
    IF(l_alert(IErr,"BlochCoefficientCalculation","allocations")) RETURN

    ! Allocations used for koch spence method development
    IF(IBlochMethodFLAG.EQ.1) THEN
      ALLOCATE( CDiagonalSgMatrix(nBeams,nBeams), STAT=IErr )
      IF(l_alert(IErr,"BlochCoefficientCalculation","allocate CDiagonalSgMatrix")) RETURN
      ALLOCATE( COffDiagonalSgMatrix(nBeams,nBeams), STAT=IErr )
      IF(l_alert(IErr,"BlochCoefficientCalculation","allocate COffDiagonalSgMatrix")) RETURN
    END IF

    ! Compute the effective Ug matrix by selecting only those beams
    ! for which IStrongBeamList has an entry
    CBeamProjectionMatrix= CZERO
    DO knd=1,nBeams
      CBeamProjectionMatrix(knd,IStrongBeamList(knd))=CONE
    ENDDO


    CUgSgMatrix = CZERO
    CBeamTranspose=TRANSPOSE(CBeamProjectionMatrix)
    ! Reduce the matrix to just include strong beams using some nifty matrix multiplication
    ! CUgMatPartial = CUgMat * CBeamTranspose
    CALL ZGEMM('N','N',INhkl,nBeams,INhkl,CONE,CUgMat, &
              INhkl,CBeamTranspose,INhkl,CZERO,CUgMatPartial,INhkl)
    ! CUgSgMatrix = CBeamProjectionMatrix * CUgMatPartial
    CALL ZGEMM('N','N',nBeams,nBeams,INhkl,CONE,CBeamProjectionMatrix, &
              nBeams,CUgMatPartial,INhkl,CZERO,CUgSgMatrix,nBeams)

    !--------------------------------------------------------------------
    ! Higher order Laue zones and weak beams
    !--------------------------------------------------------------------
	
	!!!!!! Should be able to get rid of this, as not considering higher order laue zones.
	
    IF (IHolzFLAG.EQ.1) THEN!We are considering higher order Laue Zones !?? suspect this is non-functional
      DO ind=1,nBeams
        CUgSgMatrix(ind,ind) = CUgSgMatrix(ind,ind) + TWO*RBigK*RDevPara(IStrongBeamList(ind))
      ENDDO
      DO knd =1,nBeams ! Columns
        DO ind = 1,nBeams ! Rows
          CUgSgMatrix(knd,ind) = CUgSgMatrix(knd,ind) / &
                (SQRT(1+RgDotNorm(IStrongBeamList(knd))/RKn)*&
                SQRT(1+RgDotNorm(IStrongBeamList(ind))/RKn))
        END DO
      END DO
      CUgSgMatrix = (TWOPI**2)*CUgSgMatrix/(TWO*RBigK)
    ELSE!ZOLZ only
      ! replace the diagonal parts with strong beam deviation parameters
      DO ind=1,nBeams
        CUgSgMatrix(ind,ind) = TWO*RBigK*RDevPara(IStrongBeamList(ind))/(TWOPI*TWOPI)
      ENDDO
      ! add the weak beams perturbatively for the 1st column (sumC) and
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
      !The 4pi^2 is a result of using h, not hbar, in the conversion from VG(ij) to Ug(ij).  Needs to be taken out of the weak beam calculation too 
      !Divide by 2K so off-diagonal elementa are Ug/2K, diagonal elements are Sg, Spence's (1990) 'Structure matrix'
      CUgSgMatrix = TWOPI*TWOPI*CUgSgMatrix/(TWO*RBigK)
    END IF
   
    !--------------------------------------------------------------------
    ! Diagonalize the UgMatEffective
    !--------------------------------------------------------------------

    ! If koch method - Split CUgSgMatrix into diagonal and off diagonal to speed convergence
	! Can maybe get rid of, as it is koch and spence method?
    IF(IBlochMethodFLAG.EQ.1) THEN
      COffDiagonalSgMatrix = CUgSgMatrix
      CDiagonalSgMatrix = CZERO
      DO ind = 1,SIZE(CUgSgMatrix,2)
        CDiagonalSgMatrix(ind,ind) = CUgSgMatrix(ind,ind)      
        COffDiagonalSgMatrix(ind,ind) = CZERO
      END DO
    END IF

    CALL EigenSpectrum(nBeams,CUgSgMatrix,CEigenValues(:),CEigenVectors(:,:),IErr)
    IF(l_alert(IErr,"BlochCoefficientCalculation","EigenSpectrum()")) RETURN
    ! NB destroys CUgSgMatrix
	! Is this bit needed as we aren't considering higher order laue zones.
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
    ! Fill RIndividualReflections( LACBED_ID , thickness_ID, local_pixel_ID ) 
    !--------------------------------------------------------------------
   
    ! Calculate intensities for different specimen thicknesses
    !?? Do different g-vectors have different effective thicknesses??
    DO IThicknessIndex=1,IThicknessCount,1

      RThickness = RInitialThickness + REAL((IThicknessIndex-1),RKIND)*RDeltaThickness
      IThickness = NINT(RThickness,IKIND)

      CALL CreateWaveFunctions(RThickness,RFullWaveIntensity,CFullWaveFunctions,&
                    INhkl,nBeams,IStrongBeamList,CEigenVectors,CInvertedEigenVectors,CEigenValues,IErr)
      IF(l_alert(IErr,"BlochCoefficientCalculation","CreateWaveFunctions")) RETURN

      !--------------------------------------------------------------------
      ! Optional - test koch spence prototype method
	  ! Might be able to get rid of entire IF statement.
      !--------------------------------------------------------------------
      IF(IBlochMethodFLAG.EQ.1) THEN
        RThickness = RThickness / 1000! for koch development to speed convergence
        CALL message('-----------------------------------------------------------------------')
        CALL message('-----------------------------------------------------------------------')
        CALL message('RThickness divided by 1000 to help koch series convergence')
        CALL message('CFullWaveFunctions(1:4)',CFullWaveFunctions(1:4))
        CALL CalculateElementS( CMPLX(ZERO,RThickness,CKIND), COffDiagonalSgMatrix, CDiagonalSgMatrix, 1, 1, 5, CScatteringElement )
        CALL CalculateElementS( CMPLX(ZERO,RThickness,CKIND), COffDiagonalSgMatrix, CDiagonalSgMatrix, 2, 1, 5, CScatteringElement )
        CALL CalculateElementS( CMPLX(ZERO,RThickness,CKIND), COffDiagonalSgMatrix, CDiagonalSgMatrix, 3, 1, 5, CScatteringElement )
        CALL CalculateElementS( CMPLX(ZERO,RThickness,CKIND), COffDiagonalSgMatrix, CDiagonalSgMatrix, 4, 1, 5, CScatteringElement )

!        CALL message('-----------------------------------------------------------------------')
!        CALL message('Below shows the wavefunction matrix from diagonalisation and then the koch method')
!        CALL message('The thickness has been scaled by 1/1000 to speed convergence')
!        CALL message('The matrices are also ordered differently')
!        CALL message('and in the diagonalisation method some entries are negligable and left as zero')
!        CALL message('(diagonlisation) wavefunction pixel values for this thickness and this core')
!        DO ScatterMatrixRow = 1,nBeams
!          CALL message('',CFullWaveFunctions(ScatterMatrixRow))
!        END DO
!        CALL message('(koch series) wavefunction pixel values for this thickness and this core') 
!        DO ScatterMatrixRow = 1,nBeams
!          CALL CalculateElementS( CMPLX(ZERO,RThickness,CKIND), COffDiagonalSgMatrix, &
!                CDiagonalSgMatrix, ScatterMatrixRow, 1, 4, CScatteringElement )
!          CALL message('',CScatteringElement)
!        END DO

        !test for a single pixel 
!        IF(IYPixelIndex.EQ.10.AND.IXPixelIndex.EQ.10.AND.IThicknessIndex.EQ.2) THEN 
!          CALL message('debug reached single pixel thickness test')
!          ! scale RThickness to help convergence issues
!          RThickness = RThickness / 1000

!          CALL CreateWaveFunctions(RThickness,RFullWaveIntensity,CFullWaveFunctions,&
!                        INhkl,nBeams,IStrongBeamList,CEigenVectors,CEigenValues,IErr)
!          IF(l_alert(IErr,"BlochCoefficientCalculation","CreateWaveFunctions")) RETURN
!          CALL message('CFullWaveFunctions(1:4)',CFullWaveFunctions(1:4))

!          ! calculate scattering matrix using koch spence method
!          !CALL GetCombinations()
!          CALL CalculateElementS( CMPLX(ZERO,RThickness,CKIND), COffDiagonalSgMatrix, CDiagonalSgMatrix, 1, 1, 4, CScatteringElement )

!          ! for debugging end felix here
!          CALL message('Testing koch series method, so terminate felix here.')
!          CALL message('-----------------------------------------------------------------------')
!          CALL SLEEP(1)
!          IErr = 1
!          RETURN
!        END IF

    ! deallocations used for koch spence method development
      DEALLOCATE( CDiagonalSgMatrix, COffDiagonalSgMatrix, STAT=IErr )
      IF(l_alert(IErr,"BlochCoefficientCalculation","deallocating arrays")) RETURN

      ! Testing koch series so terminate felix here
      CALL message('-----------------------------------------------------------------------')
      CALL message('Testing koch series method, so terminate felix here.')
      CALL message('-----------------------------------------------------------------------')
      CALL message('-----------------------------------------------------------------------')
      CALL message('-----------------------------------------------------------------------')
      CALL SLEEP(1)
      IErr = 1
      RETURN
    END IF

      ! Collect Intensities from all thickness for later writing
      ! Output without felix.hkl input disabled
      !IF(IHKLSelectFLAG.EQ.0) THEN ! we are not using hkl list from felix.hkl
      !  RIndividualReflections(1:INoOfLacbedPatterns,IThicknessIndex,&
      !          (IPixelNumber-IFirstPixelToCalculate)+1) = &
      !          RFullWaveIntensity(1:INoOfLacbedPatterns)
      !ELSE ! we are using hkl list from felix.hkl
        DO ind = 1,INoOfLacbedPatterns
          RIndividualReflections(ind,IThicknessIndex,&
                  (IPixelNumber-IFirstPixelToCalculate)+1) = &
                  RFullWaveIntensity(IOutputReflections(ind))
        END DO
      !END IF
    END DO


    ! DEALLOCATE eigen problem memory
    DEALLOCATE(CUgSgMatrix,CBeamTranspose, CUgMatPartial, &
         CInvertedEigenVectors, CAlphaWeightingCoefficients, &
         CEigenValues,CEigenVectors,CEigenValueDependentTerms, &
         CBeamProjectionMatrix, CDummyBeamMatrix,STAT=IErr)
    IF(l_alert(IErr,"BlochCoefficientCalculation","deallocating arrays")) RETURN
    
  END SUBROUTINE BlochCoefficientCalculation

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Calculates diffracted intensity for a specific thickness from wavefunctions. 
  !! 
  !! Also finds the excitation coefficients for bloch waves for real space wavefunction of electron.
  !! Eigenvectors and values have been solved for, as matrix has been diagonalised(see subroutine eigenspectrum)
  !! 
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE CreateWaveFunctions(RThickness,RFullWaveIntensity,CFullWaveFunctions,&
                    INhkl,nBeams,IStrongBeamList,CEigenVectors,CInvertedEigenVectors,CEigenValues,IErr)

	!Importing modules
	
    USE MyNumbers
    USE MyMPI
    USE message_mod

    IMPLICIT NONE
    
	! Defining local variables
	
    REAL(RKIND),INTENT(IN) :: RThickness
    REAL(RKIND),INTENT(OUT) :: RFullWaveIntensity(INhkl)
    COMPLEX(CKIND),INTENT(OUT) :: CFullWaveFunctions(INhkl) 
    INTEGER(IKIND),INTENT(IN) :: INhkl,nBeams,IStrongBeamList(INhkl)
    COMPLEX(CKIND),INTENT(IN) :: CEigenVectors(nBeams,nBeams),CInvertedEigenVectors(nBeams,nBeams),CEigenValues(nBeams)
    INTEGER(IKIND),INTENT(OUT) :: IErr 
    REAL(RKIND) :: RWaveIntensity(nBeams)
    COMPLEX(CKIND) :: CPsi0(nBeams),CAlphaWeightingCoefficients(nBeams),&
          CWaveFunctions(nBeams),CEigenValueDependentTerms(nBeams,nBeams)
    INTEGER(IKIND) :: ind,jnd,knd,hnd,ifullind,iuniind,gnd,ichnk
    
    IErr=0
    ! The top surface boundary conditions
    CPsi0 = CZERO ! All diffracted beams are zero
    CPsi0(1) = CONE ! The 000 beam has unit amplitude

    ! Put in the thickness
    ! From EQ 6.32 in Kirkland Advance Computing in EM
    CAlphaWeightingCoefficients = MATMUL(CInvertedEigenVectors(1:nBeams,1:nBeams),CPsi0) 
    CEigenValueDependentTerms= CZERO
    DO hnd=1,nBeams     ! This is a diagonal matrix
      CEigenValueDependentTerms(hnd,hnd) = &
            EXP(CIMAGONE*CMPLX(RThickness,ZERO,CKIND)*CEigenValues(hnd)) 
    ENDDO
    ! The diffracted intensity for each beam
    ! EQ 6.35 in Kirkland Advance Computing in EM
    ! C-1*C*alpha = alpha 
    CWaveFunctions(:) = MATMUL( &
          MATMUL(CEigenVectors(1:nBeams,1:nBeams),CEigenValueDependentTerms), & 
          CAlphaWeightingCoefficients(:) )

    !?? Possible small time saving here by only calculating the (tens of) output
    !?? Reflections rather than all strong beams (hundreds)
	
	! Calculates full wave intensity from amplitude, complex conjugate * original wavefunction.
    DO hnd=1,nBeams
       RWaveIntensity(hnd)=CONJG(CWaveFunctions(hnd)) * CWaveFunctions(hnd)
    ENDDO  
	! Whats this bit for? I thought the above bit has already found the intensity(?)
    CFullWaveFunctions=CZERO
    RFullWaveIntensity=ZERO
    DO knd=1,nBeams
       CFullWaveFunctions(IStrongBeamList(knd))=CWaveFunctions(knd)
       RFullWaveIntensity(IStrongBeamList(knd))=RWaveIntensity(knd)
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
                    IMinStrongBeams,RDevPara,CUgMat,&
                    IStrongBeamList,IWeakBeamList,nBeams,nWeakBeams,IErr)
    
    ! Select only those beams where the Ewald sphere is close to the
    ! reciprocal lattice, i.e. within RBSMaxDeviationPara
	
	
	! Check lab book for diagram of Ewald Sphere.
	
	! Importing modules
    USE MyNumbers
    USE MyMPI
    USE message_mod 
	
	! Defining local variables
	
    INTEGER(IKIND),INTENT(IN) :: INhkl
    REAL(RKIND),DIMENSION(INhkl),INTENT(IN) :: RDevPara
    COMPLEX(CKIND),DIMENSION(INhkl,INhkl),INTENT(IN) :: CUgMat
    INTEGER(IKIND),INTENT(IN) :: IMinWeakBeams, IMinStrongBeams
    INTEGER(IKIND),DIMENSION(INhkl),INTENT(OUT) :: IStrongBeamList,IWeakBeamList
    INTEGER(IKIND),INTENT(OUT) :: nBeams,nWeakBeams,IErr
    INTEGER(IKIND) :: ind,jnd
    INTEGER(IKIND),DIMENSION(:) :: IStrong(INhkl),IWeak(INhkl)
    REAL(RKIND) :: RMaxSg,RMinPertStrong,RMinPertWeak ! RMaxSg is the maximum deviation parameter from lattice points.
    REAL(RKIND),DIMENSION(:) :: RPertStrength0(INhkl)

    !----------------------------------------------------------------------------
    ! Strong beams
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

    ! Main calculation
    ! Now increase RMaxSg until we have enough strong beams
    DO WHILE (SUM(IStrong).LT.IMinStrongBeams)
      WHERE (ABS(RDevPara).LT.RMaxSg.OR.RPertStrength0.GE.RMinPertStrong)
        IStrong=1_IKIND
      END WHERE
      RMaxSg=RMaxSg+0.005
    END DO

    ! Give the strong beams a number in IStrongBeamList
    IStrongBeamList=0_IKIND
    ind=1_IKIND
    DO jnd=1,INhkl
      IF (IStrong(jnd).EQ.1) THEN
        IStrongBeamList(ind)=jnd
        ind=ind+1
      END IF
    END DO
    !The no. of strong beams gives the dimension of the Bloch wave problem. See later subroutines for solutions.
    nBeams=ind-1  

	! See message_mod for details on message.
	
    CALL message(LXL,dbg7,"Strong Beam List",IStrongBeamList)
    CALL message(LXL,dbg7,"Sg limit for strong beams = ",RMaxSg)
    CALL message(LXL,dbg7,"Smallest strong perturbation strength = ",RMinPertStrong)
    IF(SUM(IStrong)+IMinWeakBeams.GT.INhkl) IErr = 1
    IF(l_alert(IErr,"StrongAndWeakBeamsDetermination",&
          "Insufficient reflections to accommodate all Strong and Weak Beams")) RETURN
    
    !----------------------------------------------------------------------------
    ! Weak beams
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

    ! Give the weak beams a number in IWeakBeamList
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
  !! Procedure-description: Returns eigenvalues and eigenvectors of matrix. For complete electron wavefunction in crystal.
  !!
  !! Major-Authors: Rudo Roemer (2014), Richard Beanland (2016)
  !!
  SUBROUTINE EigenSpectrum(IMatrixDimension, MatrixToBeDiagonalised, EigenValues,&
                    EigenVectors, IErr)

    ! Importing modules
	! Allows access to ZGEEV
    USE MyNumbers
    USE message_mod
    USE MyMPI

    IMPLICIT NONE
	
	! Defining local variables.
	
    INTEGER(IKIND),INTENT(IN) :: IMatrixDimension
    COMPLEX(RKIND),INTENT(IN) :: MatrixToBeDiagonalised(IMatrixDimension,IMatrixDimension)
    COMPLEX(RKIND),INTENT(OUT) :: EigenValues(IMatrixDimension),&
          EigenVectors(IMatrixDimension,IMatrixDimension)
    INTEGER(IKIND),INTENT(OUT) :: IErr ! Output
    INTEGER(IKIND) :: WorkSpaceDimension
    ! Dummy vector outputs used while finding respective eigenvectors/values
    COMPLEX(CKIND),DIMENSION(:), ALLOCATABLE :: CWorkSpace 
    REAL(RKIND), DIMENSION(:), ALLOCATABLE :: WorkSpace
    EXTERNAL ZGEEV ! This module calculates the eigenvalues and eigenvectors of matrices.

    ! Find optimum size of arrays
    WorkSpaceDimension=1
    ALLOCATE(CWorkSpace(WorkSpaceDimension),STAT = IErr)
    IF(l_alert(IErr,"EigenSpectrum","allocate CWorkSpace")) RETURN
    ALLOCATE(WorkSpace(2*IMatrixDimension),STAT = IErr)
    IF(l_alert(IErr,"EigenSpectrum","allocate WorkSpace")) RETURN

    WorkSpaceDimension=-1
	! Calculating eigenvalues and eigenvectors of input matrix for diagonalisation.
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

    ! Do the actual call to get the spectrum of eigenvalues.
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
  !! Procedure-description: Inverts an M*M Complex Matrix
  !!
  !! Major-Authors: Rudo Roemer (2014), Richard Beanland (2016)
  !!
  SUBROUTINE INVERT(MatrixSize,Matrix,InvertedMatrix,IErr)  
	! Importing modules
    ! Allows access to ZGETRF and ZGETRI procedures for inversion.
    USE MyNumbers
    USE message_mod
    USE MyMPI
    
    IMPLICIT NONE
	
	! Defining local variables
    
    INTEGER(IKIND),INTENT(IN) :: MatrixSize
    COMPLEX(CKIND),INTENT(INOUT) :: Matrix(MatrixSize,MatrixSize) ! destroyed in process
    COMPLEX(CKIND),INTENT(OUT) :: InvertedMatrix(1:MatrixSize,1:MatrixSize)
    INTEGER(IKIND),INTENT(OUT) :: IErr

    INTEGER :: LWORK, INFO, I
    INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV ! Integer array containing pivot indices that define permutation(P) matrix, as solve matrix equation using LU decomposition with P matrix.
    COMPLEX(CKIND), DIMENSION(:), ALLOCATABLE :: WORK
    
    ALLOCATE(IPIV(MatrixSize),STAT=IErr) ! STAT has been assigned IErr, which contains number of error message at run time, if allocate statement passes then IErr = 0
    IF(l_alert(IErr,"EigenSpectrum","allocate IPIV")) RETURN
    
    CALL ZGETRF(MatrixSize,MatrixSize,Matrix,MatrixSize,IPIV,IErr) ! Factors matrix using gaussian elimination(LU decomposition)
    IF(l_alert(IErr,"EigenSpectrum","ZGETRF()")) RETURN

    LWORK = MatrixSize*MatrixSize
    ALLOCATE(WORK(LWORK),STAT=IErr)   
    IF(l_alert(IErr,"EigenSpectrum","WORK")) RETURN
    
    CALL ZGETRI(MatrixSize,Matrix,MatrixSize,IPIV,WORK,LWORK,IErr) ! Computes inverse of LU matrix found from ZGETRF
    IF(l_alert(IErr,"EigenSpectrum","ZGETRI()")) RETURN

    DEALLOCATE(IPIV,WORK,STAT=IErr)
    IF(l_alert(IErr,"EigenSpectrum","deallocate IPIV")) RETURN

    InvertedMatrix = Matrix  
    RETURN

  END SUBROUTINE INVERT
  
    SUBROUTINE angle_correction(IErr)

		!! Subroutine-Description:
		!! Accounts for the possibility that the incident electron beam is not parallel to the surface normal of the specimen.
		!! Proceeds to calculate the resulting wavefunction and the diffracted intensities from the beam

		!! Finds the structure matrix that accounts for absorption and the angle of the beam
		!! Finds the eigenvalues and eigenvectors of structure matrix
		!! Forms the scattering matrix
		!! Finds diffracted intensities from scattering matrix.
		!! Based off "Structure refinement using precession electron diffraction tomography
		!! and dynamical diffraction: theory and implementation" By Lukas Palatinus(2015)



		! Import modules and global variables
		
		USE MyNumbers
		USE MyMPI
		USE BlochPara ONLY : RTiltedK
		USE RPARA, ONLY : RNormDirM, RgDotNorm, RMeanInnerPotential, RgMatrix			
		USE IPara ONLY : nBeams ! nBeams after strong and weak beam determination
		USE CPara ONLY : CUgMat ! This gives Ug matrix with absorption
		
		! Imports external functions
		
		EXTERNAL NORM2
		EXTERNAL IFFT2

		IMPLICIT NONE

		! Define local variables
		
		REAL :: RKn, REigenvalues(nBeams), RThickness, &
				REigenvalue_matrix(nBeams, nBeams), RElement(nBeams), &
				RMmatrix(nBeams, nBeams), RIntensity_vector(nBeams)
		
		COMPLEX :: CStructureMatrix(nBeams, nBeams), CDiagonal_element(nBeams) &
				   CElement_off(nBeams), CEigenvector_matrix(nBeams, nBeams) &
				   CInitial_wavefunction(nBeams), CInverted_vectors(nBeams, nBeams), &
				   CWeighting_Coefficients(nBeams), CScatter_matrix(nBeams, nBeams), &
				   CAlt_scatter_matrix(nBeams, nBeams), CFinal_wavefunction_z(nBeams), &
				   CFinal_wavefunction(nBeams)

		INTEGER :: ind, jnd
		
		! Calculates components of g vectors and K vector along surface normal vector
		
		! First K vector
		RKn = DOT_PRODUCT(RTiltedK, RNormDirM)
		IF(l_alert(IErr,"angle_correction","components, K component")) CALL abort
		
		! gVectors dotted with normal vector are contained in RgDotNorm, as array
		
		
		
		! Next calculate the new structure matrix from Ug matrix, accounting for absorption
		! Finding off and on diagonal elements of structure matrix
		! First loop iterates over row
		! Second loop over columns
		K = ((NORM2(RTiltedK))**2 + RMeanInnerPotential)**0.5 ! This is mod of K, in sample, RMeanInnerPotential is U0
		DO ind = 1, nBeams
			Kg = (NORM2(K + RgMatrix(3, ind, :)))**0.5 ! This is mod of K + g
			DO jnd = 1, nBeams
				IF ind = jnd THEN
					! On diagonal elements
					CDiagonal_element = (K**2 - Kg**2)/(1 + RgDotNorm(ind)/RKn)**0.5
					CStructureMatrix(ind, ind) = CDiagonal_element
				ELSE
					! Off diagonal elements
					CElement_off = (CUgMat(ind, jnd))/(((1 + (RgDotNorm(ind))/(RKn))**0.5) * ((1 + (RgDotNorm(jnd))/(RKn))**0.5))
					CStructureMatrix(ind, jnd) = CElement_off
			END DO
		END Do
		
		IF(l_alert(IErr,"angle_correction","structure_matrix")) CALL abort
		
		
		! Now find eigenvalues and eigenvectors of new structure matrix, so as to create scatter matrices
		! Matrices are hermitian, so eigenvalues are real
		! All matrices are square
		
		CALL EigenSpectrum(nBeams, CStructureMatrix)
		CEigenvector_matrix = VR ! VR is output from ZGEEV(in EigenSpectrum) containing eigenvectors along columns of matrix
		IF(l_alert(IErr,"angle_correction","eigen_structure_matrix, eigenvectors")) CALL abort
		REigenvalues = W ! W is output from ZGEEV(in EigenSpectrum) containing eigenvalues as an array
		IF(l_alert(IErr,"angle_correction","eigen_structure_matrix, eigenvalues")) CALL abort
		
		! Form matrix of eigenvalues, with exponent, thickness, K component, along diagonal
		DO ind = 1, nBeams
			REigenvalue_matrix(ind, ind) = EXP((((TWOPI * CIMAGONE * RThickness)/ (2 * RKn))) * REigenvalues(ind))
		END DO
		
		
		! Next form the M matrix for calculating the scattering matrix
		! This diagonal matrix takes into account the orientation of incident beam compared to surface normal
		! See equation (5) in Palatinus
		! Becomes identity matrix if g vectors are parallel to surface of crystal(dot product goes to 0)
		! Scattering matrix then becomes same as Kirkland
		
		
		! Calculates each element on diagonal and inputs to matrix
		DO ind = 1, nBeams
			RElement = 1/((1 + (RgDotNorm(ind)/RKn))**0.5)
			RMmatrix(ind, ind) = RElement
		END DO
		
		IF(l_alert(IErr,"angle_correction","M_matrix_formation, form matrix")) CALL abort
		
		
		! Determining weighting coefficients for bloch waves
		! Takes eigenvectors matrix
		! Inverts it
		! Weighting coefficients can then be found from
		! Product of matrix and initial conditions on wavefunction
		
		
		! Boundary conditions on surface of sample
		! These particular conditions are only valid for singular incident plane wave
		CInitial_wavefunction = CZERO ! All diffracted beams are zero
		CInitial_wavefunction(1) = CONE ! 000 beam has unit amplitude
		
		! Invert matrix of eignvectors
		CInverted_vectors = CALL INVERT(nBeams, CEigenvector_matrix)
		IF(l_alert(IErr,"angle_correction","weighting_coefficients, inverted_vectors")) CALL abort
		
		! Finds weighting coefficients
		CWeighting_Coefficients = MATMUL(CInverted_vectors, CInitial_wavefunction)
		IF(l_alert(IErr,"angle_correction","weighting_coefficients, alpha")) CALL abort
		
		
		
		! Compute scattering matrix from the eigenvalue, eignvector and M matrix
		
		! Eigenvector and M matrix need to be inverted for scattering matrix
		
		! First eigenvector matrix
		CInverted_eigen = CALL INVERT(nBeams, CEigenvector_matrix)
		IF(l_alert(IErr,"angle_correction","scattering_matrix, inverted_eigenvectors")) CALL abort
		
		! Now M matrix
		RInverted_M = CALL INVERT(nBeams, RMmatrix)
		IF(l_alert(IErr,"angle_correction","scattering_matrix, inverted_M")) CALL abort
		
		!Calculating scattering matrix
		CScatter_matrix = MATMUL(MATMUL(MATMUL(RMmatrix, CEigenvector_matrix), &
						  REigenvalue_matrix), MATMUL(CInverted_eigen, RInverted_M))
		IF(l_alert(IErr,"angle_correction","scattering_matrix, scattering")) CALL abort
		
		!!!! Alternative Scattering matrix that operates on weighting coefficients instead of boundary conditions
		
		! Forming alternative scatter matrix
		CAlt_scatter_matrix = MATMUL(RMmatrix, MATMUL(CEigenvector_matrix, &
										MATMUL(REigenvalue_matrix, &
										MATMUL(MATMUL(CInverted_eigen, RInverted_M), &
										CEigenvector_matrix)))
		IF(l_alert(IErr,"angle_correction","scattering_matrix, alt_scattering")) CALL abort
		
		
		! Now find the diffracted intensities from scattering matrices
		
		! 2 Methods can be used here
		! Calculates wavefunction column vector from
		! Product of scattering matrix and wavefunction on surface(boundary condition)
		! Take column vector of wavefunction
		! Inverse fourier transform of wavefunction to get in terms of all space
		! Get intensity from modulus squared of each element
		! Or use the weighting coefficients for bloch waves
		
		! Boundary conditions on surface of sample
		! These particular conditions are only valid for single incident plane wave
		CInitial_wavefunction = CZERO ! All diffracted beams are zero
		CInitial_wavefunction(1) = CONE ! 000 beam has amplitude of unity
		
		! First find wavefunction at any thickness in specimen, psi(z).
		CFinal_wavefunction_z = MATMUL(CScatter_matrix, CInitial_wavefunction)
		IF(l_alert(IErr,"angle_correction","diffracted_intensities, final_wavefunction_z")) CALL abort
		
		! Now need to find total wavefunction as a function of all space
		! Do 2D inverse fourier transform on CFinal_wavefunction_z
		CFinal_wavefunction = CALL IFFT2(CFinal_wavefunction_z)
		IF(l_alert(IErr,"angle_correction","diffracted_intensities, final_wavefunction")) CALL abort
		
		! Modulus squared of each element of final wavefunction for intensity
		! Could maybe do this without a loop(?)
		DO ind = 1, nBeams
			RIntensity_vector(ind) = CFinal_wavefunction(ind) * CONJG(CFinal_wavefunction(ind))
		END DO
		
		IF(l_alert(IErr,"angle_correction","diffracted_intensities, RIntensity_vector")) CALL abort
		
		!!!! Alternative way to find wavefunction as function of depth
		
		! Wavefunction at any thickness, psi(z)
		CFinal_wavefunction_z = MATMUL(CAlt_scatter_matrix, CWeighting_Coefficients)
		
		! Total wavefunction is then found from inverse fourier transform
		CFinal_wavefunction = CALL IFFT2(CFinal_wavefunction_z)
		IF(l_alert(IErr,"angle_correction","diffracted_intensities, alt_final_wavefunction")) CALL abort
		
		! Modulus squared of each element of final wavefunction for intensity
		! Again, could maybe do this without a loop(?)
		DO ind = 1, nBeams
			RIntensity_vector(ind) = CFinal_wavefunction(ind) * CONJG(CFinal_wavefunction(ind))
		END DO
		
		RETURN RIntensity_vector

	END SUBROUTINE angle_correction
	
END MODULE bloch_mod
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
