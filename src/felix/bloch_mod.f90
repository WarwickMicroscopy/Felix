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
                    IFirstPixelToCalculate,IErr)

    ! accesses procedure ZGEMM
    USE MyNumbers
    USE MyMPI
    USE message_mod
    
    USE test_koch_mod
  
    ! globals - output
    USE RPara, ONLY : RIndividualReflections ! RIndividualReflections( LACBED_ID, thickness_ID, local_pixel_ID )
    USE CPara, ONLY : CAmplitudeandPhase
    USE IPara, ONLY : IPixelComputed

    ! globals - input  
    USE CPara, ONLY : CUgMat
    USE RPara, ONLY : RDeltaK,RDeltaThickness,RInitialThickness,RNormDirM,RgDotNorm,RgPool,&
                      RgPoolMag,Rhkl
    USE IPara, ONLY : IHKLSelectFLAG,IHolzFLAG,IImageFLAG,IMinStrongBeams,IMinWeakBeams,&
                      INoOfLacbedPatterns,IPixelCount,IThicknessCount,nReflections,&
                      IOutputReflections
    USE BlochPara, ONLY : RBigK            
    
    IMPLICIT NONE
    
    INTEGER(IKIND),INTENT(IN) :: IYPixelIndex,IXPixelIndex,IPixelNumber,&
          IFirstPixelToCalculate
    INTEGER(IKIND),INTENT(OUT) :: IErr
    
    COMPLEX(CKIND),ALLOCATABLE :: CBeamProjectionMatrix(:,:),&
          CDummyBeamMatrix(:,:),CUgSgMatrix(:,:),CEigenVectors(:,:),CEigenValues(:),&
          CInvertedEigenVectors(:,:),CAlphaWeightingCoefficients(:),&
          CEigenValueDependentTerms(:,:)
    COMPLEX(CKIND) :: CFullWaveFunctions(nReflections)
    REAL(RKIND) :: RFullWaveIntensity(nReflections),RDevPara(nReflections),&
          RTiltedK(ITHREE)
    INTEGER(IKIND) :: IStrongBeamList(nReflections),IWeakBeamList(nReflections),&
          nBeams,nWeakBeams
    INTEGER(IKIND) :: ind,knd,pnd,IThickness,IThicknessIndex,ILowerLimit,&
          IUpperLimit       
    REAL(RKIND) :: RThickness,RKn
    COMPLEX(CKIND) sumC,sumD
    COMPLEX(CKIND), DIMENSION(:,:), ALLOCATABLE :: CGeneralSolutionMatrix, &
         CGeneralEigenSpectrumEigenVectors,CBeamTranspose,CUgMatPartial
    COMPLEX(CKIND),DIMENSION(:),ALLOCATABLE :: CGeneralEigenValues
    CHARACTER*40 surname
    CHARACTER*100 SindString,SjndString,SPixelCount,SnBeams,SWeakBeamIndex,SPrintString

    ! variables used for koch spence method development
    LOGICAL,PARAMETER :: TestKochMethod = .FALSE.
    COMPLEX(CKIND),ALLOCATABLE :: CDiagonalSgMatrix(:,:), COffDiagonalSgMatrix(:,:)
    COMPLEX(CKIND) :: CScatteringElement
    INTEGER(IKIND) :: ScatterMatrixRow
      
    ! we are inside the mask
    IPixelComputed= IPixelComputed + 1

    ! TiltedK is the vector of the incoming tilted beam
    ! in units of (1/A), in the microscope ref frame(NB exp(i*k.r), physics convention)
    ! x-position in k-space
    RTiltedK(1)= (REAL(IYPixelIndex,RKIND)-REAL(IPixelCount,RKIND)-0.5_RKIND)*RDeltaK
    ! y-position in k-space
    RTiltedK(2)= (REAL(IXPixelIndex,RKIND)-REAL(IPixelCount,RKIND)-0.5_RKIND)*RDeltaK 
    RTiltedK(3)= SQRT(RBigK**2 - RTiltedK(1)**2 - RTiltedK(2)**2) 
    RKn = DOT_PRODUCT(RTiltedK,RNormDirM)
    
    ! Compute the deviation parameter for reflection pool
    ! NB RDevPara is in units of (1/A)
    ! in the microscope ref frame(NB exp(i*s.r), physics convention)
    DO knd=1,nReflections
      ! Sg parallel to z: Sg=-[k'z+gz-sqrt( (k'z+gz)^2-2k'.g-g^2)]
      RDevPara(knd)= -RTiltedK(3)-RgPool(knd,3)+&
	    SQRT( (RTiltedK(3)+RgPool(knd,3))**2 - &
            2*DOT_PRODUCT(RgPool(knd,:),RTiltedK(:)) - RgPoolMag(knd)**2 )
      !??  remove below commented out
      !Keith's old version, Sg parallel to k'
      !RDevPara(knd)= -( RBigK + DOT_PRODUCT(RgPool(knd,:),RTiltedK(:)) /RBigK) + &
      !  SQRT( ( RBigK**2 + DOT_PRODUCT(RgPool(knd,:),RTiltedK(:)) )**2 /RBigK**2 - &
      !  (RgPoolMag(knd)**2 + TWO*DOT_PRODUCT(RgPool(knd,:),RTiltedK(:))) )
      IF(knd.EQ.2.AND.IYPixelIndex.EQ.10.AND.IXPixelIndex.EQ.10) THEN
        CALL message(LM,dbg7,"RBigK ",RBigK)
        CALL message(LM,dbg7,"Rhkl(knd) ",Rhkl(knd:knd,:))
        CALL message(LM,dbg7,"RgPool(knd) ",RgPool(knd:knd,:))
        CALL message(LM,dbg7,"RTiltedK ",RTiltedK)
        CALL message(LM,dbg7,"RDevPara ",RDevPara(knd))
      END IF
    END DO

    ! select only those beams where the Ewald sphere is close to the
    ! reciprocal lattice, i.e. within RBSMaxDeviationPara
    CALL StrongAndWeakBeamsDetermination(nReflections,IMinWeakBeams,&
                    IMinStrongBeams,RDevPara,CUgMat,&
                    IStrongBeamList,IWeakBeamList,nBeams,nWeakBeams,IErr)
    IF(l_alert(IErr,"BlochCoefficientCalculation",&
          "StrongAndWeakBeamsDetermination()")) RETURN
    CALL message(LXL,dbg7,"strong beams",nBeams)
    CALL message(LXL,dbg7,"weak beams",nWeakBeams)
    CALL message(LXL,dbg7,"nReflections",nReflections)

    !--------------------------------------------------------------------
    ! ALLOCATE memory for eigen problem
    !--------------------------------------------------------------------

    ! now nBeams determined, allocate complex arrays
    ALLOCATE( CBeamProjectionMatrix(nBeams,nReflections), STAT=IErr )
    IF(l_alert(IErr,"BlochCoefficientCalculation","allocate CBeamProjectionMatrix")) RETURN
    ALLOCATE( CDummyBeamMatrix(nBeams,nReflections), STAT=IErr )
    IF(l_alert(IErr,"BlochCoefficientCalculation","allocate CBeamProjectionMatrix")) RETURN
    ALLOCATE( CUgSgMatrix(nBeams,nBeams), STAT=IErr )
    IF(l_alert(IErr,"BlochCoefficientCalculation","allocate CBeamProjectionMatrix")) RETURN
    ALLOCATE( CEigenValues(nBeams), STAT=IErr )
    IF(l_alert(IErr,"BlochCoefficientCalculation","allocate CBeamProjectionMatrix")) RETURN
    ALLOCATE( CEigenVectors(nBeams,nBeams), STAT=IErr )
    IF(l_alert(IErr,"BlochCoefficientCalculation","allocate CBeamProjectionMatrix")) RETURN
    ALLOCATE( CInvertedEigenVectors(nBeams,nBeams), STAT=IErr )
    IF(l_alert(IErr,"BlochCoefficientCalculation","allocate CBeamProjectionMatrix")) RETURN
    ALLOCATE( CBeamTranspose(nReflections,nBeams), STAT=IErr )
    IF(l_alert(IErr,"BlochCoefficientCalculation","allocate CBeamProjectionMatrix")) RETURN
    ALLOCATE( CUgMatPartial(nReflections,nBeams), STAT=IErr )
    IF(l_alert(IErr,"BlochCoefficientCalculation","allocate CBeamProjectionMatrix")) RETURN
    ALLOCATE( CAlphaWeightingCoefficients(nBeams), STAT=IErr )
    IF(l_alert(IErr,"BlochCoefficientCalculation","allocate CBeamProjectionMatrix")) RETURN
    ALLOCATE( CEigenValueDependentTerms(nBeams,nBeams), STAT=IErr )
    IF(l_alert(IErr,"BlochCoefficientCalculation","allocate CBeamProjectionMatrix")) RETURN

    ! allocations used for koch spence method development
    IF(TestKochMethod) THEN
      ALLOCATE( CDiagonalSgMatrix(nBeams,nBeams), STAT=IErr )
      IF(l_alert(IErr,"BlochCoefficientCalculation","allocate CDiagonalSgMatrix")) RETURN
      ALLOCATE( COffDiagonalSgMatrix(nBeams,nBeams), STAT=IErr )
      IF(l_alert(IErr,"BlochCoefficientCalculation","allocate COffDiagonalSgMatrix")) RETURN
    END IF

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
    CALL ZGEMM('N','N',nReflections,nBeams,nReflections,CONE,CUgMat, &
              nReflections,CBeamTranspose,nReflections,CZERO,CUgMatPartial,nReflections)
    ! CUgSgMatrix = CBeamProjectionMatrix * CUgMatPartial
    CALL ZGEMM('N','N',nBeams,nBeams,nReflections,CONE,CBeamProjectionMatrix, &
              nBeams,CUgMatPartial,nReflections,CZERO,CUgSgMatrix,nBeams)

    !--------------------------------------------------------------------
    ! higher order Laue zones
    !--------------------------------------------------------------------

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

          !??  remove commented Keith's old version
          !REAL(CUgMat(IStrongBeamList(knd),IWeakBeamList(ind))) * &
          !REAL(CUgMat(IWeakBeamList(ind),1)) / &
          !(4*RBigK*RBigK*RDevPara(IWeakBeamList(ind)))
          sumD = sumD + &
		      ! Zuo&Weickenmeier Ultramicroscopy 57 (1995) 375-383 eq.5
          CUgMat(IStrongBeamList(knd),IWeakBeamList(ind))*&
          CUgMat(IWeakBeamList(ind),IStrongBeamList(knd))/&
          (TWO*RBigK*RDevPara(IWeakBeamList(ind)))
          !?? remove commented Keith's old version
          !REAL(CUgMat(IStrongBeamList(knd),IWeakBeamList(ind))) * &
          !REAL(CUgMat(IWeakBeamList(ind),IStrongBeamList(knd))) / &
          !(4*RBigK*RBigK*RDevPara(IWeakBeamList(ind)))
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
    ! diagonalize the UgMatEffective
    !--------------------------------------------------------------------

    ! If koch method - Split CUgSgMatrix into diagonal and off diagonal to speed convergence
    IF(TestKochMethod) THEN
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

    IF (IHolzFLAG.EQ.1) THEN ! higher order laue zone included so adjust Eigen values/vectors
      CEigenValues = CEigenValues * RKn/RBigK
      DO knd = 1,nBeams
        CEigenVectors(knd,:) = CEigenVectors(knd,:) / &
              SQRT(1+RgDotNorm(IStrongBeamList(knd))/RKn)
      END DO
    END IF

    !--------------------------------------------------------------------
    ! fill RIndividualReflections( LACBED_ID , thickness_ID, local_pixel_ID ) 
    !--------------------------------------------------------------------
   
    ! Calculate intensities for different specimen thicknesses
    !?? ADD VARIABLE PATH LENGTH HERE !?? what does this comment mean?
    DO IThicknessIndex=1,IThicknessCount,1

      RThickness = RInitialThickness + REAL((IThicknessIndex-1),RKIND)*RDeltaThickness
      IThickness = NINT(RThickness,IKIND)

      ! optional - for koch development to speed convergence
      IF(TestKochMethod) RThickness = RThickness / 1000

      CALL CreateWaveFunctions(RThickness,RFullWaveIntensity,CFullWaveFunctions,&
                    nReflections,nBeams,IStrongBeamList,CEigenVectors,CEigenValues,IErr)
      IF(l_alert(IErr,"BlochCoefficientCalculation","CreateWaveFunctions")) RETURN

      !--------------------------------------------------------------------
      ! Optional - test koch spence prototype method
      !--------------------------------------------------------------------

      IF(TestKochMethod) THEN

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
!                        nReflections,nBeams,IStrongBeamList,CEigenVectors,CEigenValues,IErr)
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
      IF(IHKLSelectFLAG.EQ.0) THEN ! we are not using hkl list from felix.hkl
        IF(IImageFLAG.LE.2) THEN ! output is 0=montage, 1=individual images
          RIndividualReflections(1:INoOfLacbedPatterns,IThicknessIndex,&
                (IPixelNumber-IFirstPixelToCalculate)+1) = &
                RFullWaveIntensity(1:INoOfLacbedPatterns)
        ELSE ! output is 2=amplitude+phase images
          CAmplitudeandPhase(1:INoOfLacbedPatterns,IThicknessIndex,&
                (IPixelNumber-IFirstPixelToCalculate)+1) = &
                CFullWavefunctions(1:INoOfLacbedPatterns)
        END IF
      ELSE ! we are using hkl list from felix.hkl
        IF(IImageFLAG.LE.2) THEN
          DO pnd = 1,INoOfLacbedPatterns
            RIndividualReflections(pnd,IThicknessIndex,&
                  (IPixelNumber-IFirstPixelToCalculate)+1) = &
                  RFullWaveIntensity(IOutputReflections(pnd))
          END DO
        ELSE
          DO pnd = 1,INoOfLacbedPatterns
            CAmplitudeandPhase(pnd,IThicknessIndex,(IPixelNumber-IFirstPixelToCalculate)+1)=&
                 CFullWavefunctions(IOutputReflections(pnd))
          END DO
        END IF
      END IF
    END DO

    ! deallocations used for koch spence method development
    IF(TestKochMethod) THEN
      DEALLOCATE( CDiagonalSgMatrix, COffDiagonalSgMatrix, STAT=IErr )
      IF(l_alert(IErr,"BlochCoefficientCalculation","deallocating arrays")) RETURN
    END IF
    
    ! DEALLOCATE eigen problem memory
    DEALLOCATE(CUgSgMatrix,CBeamTranspose, CUgMatPartial, &
         CInvertedEigenVectors, CAlphaWeightingCoefficients, &
         CEigenValues,CEigenVectors,CEigenValueDependentTerms, &
         CBeamProjectionMatrix, CDummyBeamMatrix,STAT=IErr)
    IF(l_alert(IErr,"BlochCoefficientCalculation","deallocating arrays")) RETURN
    
  END SUBROUTINE BlochCoefficientCalculation

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Calculates diffracted intensity for a specific thickness
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE CreateWaveFunctions(RThickness,RFullWaveIntensity,CFullWaveFunctions,&
                    nReflections,nBeams,IStrongBeamList,CEigenVectors,CEigenValues,IErr)

    USE MyNumbers
    USE MyMPI
    USE message_mod

    IMPLICIT NONE
    
    REAL(RKIND),INTENT(IN) :: RThickness
    REAL(RKIND),INTENT(OUT) :: RFullWaveIntensity(nReflections)
    COMPLEX(CKIND),INTENT(OUT) :: CFullWaveFunctions(nReflections) 
    INTEGER(IKIND),INTENT(IN) :: nReflections,nBeams,IStrongBeamList(nReflections)
    COMPLEX(CKIND),INTENT(IN) :: CEigenVectors(nBeams,nBeams),CEigenValues(nBeams)
    INTEGER(IKIND),INTENT(OUT) :: IErr 
    REAL(RKIND) :: RWaveIntensity(nBeams)
    COMPLEX(CKIND) :: CInvertedEigenVectors(nBeams,nBeams),CPsi0(nBeams),&
          CWaveFunctions(nBeams),CEigenValueDependentTerms(nBeams,nBeams),&
          CAlphaWeightingCoefficients(nBeams)
    INTEGER(IKIND) :: ind,jnd,knd,hnd,ifullind,iuniind,gnd,ichnk
    COMPLEX(CKIND) :: CDummyEigenVectors(nBeams,nBeams)
    
    ! The top surface boundary conditions
    CPsi0 = CZERO ! all diffracted beams are zero
    CPsi0(1) = CONE ! the 000 beam has unit amplitude
    
    ! Invert the EigenVector matrix
    CDummyEigenVectors = CEigenVectors
    CALL INVERT(nBeams,CDummyEigenVectors(:,:),CInvertedEigenVectors,IErr)

    ! put in the thickness
    ! From EQ 6.32 in Kirkland Advance Computing in EM
    CAlphaWeightingCoefficients = MATMUL(CInvertedEigenVectors(1:nBeams,1:nBeams),CPsi0) 
    CEigenValueDependentTerms= CZERO
    DO hnd=1,nBeams     ! This is a diagonal matrix
      CEigenValueDependentTerms(hnd,hnd) = &
            EXP(CIMAGONE*CMPLX(RThickness,ZERO,CKIND)*CEigenValues(hnd)) 
    ENDDO
    ! The diffracted intensity for each beam
    ! EQ 6.35 in Kirkland Advance Computing in EM
    ! C-1*C*alpha 
    CWaveFunctions(:) = MATMUL( &
          MATMUL(CEigenVectors(1:nBeams,1:nBeams),CEigenValueDependentTerms), & 
          CAlphaWeightingCoefficients(:) )

    !?? possible small time saving here by only calculating the (tens of) output
    !?? reflections rather than all strong beams (hundreds)
    DO hnd=1,nBeams
       RWaveIntensity(hnd)=CONJG(CWaveFunctions(hnd)) * CWaveFunctions(hnd)
    ENDDO  

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
  !! pertubation strengths and iterates over the number of weak and strong until
  !! there are enough.
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE StrongAndWeakBeamsDetermination(nReflections,IMinWeakBeams,&
                    IMinStrongBeams,RDevPara,CUgMat,&
                    IStrongBeamList,IWeakBeamList,nBeams,nWeakBeams,IErr)
    
    ! select only those beams where the Ewald sphere is close to the
    ! reciprocal lattice, i.e. within RBSMaxDeviationPara

    USE MyNumbers
    USE MyMPI
    USE message_mod 

    INTEGER(IKIND),INTENT(IN) :: nReflections
    REAL(RKIND),DIMENSION(nReflections),INTENT(IN) :: RDevPara
    COMPLEX(CKIND),DIMENSION(nReflections,nReflections),INTENT(IN) :: CUgMat
    INTEGER(IKIND),INTENT(IN) :: IMinWeakBeams, IMinStrongBeams
    INTEGER(IKIND),DIMENSION(nReflections),INTENT(OUT) :: IStrongBeamList,IWeakBeamList
    INTEGER(IKIND),INTENT(OUT) :: nBeams,nWeakBeams,IErr
    INTEGER(IKIND) :: ind,jnd
    INTEGER(IKIND),DIMENSION(:) :: IStrong(nReflections),IWeak(nReflections)
    REAL(RKIND) :: RMaxSg,RMinPertStrong,RMinPertWeak
    REAL(RKIND),DIMENSION(:) :: RPertStrength0(nReflections)

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
    DO jnd=1,nReflections
      IF (IStrong(jnd).EQ.1) THEN
	      IStrongBeamList(ind)=jnd
        ind=ind+1
	    END IF
    END DO

    ! this is used to give the dimension of the Bloch wave problem
    nBeams=ind-1  

    CALL message(LXL,dbg7,"Strong Beam List",IStrongBeamList)
    CALL message(LXL,dbg7,"Sg limit for strong beams = ",RMaxSg)
    CALL message(LXL,dbg7,"Smallest strong perturbation strength = ",RMinPertStrong)
    IF(SUM(IStrong)+IMinWeakBeams.GT.nReflections) IErr = 1
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
    DO jnd=1,nReflections
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
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
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
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
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
