!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all rights reserved
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
!  This file is part of felixsim.
!
!  felixsim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  felixsim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with felixsim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE BlochCoefficientCalculation(IYPixelIndex,IXPixelIndex,IPixelNumber,IFirstPixelToCalculate,IErr)
  
  USE WriteToScreen
  USE MyNumbers
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IYPixelIndex,IXPixelIndex,ind,knd,IPixelNumber,pnd,&
       ierr,IThickness,IThicknessIndex,ILowerLimit,IUpperLimit,IFirstPixelToCalculate       
  REAL(RKIND) :: RThickness,RKn
  COMPLEX(CKIND) sumC,sumD
  COMPLEX(CKIND), DIMENSION(:,:), ALLOCATABLE :: CGeneralSolutionMatrix, &
       CGeneralEigenVectors,CBeamTranspose,CUgMatPartial
  COMPLEX(CKIND),DIMENSION(:),ALLOCATABLE :: CGeneralEigenValues
  CHARACTER*40 surname
  CHARACTER*200 SindString,SjndString,SPixelCount,SnBeams,SWeakBeamIndex,SPrintString
   
  IF (my_rank.EQ.0) THEN
    DO WHILE (IMessageCounter .LT.1)
      CALL Message("BlochCoefficientCalculation",IMust,IErr)
      CALL Message("BlochCoefficientCalculation",IMust+IDebug,IErr, & 
              MessageString = "is looping, and calling subroutines itself, They are:")
      IMessageCounter = IMessageCounter +1
    END DO
  END IF
    
  ! we are inside the mask
  IPixelComputed= IPixelComputed + 1

  !!$   Displays Pixel currently working on
  WRITE(SindString,'(I6.1)') IYPixelIndex
  WRITE(SjndString,'(I6.1)') IXPixelIndex
  WRITE(SPixelCount,'(I6.1)') 2*IPixelCount
  CALL Message("BlochCoefficientCalculation",IAllInfo,IErr, &
       MessageString="working on pixel("//TRIM(ADJUSTL(SindString))//",&
       &"//TRIM(ADJUSTL(SjndString))//") of ("//TRIM(ADJUSTL(SPixelCount))//",&
       &"//TRIM(ADJUSTL(SPixelCount))//") in total")

  !--------------------------------------------------------------------
  ! TiltedK is the vector of the incoming tilted beam
  ! in units of (1/A), in the microscope ref frame(NB exp(i*k.r), physics convention)
  RTiltedK(1)= (REAL(IYPixelIndex,RKIND)-REAL(IPixelCount,RKIND)-0.5_RKIND)*RDeltaK ! x-position in k-space
  RTiltedK(2)= (REAL(IXPixelIndex,RKIND)-REAL(IPixelCount,RKIND)-0.5_RKIND)*RDeltaK ! y-position in k-space
  RTiltedK(3)= SQRT(RBigK**2 - RTiltedK(1)**2 - RTiltedK(2)**2) 
  RKn = DOT_PRODUCT(RTiltedK,RNormDirM)
  
  !Compute the deviation parameter for reflection pool
  !NB RDevPara is in units of (1/A), in the microscope ref frame(NB exp(i*s.r), physics convention)
  DO knd=1,nReflections
    !Sg parallel to z: Sg=-[k'z+gz-sqrt( (k'z+gz)^2-2k'.g-g^2)]
    RDevPara(knd)= -RTiltedK(3)-RgPool(knd,3)+&
	SQRT( (RTiltedK(3)+RgPool(knd,3))**2-2*DOT_PRODUCT(RgPool(knd,:),RTiltedK(:))-RgPoolMag(knd)**2)
	!Keith's old version, Sg parallel to k'
    !RDevPara(knd)= -( RBigK + DOT_PRODUCT(RgPool(knd,:),RTiltedK(:)) /RBigK) + &
    !  SQRT( ( RBigK**2 + DOT_PRODUCT(RgPool(knd,:),RTiltedK(:)) )**2 /RBigK**2 - &
    !  (RgPoolMag(knd)**2 + TWO*DOT_PRODUCT(RgPool(knd,:),RTiltedK(:))) )
    IF(IWriteFLAG.EQ.6.AND.knd.EQ.2.AND.IYPixelIndex.EQ.10.AND.IXPixelIndex.EQ.10) THEN
      PRINT*,"RBigK",RBigK
      PRINT*,"Rhkl(knd)",Rhkl(knd,:)
      PRINT*,"RgPool(knd)",RgPool(knd,:)
      PRINT*,"RTiltedK",RTiltedK
      PRINT*,"RDevPara",RDevPara(knd)
    END IF
  END DO
  
  ! select only those beams where the Ewald sphere is close to the
  ! reciprocal lattice, i.e. within RBSMaxDeviationPara
  CALL StrongAndWeakBeamsDetermination(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(",my_rank,") error in Determination of Strong and Weak beams"
     RETURN
  END IF
  IF(IWriteFLAG.EQ.3.AND.IYPixelIndex.EQ.10.AND.IXPixelIndex.EQ.10) THEN
    PRINT*, IStrongBeamIndex,"strong beams"
    PRINT*, IWeakBeamIndex,"weak beams"
  END IF  ! select the highest reflection that corresponds to a strong beam
  nBeams= IStrongBeamIndex

  !--------------------------------------------------------------------
  ! ALLOCATE memory for eigen problem
  !--------------------------------------------------------------------
  ALLOCATE(CBeamProjectionMatrix(nBeams,nReflections),STAT=IErr)
  ALLOCATE(CDummyBeamMatrix(nBeams,nReflections),STAT=IErr)
  ALLOCATE(CUgSgMatrix(nBeams,nBeams),STAT=IErr)
  ALLOCATE(CEigenVectors(nBeams,nBeams),STAT=IErr)
  ALLOCATE(CEigenValues(nBeams),STAT=IErr)
  ALLOCATE(CInvertedEigenVectors(nBeams,nBeams),STAT=IErr)
  ALLOCATE(CBeamTranspose(nReflections,nBeams),STAT=IErr)
  ALLOCATE(CUgMatPartial(nReflections,nBeams),STAT=IErr)
  ALLOCATE(CAlphaWeightingCoefficients(nBeams),STAT=IErr)
  ALLOCATE(CEigenValueDependentTerms(nBeams,nBeams),STAT=IErr)
  ALLOCATE(CWaveFunctions(nBeams),STAT=IErr)
  ALLOCATE(RWaveIntensity(nBeams),STAT=IErr)
  ALLOCATE(CPsi0(nBeams),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(",my_rank,")error in allocations"
     RETURN
  END IF
  
  !--------------------------------------------------------------------
  ! construct the effective UgMat (for strong beams only at the moment)
  WRITE(SnBeams,"(I6.1)") nBeams
  WRITE(SWeakBeamIndex,"(I6.1)") IWeakBeamIndex
  CALL Message("BlochCoefficientCalculation",IAllInfo,IErr, &
       MessageString="using n(Strong) Beams = "//ADJUSTL(TRIM(SnBeams))// &
       "with nWeakBeams = "// ADJUSTL(TRIM(SWeakBeamIndex)))

  ! compute the effective Ug matrix by selecting only those beams
  ! for which IStrongBeamList has an entry
  CBeamProjectionMatrix= CZERO
  DO knd=1,nBeams
     CBeamProjectionMatrix(knd,IStrongBeamList(knd))=CONE
  ENDDO
  CUgSgMatrix = CZERO
  CBeamTranspose=TRANSPOSE(CBeamProjectionMatrix)

  CALL ZGEMM('N','N',nReflections,nBeams,nReflections,CONE,CUgMat, &
       nReflections,CBeamTranspose,nReflections,CZERO,CUgMatPartial,nReflections)
  CALL ZGEMM('N','N',nBeams,nBeams,nReflections,CONE,CBeamProjectionMatrix, &
       nBeams,CUgMatPartial,nReflections,CZERO,CUgSgMatrix,nBeams)

  IF (IZolzFLAG.EQ.0) THEN
    DO ind=1,nBeams
      CUgSgMatrix(ind,ind) = CUgSgMatrix(ind,ind) + TWO*RBigK*RDevPara(IStrongBeamList(ind))
    ENDDO
    DO knd =1,nBeams ! Columns
      DO ind = 1,nBeams ! Rows
        CUgSgMatrix(knd,ind) = CUgSgMatrix(knd,ind) / &
         (SQRT(1+RgDotNorm(IStrongBeamList(knd))/RKn)*SQRT(1+RgDotNorm(IStrongBeamList(ind))/RKn))
      END DO
    END DO
    CUgSgMatrix = (TWOPI**2)*CUgSgMatrix/(TWO*RBigK)
  ELSE
    ! replace the diagonal parts with strong beam deviation parameters
    DO ind=1,nBeams
      CUgSgMatrix(ind,ind) = TWO*RBigK*RDevPara(IStrongBeamList(ind))/(TWOPI*TWOPI)
    ENDDO
    ! add the weak beams perturbatively for the 1st column (sumC) and
    ! the diagonal elements (sumD)
    DO knd=2,nBeams
      sumC=CZERO
      sumD=CZERO
      DO ind=1,IWeakBeamIndex
       sumC=sumC + &
		!Zuo&Weickenmeier Ultramicroscopy 57 (1995) 375-383 eq.4
        CUgMat(IStrongBeamList(knd),IWeakBeamList(ind))*&
        CUgMat(IWeakBeamList(ind),1)/(TWO*RBigK*RDevPara(IWeakBeamList(ind)))
!Keith's old version
!          REAL(CUgMat(IStrongBeamList(knd),IWeakBeamList(ind))) * &
!          REAL(CUgMat(IWeakBeamList(ind),1)) / &
!         (4*RBigK*RBigK*RDevPara(IWeakBeamList(ind)))
        sumD = sumD + &
		!Zuo&Weickenmeier Ultramicroscopy 57 (1995) 375-383 eq.5
        CUgMat(IStrongBeamList(knd),IWeakBeamList(ind))*&
        CUgMat(IWeakBeamList(ind),IStrongBeamList(knd))/&
        (TWO*RBigK*RDevPara(IWeakBeamList(ind)))
!Keith's old version
!          REAL(CUgMat(IStrongBeamList(knd),IWeakBeamList(ind))) * &
!          REAL(CUgMat(IWeakBeamList(ind),IStrongBeamList(knd))) / &
!         (4*RBigK*RBigK*RDevPara(IWeakBeamList(ind)))
      ENDDO
	  !Replace the Ug's
	  WHERE (CUgSgMatrix.EQ.CUgSgMatrix(knd,1))
        CUgSgMatrix= CUgSgMatrix(knd,1) - sumC
	  END WHERE
	  !Replace the Sg's
      CUgSgMatrix(knd,knd)= CUgSgMatrix(knd,knd) - TWO*RBigK*sumD/(TWOPI*TWOPI)
    ENDDO
	!Divide by 2K so off-diagonal elementa are Ug/2K, diagonal elements are Sg
	!DON'T KNOW WHERE THE 4pi^2 COMES FROM!! 
    CUgSgMatrix = TWOPI*TWOPI*CUgSgMatrix/(TWO*RBigK)
  END IF

  IF(IWriteFLAG.EQ.3.AND.IYPixelIndex.EQ.10.AND.IXPixelIndex.EQ.10) THEN
   PRINT*,"Ug/2K + {Sg} matrix (nm^-2)"
	DO ind =1,6
     WRITE(SPrintString,FMT='(3(1X,I3),A1,8(1X,F7.3,F7.3))') NINT(Rhkl(ind,:)),":",100*CUgSgMatrix(ind,1:6)
     PRINT*,TRIM(SPrintString)
    END DO
  END IF	
  
  !--------------------------------------------------------------------
  ! diagonalize the UgMatEffective
  IF (IZolzFLAG.EQ.0) THEN
     CALL EigenSpectrum(nBeams,CUgSgMatrix,CEigenValues(:), CEigenVectors(:,:),IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in EigenSpectrum()"
        RETURN
     END IF
     CEigenValues = CEigenValues * RKn/RBigK
     DO knd = 1,nBeams
        CEigenVectors(knd,:) = CEigenVectors(knd,:) / &
             SQRT(1+RgDotNorm(IStrongBeamList(knd))/RKn)
     END DO
  ELSE
     CALL EigenSpectrum(nBeams,CUgSgMatrix,CEigenValues(:),CEigenVectors(:,:),IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in EigenSpectrum()"
        RETURN
     END IF
  END IF
 
  DO IThicknessIndex=1,IThicknessCount,1
     RThickness = RInitialThickness + REAL((IThicknessIndex-1),RKIND)*RDeltaThickness 
     IThickness = NINT(RInitialThickness + REAL((IThicknessIndex-1),RKIND)*RDeltaThickness,IKIND) 
     CALL CreateWaveFunctions(RThickness,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in CreateWavefunction()"
        RETURN
     END IF
     !Collection Wave Intensities from all thickness for later writing
     IF(IHKLSelectFLAG.EQ.0) THEN
        IF(IImageFLAG.LE.2) THEN
           RIndividualReflections(1:INoOfLacbedPatterns,IThicknessIndex,(IPixelNumber-IFirstPixelToCalculate)+1) = &
                RFullWaveIntensity(1:INoOfLacbedPatterns)
        ELSE
           CAmplitudeandPhase(1:INoOfLacbedPatterns,IThicknessIndex,(IPixelNumber-IFirstPixelToCalculate)+1) = &
                CFullWavefunctions(1:INoOfLacbedPatterns)
        END IF
     ELSE
        IF(IImageFLAG.LE.2) THEN
           DO pnd = 1,INoOfLacbedPatterns
              RIndividualReflections(pnd,IThicknessIndex,(IPixelNumber-IFirstPixelToCalculate)+1) = &
                   RFullWaveIntensity(IOutputReflections(pnd))
           END DO
        ELSE
           DO pnd = 1,INoOfLacbedPatterns
              CAmplitudeandPhase(pnd,IThicknessIndex,(IPixelNumber-IFirstPixelToCalculate)+1) = &
                   CFullWavefunctions(IOutputReflections(pnd))
           END DO
        END IF
     END IF
  END DO
  
  !--------------------------------------------------------------------
  ! DEALLOCATE eigen problem memory
  DEALLOCATE(CUgSgMatrix,CPsi0,CBeamTranspose, CUgMatPartial, &
       CInvertedEigenVectors, CAlphaWeightingCoefficients, &
       CEigenValues,CEigenVectors,CEigenValueDependentTerms, &
       CBeamProjectionMatrix, CDummyBeamMatrix,CWavefunctions, &
       RWaveIntensity,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in Deallocation()"
     RETURN
  ENDIF
  
END SUBROUTINE BlochCoefficientCalculation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CreateWaveFunctions(RThickness,IErr)

  USE WriteToScreen
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,jnd,knd,hnd,IErr, ifullind, iuniind,gnd,ichnk
  REAL(RKIND) :: RThickness 
  COMPLEX(CKIND),DIMENSION(:,:),ALLOCATABLE :: CDummyEigenVectors

  IF (my_rank.EQ.0) THEN
     DO WHILE (IMessageCounter .LT.6)
        CALL Message("CreateWaveFunctions",IMust,IErr)
        IMessageCounter = IMessageCounter +1
     END DO
  END IF
  
  !--------------------------------------------------------------------
  ! calculate wavefunctions
  !--------------------------------------------------------------------
  CPsi0 = CZERO
  IF(nBeams .GE. 0) CPsi0(1) = CONE
  
  ALLOCATE(CDummyEigenVectors(nBeams,nBeams),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateWavefunctions(",my_rank,")error allocating CDummyEigenVectors"
     RETURN
  ENDIF
  
  ! Invert the EigenVector matrix
  CDummyEigenVectors = CEigenVectors
  CALL INVERT(nBeams,CDummyEigenVectors(:,:),CInvertedEigenVectors,IErr)

  !From EQ 6.32 in Kirkland Advance Computing in EM
  CAlphaWeightingCoefficients = MATMUL(CInvertedEigenVectors(1:nBeams,1:nBeams),CPsi0) 
  CEigenValueDependentTerms= CZERO
  DO hnd=1,nBeams     ! This is a diagonal matrix
    CEigenValueDependentTerms(hnd,hnd)=EXP(CIMAGONE*CMPLX(RThickness,ZERO,CKIND)*CEigenValues(hnd)) 
  ENDDO
  
  ! EQ 6.35 in Kirkland Advance Computing in EM
  ! C-1*C*alpha 
  CWaveFunctions(:)=MATMUL(MATMUL(CEigenVectors(1:nBeams,1:nBeams),CEigenValueDependentTerms), & 
       CAlphaWeightingCoefficients(:) )
  DO hnd=1,nBeams
     RWaveIntensity(hnd)=CONJG(CWaveFunctions(hnd)) * CWaveFunctions(hnd)
  ENDDO  
  
  !--------------------------------------------------------------------
  ! rePADDing of wave function and intensities with zero's 
  !--------------------------------------------------------------------
  CFullWaveFunctions=CZERO
  RFullWaveIntensity=ZERO
  DO knd=1,nBeams
     CFullWaveFunctions(IStrongBeamList(knd))=CWaveFunctions(knd)
     RFullWaveIntensity(IStrongBeamList(knd))=RWaveIntensity(knd)
  ENDDO
  
  DEALLOCATE(CDummyEigenVectors,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateWavefunctions(",my_rank,")error deallocating CDummyEigenVectors"
     RETURN
  ENDIF
  
END SUBROUTINE CreateWavefunctions

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE StrongAndWeakBeamsDetermination(IErr)
  
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI
  
  INTEGER(IKIND) :: ind,knd,IErr,IMinimum,IMaximum,ICheck,jnd,hnd, &
       IAdditionalBmaxStrongBeams,IAdditionalPmaxStrongBeams,&
       IBeamIterationCounter,IFound
  REAL(RKIND) :: RDummySg(nReflections),RPertStrength,RMinPertStrength,sumC
  INTEGER(IKIND), DIMENSION(:),ALLOCATABLE  :: IAdditionalBmaxStrongBeamList,IAdditionalPmaxStrongBeamList

  IF (my_rank.EQ.0) THEN
    DO WHILE (IMessageCounter .LT.4)
      CALL Message("StrongAndWeakBeamsDetermination",IMust,IErr)
      IMessageCounter = IMessageCounter +1
     END DO
  END IF

  !----------------------------------------------------------------------------
  !STRONG BEAMS
  !Use perturbation strength to define strong beams
  !PerturbationStrength Eq. 8 Zuo Ultramicroscopy 57 (1995) 375, |Ug/2KSg|
  !Use IMinStrongBeams to calculate the maximum value of Sg, RBSMaxDeviationPara
  IStrongBeamList = 0_IKIND
  RBSMaxDeviationPara = ZERO
  RMinPertStrength=0.00001!1e-04 gives about 140 strong beams In Ca3Mn2O7
  RDummySg = ABS(RDevPara)!list of ABS(Sg) for the different reflections for this pixel
!  DO ind=1,IMinStrongBeams
!    !find the smallest Sg
!    IMinimum = MINLOC(RDummySg,1)
!	!does it have the largest Sg so far?
!    IF(RDummySg(IMinimum).GT.RBSMaxDeviationPara) THEN
!      RBSMaxDeviationPara = ABS(RDevPara(IMinimum))
!	END IF
!	!does it have the smallest perturbation strength so far?
!    RPertStrength = ABS(CUgMat(IMinimum,1)/(TWO*RBigK*RDevPara(IMinimum)))
!    IF(RPertStrength.LT.RMinPertStrength) THEN
!      RMinPertStrength = RPertStrength
!    END IF
!    RDummySg(IMinimum)=-NEGHUGE!finished with this one, on to the next
!  END DO
  IF(IWriteFLAG.EQ.6.AND.my_rank.EQ.0) THEN
!    PRINT*, "Sg limit for strong beams=",RBSMaxDeviationPara
    PRINT*, "Minimum perturbation for strong beams=",RMinPertStrength
  END IF
  !Count the reflectionswith Sg <= RBSMaxDeviationPara or RPertStrength>RMinPertStrength
  IStrongBeamIndex=0
  DO ind=1,nReflections
    RPertStrength = ABS(CUgMat(ind,1)/(TWO*RBigK*RDevPara(ind)))
    IF(RPertStrength.GE.RMinPertStrength) THEN!it's a strong beam, add it to the list
      IStrongBeamIndex= IStrongBeamIndex +1
      IStrongBeamList(IStrongBeamIndex)= ind!list of reflection ID no.s up to nStrongBeams
    ENDIF
  ENDDO
  IF(IWriteFLAG.EQ.6.AND.my_rank.EQ.0) THEN
    PRINT*, IStrongBeamIndex,"strong beams"
  ENDIF
  IF(IStrongBeamIndex+IMinWeakBeams.GT.nReflections) IErr = 1
  IF( IErr.NE.0 ) THEN
    PRINT*,"StrongAndWeakBeamDetermination(", my_rank, ") error ", IErr, &
          " Insufficient reflections to accommodate all Strong and Weak Beams"
    RETURN
  ENDIF
  
  !----------------------------------------------------------------------------
  !WEAK BEAMS
  !Weak beams must have a perturbation strength greater than 1/20 of the weakest strong beam
  !This value based on a convergence test 30 June 2016 using GaAs[110], 000 beam
  RMinPertStrength=0.05*RMinPertStrength
  IWeakBeamIndex = 0
  IWeakBeamList = 0
  DO ind=1,nReflections
    IF(MINVAL(ABS(IStrongBeamList-ind)).NE.0) THEN!it's not a strong beam
      RPertStrength = ABS(CUgMat(ind,1)/(TWO*RBigK*RDevPara(ind)))
	  IF(RPertStrength.GE.RMinPertStrength) THEN!but it is strong enough to be a weak beam
        IWeakBeamIndex= IWeakBeamIndex +1
        IWeakBeamList(IWeakBeamIndex)= ind	    
      ENDIF
    ENDIF
  ENDDO
  IF(IWriteFLAG.EQ.6.AND.my_rank.EQ.0) THEN
    PRINT*, IWeakBeamIndex,"weak beams"
  ENDIF

END SUBROUTINE StrongAndWeakBeamsDetermination
