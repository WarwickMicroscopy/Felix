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
  
  INTEGER(IKIND) :: IYPixelIndex,IXPixelIndex,hnd,knd,IPixelNumber,pnd,&
       ierr,IThickness,IThicknessIndex,ILowerLimit,IUpperLimit,IFirstPixelToCalculate       
  REAL(RKIND) :: RPixelGVectorXPosition,RPixelGVectorYPosition, RThickness,RKn
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

  RPixelGVectorXPosition=(REAL(IYPixelIndex,RKIND)-REAL(IPixelCount,RKIND)-0.5_RKIND)*RDeltaK ! x-position in the disk
  RPixelGVectorYPosition=(REAL(IXPixelIndex,RKIND)-REAL(IPixelCount,RKIND)-0.5_RKIND)*RDeltaK ! y-position in the disk
    
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
  ! calculate deviation parameter Sg for the tilted Ewald spheres
  ! TiltedK is the vector of the incoming tilted beam
  !TiltedK is in units of (1/A), in the microscope ref frame(NB exp(2*pi*i*k.r), optical convention)
  RTiltedK(1)= RPixelGVectorXPosition  !!$  kx - based on crystal orientation
  RTiltedK(2)= RPixelGVectorYPosition  !!$  ky - based on crystal orientation
  RTiltedK(3)= SQRT(RBigK**2 - RPixelGVectorXPosition**2 - RPixelGVectorYPosition**2)  !!$  kz -  from: ky^2 +kx^2 +kz^2 = K^2  
  RKn = DOT_PRODUCT(RTiltedK,RNormDirM)
  
  !Compute the deviation parameter for reflection pool
  !NB RDevPara is in units of (1/A), in the microscope ref frame(NB exp(2*pi*i*k.r), optical convention)
  DO knd=1,nReflections
    RDevPara(knd)= -( RBigK + DOT_PRODUCT(RgPool(knd,:),RTiltedK(:)) /RBigK) + &
      SQRT( ( RBigK**2 + DOT_PRODUCT(RgPool(knd,:),RTiltedK(:)) )**2 /RBigK**2 - &
      (RgPoolMag(knd)**2 + TWO*DOT_PRODUCT(RgPool(knd,:),RTiltedK(:))) )
    IF(IWriteFLAG.EQ.6.AND.my_rank.EQ.0.AND.knd.EQ.2.AND.IYPixelIndex.EQ.10.AND.IXPixelIndex.EQ.10) THEN
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

  ! select the highest reflection that corresponds to a strong beam
  nBeams= IStrongBeamIndex

  !--------------------------------------------------------------------
  ! ALLOCATE memory for eigen problem
  !--------------------------------------------------------------------
  ALLOCATE(CBeamProjectionMatrix(nBeams,nReflections),STAT=IErr)
  ALLOCATE(CDummyBeamMatrix(nBeams,nReflections),STAT=IErr)
  ALLOCATE(CUgMatEffective(nBeams,nBeams),STAT=IErr)
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
  CUgMatEffective = CZERO
  CBeamTranspose=TRANSPOSE(CBeamProjectionMatrix)

  CALL ZGEMM('N','N',nReflections,nBeams,nReflections,CONE,CUgMat, &
       nReflections,CBeamTranspose,nReflections,CZERO,CUgMatPartial,nReflections)
  CALL ZGEMM('N','N',nBeams,nBeams,nReflections,CONE,CBeamProjectionMatrix, &
       nBeams,CUgMatPartial,nReflections,CZERO,CUgMatEffective,nBeams)

  IF (IZolzFLAG.EQ.0) THEN
    DO hnd=1,nBeams
      !CUgMatEffective(hnd,hnd) = CUgMatEffective(hnd,hnd) + TWO*RBigK*RDevPara(IStrongBeamList(hnd))
      CUgMatEffective(hnd,hnd) = CUgMatEffective(hnd,hnd) + TWO*TWOPI*RBigK*TWOPI*RDevPara(IStrongBeamList(hnd))
    ENDDO
    DO knd =1,nBeams ! Columns
      DO hnd = 1,nBeams ! Rows
        CUgMatEffective(knd,hnd) = CUgMatEffective(knd,hnd) / &
         !(SQRT(1+RgDotNorm(IStrongBeamList(knd))/RKn)*SQRT(1+RgDotNorm(IStrongBeamList(hnd))/RKn))
         (SQRT(1+TWOPI*RgDotNorm(IStrongBeamList(knd))/RKn)*SQRT(1+TWOPI*RgDotNorm(IStrongBeamList(hnd))/RKn))
      END DO
    END DO
    !CUgMatEffective = CUgMatEffective/(TWO*RBigK)
    CUgMatEffective = CUgMatEffective/(TWO*TWOPI*RBigK)
  ELSE
    !CUgMatEffective = CUgMatEffective/(TWO*RBigK)
    CUgMatEffective = CUgMatEffective/(TWO*TWOPI*RBigK)
    ! set the diagonal parts of the matrix to be equal to 
    ! strong beam deviation parameters (*2 BigK) 
    DO hnd=1,nBeams
      !CUgMatEffective(hnd,hnd) = CUgMatEffective(hnd,hnd)+RDevPara(IStrongBeamList(hnd))
      CUgMatEffective(hnd,hnd) = CUgMatEffective(hnd,hnd)+TWOPI*RDevPara(IStrongBeamList(hnd))
    ENDDO
    ! add the weak beams perturbatively for the 1st column (sumC) and
    ! the diagonal elements (sumD)
    DO knd=2,nBeams
      sumC= CZERO
      sumD= CZERO
      DO hnd=1,IWeakBeamIndex
        sumC = sumC + &
          REAL(CUgMat(IStrongBeamList(knd),IWeakBeamList(hnd))) * &
          REAL(CUgMat(IWeakBeamList(hnd),1)) / &
         !(4*RBigK*RBigK*RDevPara(IWeakBeamList(hnd)))
         (4*TWOPI*RBigK*TWOPI*RBigK*RDevPara(IWeakBeamList(hnd)))
        sumD = sumD + &
          REAL(CUgMat(IStrongBeamList(knd),IWeakBeamList(hnd))) * &
          REAL(CUgMat(IWeakBeamList(hnd),IStrongBeamList(knd))) / &
         !(4*TRBigK*RBigK*TWOPI*RDevPara(IWeakBeamList(hnd)))
         (4*TWOPI*RBigK*TWOPI*RBigK*TWOPI*RDevPara(IWeakBeamList(hnd)))
      ENDDO
      CUgMatEffective(knd,1)= CUgMatEffective(knd,1) - sumC
      CUgMatEffective(knd,knd)= CUgMatEffective(knd,knd) - sumD
    ENDDO
  END IF

  IF(IWriteFLAG.EQ.7.AND.my_rank.EQ.0) THEN
   PRINT*,"Effective Ug matrix"
	DO hnd =1,8
     WRITE(SPrintString,FMT='(16(1X,F5.2))') CUgMatEffective(hnd,1:8)
     PRINT*,TRIM(ADJUSTL(SPrintString))
    END DO
  END IF	
  
  !--------------------------------------------------------------------
  ! diagonalize the UgMatEffective
  IF (IZolzFLAG.EQ.0) THEN
     CALL EigenSpectrum(nBeams,CUgMatEffective,CEigenValues(:), CEigenVectors(:,:),IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in EigenSpectrum()"
        RETURN
     END IF
!     CEigenValues = CEigenValues * RKn/RBigK
     CEigenValues = CEigenValues * RKn/RBigK
     DO knd = 1,nBeams
        CEigenVectors(knd,:) = CEigenVectors(knd,:) / &
!             SQRT(1+RgDotNorm(IStrongBeamList(knd))/RKn)
             SQRT(1+TWOPI*RgDotNorm(IStrongBeamList(knd))/RKn)
     END DO
  ELSE
     CALL EigenSpectrum(nBeams,CUgMatEffective,CEigenValues(:),CEigenVectors(:,:),IErr)
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
  !--------------------------------------------------------------------
  
  DEALLOCATE(CUgMatEffective,CPsi0, &
       CBeamTranspose, CUgMatPartial, &
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

!!$Calculates the x,y,z components of the incident tilted k_vector
SUBROUTINE KVectorsCalculation(RPixelGVectorXPosition,RPixelGVectorYPosition,IErr)
!This is redundant
  USE WriteToScreen
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  REAL(RKIND) :: RPixelGVectorXPosition,RPixelGVectorYPosition
  INTEGER(IKIND) :: IErr

  IF (my_rank.EQ.0) THEN
     DO WHILE (IMessageCounter .LT.2)
        CALL Message("KVectorsCalculation",IMust,IErr)
        IMessageCounter = IMessageCounter +1
     END DO
  END IF

  !!$  k_x - based on crystal orientation
  RTiltedK(1)= RPixelGVectorXPosition
  !!$  k_y - based on crystal orientation
  RTiltedK(2)= RPixelGVectorYPosition
  !!$  k_z - taken from: k_z = (k_y)^2 + (k_x)^2 + (k_z)^2 = K^2  
  RTiltedK(3)= SQRT(RBigK**2 - RPixelGVectorXPosition**2 - RPixelGVectorYPosition**2)
  
END SUBROUTINE KVectorsCalculation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DeviationParameterCalculation(IErr)
!This is redundant
  USE WriteToScreen
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  INTEGER(IKIND) :: knd,IErr
  
  IF (my_rank.EQ.0) THEN
     DO WHILE (IMessageCounter .LT.3)
        CALL Message("DeviationParameterCalculation",IMust,IErr)
        IMessageCounter = IMessageCounter +1
     END DO
  END IF

  DO knd=1,nReflections
     ! DevPara is devitaion parameter, also known as Sg 
     RDevPara(knd)= -( RBigK + DOT_PRODUCT(RgPool(knd,:),RTiltedK(:)) /RBigK) + &
          SQRT( ( RBigK**2 + DOT_PRODUCT(RgPool(knd,:),RTiltedK(:)) )**2 /RBigK**2 - &
          (RgPoolMag(knd)**2 + &
          TWO* DOT_PRODUCT(RgPool(knd,:),RTiltedK(:))) )
  END DO

END SUBROUTINE DeviationParameterCalculation

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
  REAL(RKIND) :: RDummySg(nReflections), sumC
  INTEGER(IKIND), DIMENSION(:),ALLOCATABLE  :: IAdditionalBmaxStrongBeamList,IAdditionalPmaxStrongBeamList

  IF (my_rank.EQ.0) THEN
    DO WHILE (IMessageCounter .LT.4)
      CALL Message("StrongAndWeakBeamsDetermination",IMust,IErr)
      IMessageCounter = IMessageCounter +1
     END DO
  END IF

  !----------------------------------------------------------------------------
  ! Determine RBSMaxDeviationPara
  IStrongBeamList = 0
  RDummySg = ABS(RDevPara)
  DO ind=1,IMinStrongBeams
    IMinimum = MINLOC(RDummySg,1)
    IF(ind.EQ.IMinStrongBeams) THEN
       RBSMaxDeviationPara = ABS(RDummySg(IMinimum))
    ELSE
      RDummySg(IMinimum) = 1000000 !Large number
    END IF
  END DO
  IStrongBeamIndex=0
  IWeakBeamIndex=0
  DO knd=1,nReflections
    IF( ABS(RDevPara(knd)) .LE. RBSMaxDeviationPara ) THEN
      IStrongBeamIndex= IStrongBeamIndex +1
      IStrongBeamList(IStrongBeamIndex)= knd
    ENDIF
  ENDDO
  RDummySg = ABS(RMeanInnerPotential/RDevPara)
  
  !----------------------------------------------------------------------------
  ! Apply Bmax Criteria 
  IAdditionalBmaxStrongBeams = 0
  IAdditionalPmaxStrongBeams = 0
  IF(IStrongBeamIndex+IMinWeakBeams.GT.nReflections) IErr = 1
  IF( IErr.NE.0 ) THEN
    PRINT*,"StrongAndWeakBeamDetermination(", my_rank, ") error ", IErr, &
          " Insufficient reflections to accommodate all Strong and Weak Beams"
    RETURN
  ENDIF
    
  jnd=0
  DO ind=1,nReflections
    ICheck = 0
    IMaximum = MAXLOC(RDummySg,1)
    DO knd = 1,IStrongBeamIndex
      IF(IMaximum.EQ.IStrongBeamList(knd)) THEN
        ICheck = 1
        EXIT
      END IF
    END DO
    IF(ICheck.EQ.0) THEN
      jnd = jnd+1
    END IF
    IF(jnd.EQ.IMinWeakBeams) THEN
      RBSBethePara = (RDummySg(IMaximum))
    ELSE
      RDummySg(IMaximum) = 0.D0 !Large number
    END IF
  END DO

  IWeakBeamIndex=0
  IWeakBeamList = 0
  DO knd=1,nReflections
    IFound = 0
    IF( (ABS(RDevPara(knd)) .GT. RBSMaxDeviationPara).AND. &
        (ABS(RMeanInnerPotential/RDevPara(knd)) .GE. RBSBethePara)) THEN
      DO ind = 1,IStrongBeamIndex
        IF(IStrongBeamList(ind).EQ.knd) THEN
          IFound=IFound+1
        END IF
      END DO
      IF(IFound.EQ.0) THEN
        IWeakBeamIndex= IWeakBeamIndex +1
        IWeakBeamList(IWeakBeamIndex)= knd
      END IF
      IFound = 0
    ENDIF
  ENDDO

  ALLOCATE(IAdditionalBmaxStrongBeamList(IWeakBeamIndex),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StrongAndWeakBeamsDetermination(",my_rank,")error allocating IAdditionalBmaxStrongBeamList"
     RETURN
  ENDIF
  
  IAdditionalBmaxStrongBeamList = 0
  DO knd = 2,IStrongBeamIndex
    DO hnd = 1,IWeakBeamIndex
      IFound = 0
      sumC = ZERO
      sumC = sumC + REAL(CUgMatNoAbs(IStrongBeamList(knd),IWeakBeamList(hnd)))* &
             REAL(CUgMatNoAbs(IWeakBeamList(hnd),1)) / &
             (2*RBigK*RDevPara(IWeakBeamList(hnd)))
      sumC = sumC/REAL(CUgMatNoAbs(IStrongBeamList(knd),1))
      IF(ABS(sumC).GE.RBSBmax) THEN
        DO ind =1,IWeakBeamIndex
          IF(IAdditionalBmaxStrongBeamList(ind).EQ.IWeakBeamList(hnd)) THEN
            IFound = IFound+1
          END IF
        END DO
        IF(IFound.EQ.0) THEN
          IAdditionalBmaxStrongBeams = IAdditionalBmaxStrongBeams + 1
          IAdditionalBmaxStrongBeamList(IAdditionalBmaxStrongBeams) = IWeakBeamList(hnd)
        END IF
      END IF
    END DO
  END DO
  
  IF(IAdditionalBmaxStrongBeams.NE.0) THEN
     IStrongBeamList((IStrongBeamIndex+1):(IStrongBeamIndex+IAdditionalBmaxStrongBeams)) = &
          IAdditionalBmaxStrongBeamList(:IAdditionalBmaxStrongBeams)
     IStrongBeamIndex = IStrongBeamIndex  +IAdditionalBmaxStrongBeams
  END IF

  !----------------------------------------------------------------------------
  ! Apply Pmax Criteria 
  IWeakBeamIndex=0
  IWeakBeamList = 0
  DO knd=1,nReflections
    IFound = 0
    IF( (ABS(RDevPara(knd)) .GT. RBSMaxDeviationPara).AND. &
        (ABS(RMeanInnerPotential/RDevPara(knd)) .GE. RBSBethePara)) THEN
      DO ind = 1,IStrongBeamIndex
        IF(IStrongBeamList(ind).EQ.knd) THEN
          IFound = IFound + 1
        END IF
      END DO
      IF(IFound.EQ.0) THEN
        IWeakBeamIndex= IWeakBeamIndex +1
        IWeakBeamList(IWeakBeamIndex)= knd
      END IF
      IFound = 0
    ENDIF
  ENDDO

  ALLOCATE(IAdditionalPmaxStrongBeamList(IWeakBeamIndex),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StrongAndWeakBeamsDetermination(",my_rank,")error allocating PmaxStrongBeamList"
     RETURN
  ENDIF
  
  IAdditionalPmaxStrongBeamList = 0
  DO knd = 1,IAdditionalBmaxStrongBeams
    DO hnd = 1,IWeakBeamIndex
      IFound = 0
      sumC = ZERO
      sumC = sumC +REAL(CUgMatNoAbs(IAdditionalBmaxStrongBeamList(knd),IWeakBeamList(hnd)))* &
             REAL(CUgMatNoAbs(IWeakBeamList(hnd),1)) / &
             (2*RBigK*RDevPara(IWeakBeamList(hnd)))
      sumC = sumC/REAL(CUgMatNoAbs(IAdditionalBmaxStrongBeamList(knd),1))
      IF(ABS(REAL(sumC)).GE.RBSPmax) THEN
        DO ind =1,IWeakBeamIndex
          IF(IAdditionalPmaxStrongBeamList(ind).EQ.IWeakBeamList(hnd)) THEN
            IFound = IFound+1
          END IF
        END DO
        IF(IFound.EQ.0) THEN
          IAdditionalPmaxStrongBeams = IAdditionalPmaxStrongBeams + 1
          IAdditionalPmaxStrongBeamList(IAdditionalPmaxStrongBeams) = IWeakBeamList(hnd)
        END IF
      END IF
    END DO
  END DO
  
  IF(IAdditionalPmaxStrongBeams.NE.0) THEN
    IStrongBeamList((IStrongBeamIndex+1):(IStrongBeamIndex+IAdditionalPmaxStrongBeams)) = &
          IAdditionalPmaxStrongBeamList(:IAdditionalPmaxStrongBeams)
    IStrongBeamIndex = IStrongBeamIndex  + IAdditionalPmaxStrongBeams
  END IF

  IWeakBeamIndex = 0
  IWeakBeamList = 0
  DO knd=1,nReflections
    IFound = 0
    IF( (ABS(RDevPara(knd)) .GT. RBSMaxDeviationPara).AND. &
        (ABS(RMeanInnerPotential/RDevPara(knd)) .GE. RBSBethePara)) THEN
       DO ind = 1,IStrongBeamIndex
         IF(IStrongBeamList(ind).EQ.knd) THEN
           IFound = IFound + 1
         END IF
      END DO
      IF(IFound.EQ.0) THEN
        IWeakBeamIndex= IWeakBeamIndex +1
        IWeakBeamList(IWeakBeamIndex)= knd
      END IF
      IFound = 0
    ENDIF
  ENDDO

  DEALLOCATE(IAdditionalBmaxStrongBeamList,STAT=IErr)
  DEALLOCATE(IAdditionalPmaxStrongBeamList,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StrongAndWeakBeamsDetermination(",my_rank,")error in deallocations"
     RETURN
  ENDIF

END SUBROUTINE StrongAndWeakBeamsDetermination
