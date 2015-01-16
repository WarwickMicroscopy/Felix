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

SUBROUTINE BlochCoefficientCalculation(ind,jnd,gnd,ILocalPixelCountMin,IErr)
  
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
  
  INTEGER(IKIND) ind,jnd,hnd,knd,gnd,pnd,&
       ierr,IThickness, &
       IThicknessIndex, ILowerLimit, &
       IUpperLimit,ILocalPixelCountMin
       
  REAL(RKIND) :: &
       Rx0,Ry0, RThickness,RKn

  CHARACTER*40 surname

  COMPLEX(CKIND) sumC,sumD
  
  COMPLEX(CKIND), DIMENSION(:,:), ALLOCATABLE :: &
       CGeneralSolutionMatrix, CGeneralEigenVectors
  COMPLEX(CKIND),DIMENSION(:),ALLOCATABLE :: &
       CGeneralEigenValues
  REAL(RKIND),DIMENSION(IReflectOut) :: &
       RPreviousWaveIntensity

   IF (my_rank.EQ.0) THEN
      DO WHILE (IMessageCounter .LT.1)
         CALL Message("BlochCoefficientCalculation",IMust,IErr)
         CALL Message("BlochCoefficientCalculation",IMust+IDebug,IErr, & 
              MessageString = "is looping, and calling subroutines itself, They are:")
         IMessageCounter = IMessageCounter +1
      END DO
   END IF

  Rx0=(ind-IPixelCount-0.5D0)*RDeltaK ! x-position in the disk
  
  Ry0=(jnd-IPixelCount-0.5D0)*RDeltaK ! y-position in the disk
    
  ! we are inside the mask
  IPixelComputed= IPixelComputed + 1

  !--------------------------------------------------------------------
  ! protocol progress
  !--------------------------------------------------------------------
  
  !!$   Displays Pixel currently working on
  IF((IWriteFLAG.GE.10.AND.IWriteFLAG.LT.100).OR.IWriteFLAG.GE.110) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, "): working on pixel (", ind, ",", jnd,") of (", &
          2*IPixelCount, ",", 2*IPixelCount, ") in total."
  ENDIF
       
  !--------------------------------------------------------------------
  ! calculate deviation parameter Sg for the tilted Ewald spheres
  !--------------------------------------------------------------------
  
  ! TiltedK used to be called Kprime2
  ! the vector of the incoming tilted beam

  CALL KVectorsCalculation(Rx0,Ry0,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " In Calculation of KVectors"
     RETURN
  ENDIF

  RKn = DOT_PRODUCT(RTiltedK,RNormDirM)
  !PRINT*,"RKn = ",RKn

  ! Compute the deviation parameter for ALL reflections
  ! within RBSMaxGVecAmp

  CALL DeviationParameterCalculation(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in Calculation of Deviation Parameter"
     RETURN
  ENDIF

  ! select only those beams where the Ewald sphere is close to the
  ! reciprocal lattice, i.e. within RBSMaxDeviationPara

  CALL StrongAndWeakBeamsDetermination(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in Determination of Strong and Weak beams"
     RETURN
  ENDIF
 
  ! select the highest reflection that corresponds to a strong beam
  nBeams= IStrongBeamIndex

  !PRINT*,nbeams

  !--------------------------------------------------------------------
  ! ALLOCATE memory for eigen problem
  !--------------------------------------------------------------------
  
  !Eigen Problem Solving
  ALLOCATE( &
       CBeamProjectionMatrix(nBeams,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CBeamProjectionMatrix"
     RETURN
  ENDIF
  ALLOCATE( &
       CDummyBeamMatrix(nBeams,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CDummyBeamMatrix"
     RETURN
  ENDIF
  ALLOCATE( &
       CUgMatEffective(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CUgMatEffective"
     RETURN
  ENDIF
  
  ALLOCATE( & 
       CEigenVectors(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CEigenVectors"
     RETURN
  ENDIF

  ALLOCATE( &
       CEigenValues(nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CEigenValues"
     RETURN
  ENDIF
  ALLOCATE( &
       CInvertedEigenVectors(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CInvertedEigenVectors"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       CAlphaWeightingCoefficients(nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CAlphaWeightingCoefficients"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       CEigenValueDependentTerms(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CEigenValueDependentTerms"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       CWaveFunctions(nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CWaveFunctions"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       RWaveIntensity(nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RWaveIntensity"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       CPsi0(nBeams), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CPsi0"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     RETURN
  ENDIF

  !--------------------------------------------------------------------
  ! construct the effective UgMat (for strong beams only at the moment)
  !--------------------------------------------------------------------
  
  IF((IWriteFLAG.GE.10.AND.IWriteFLAG.LT.100).OR.IWriteFLAG.GE.110) THEN 
     PRINT*,"BlochCoefficientCalculation(", my_rank, &
          ") using n(Strong)Beams= ", nBeams, " beams", &
          " with nWeakBeams=", IWeakBeamIndex
  ENDIF
  
  !--------------------------------------------------------------------
  ! back to eigen problem solution
  !--------------------------------------------------------------------
  
  ! compute the effective Ug matrix by selecting only those beams
  ! for which IStrongBeamList has an entry
  
  CBeamProjectionMatrix= CZERO
  DO knd=1,nBeams
     CBeamProjectionMatrix(knd,IStrongBeamList(knd))=CONE
  ENDDO

  CUgMatEffective = CZERO
  
  CUgMatEffective= &
       MATMUL( &
       CBeamProjectionMatrix, &
       MATMUL(CUgMat,TRANSPOSE(CBeamProjectionMatrix)) &
       )
    
  IF (IZolzFLAG.EQ.0) THEN

     DO hnd=1,nBeams
        CUgMatEffective(hnd,hnd) = CUgMatEffective(hnd,hnd) + TWO*RBigK*RDevPara(IStrongBeamList(hnd))
     ENDDO
     DO knd =1,nBeams ! Columns
        
        DO hnd = 1,nBeams ! Rows

           CUgMatEffective(knd,hnd) = CUgMatEffective(knd,hnd) / &
                (SQRT(1+RGn(IStrongBeamList(knd))/RKn)*SQRT(1+RGn(IStrongBeamList(hnd))/RKn))
           
        END DO
     END DO
     
     CUgMatEffective = CUgMatEffective/(TWO*RBigK)
  ELSE
     
     CUgMatEffective = CUgMatEffective/(TWO*RBigK)
     
      ! set the diagonal parts of the matrix to be equal to 
     ! strong beam deviation parameters (*2 BigK) 
     DO hnd=1,nBeams
        CUgMatEffective(hnd,hnd) = CUgMatEffective(hnd,hnd)+RDevPara(IStrongBeamList(hnd))
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
                (4*RBigK*RBigK*RDevPara(IWeakBeamList(hnd)))
           
           sumD = sumD + &
                REAL(CUgMat(IStrongBeamList(knd),IWeakBeamList(hnd))) * &
                REAL(CUgMat(IWeakBeamList(hnd),IStrongBeamList(knd))) / &
                (4*RBigK*RBigK*RDevPara(IWeakBeamList(hnd)))
           
        ENDDO
        CUgMatEffective(knd,1)= CUgMatEffective(knd,1) - sumC
        CUgMatEffective(knd,knd)= CUgMatEffective(knd,knd) - sumD
     ENDDO
     
  END IF
   
  !--------------------------------------------------------------------
  ! diagonalize the UgMatEffective
  !--------------------------------------------------------------------
   
  IF (IZolzFLAG.EQ.0) THEN
     CALL EigenSpectrum(nBeams, &
          CUgMatEffective, &
          CEigenValues(:), CEigenVectors(:,:), &
          IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in EigenSpectrum()"
        RETURN
     ENDIF
     CEigenValues = CEigenValues * RKn/RBigK
     DO knd = 1,nBeams
        CEigenVectors(knd,:) = CEigenVectors(knd,:) / &
             SQRT(1+RGn(IStrongBeamList(knd))/RKn)
     END DO
  ELSE
     CALL EigenSpectrum(nBeams, &
          CUgMatEffective, &
          CEigenValues(:), CEigenVectors(:,:), &
          IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in EigenSpectrum()"
        RETURN
     ENDIF
  END IF
 
  DO IThicknessIndex=1,IThicknessCount,1
     
     RThickness = RInitialThickness + (IThicknessIndex-1)*RDeltaThickness 
     IThickness = RInitialThickness + (IThicknessIndex-1)*RDeltaThickness 
     
     !PRINT*,"RThickness = ",RThickness

     CALL CreateWaveFunctions(rthickness,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in CreateWavefunction()"
        RETURN
     ENDIF

     IF (IOutputFLAG.GE.3) THEN
        IF(IPixelComputed.EQ.1) THEN
           ! wave functions
           CALL OpenData_MPI(IChOutWF_MPI, "WF", surname, IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in OpenDataMPI()"
              RETURN
           ENDIF
        ELSE
           
           ! wave functions
           CALL OpenDataForAppend_MPI(IChOutWF_MPI, "WF", surname, IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in OpenDataForAppend_ MPI()"
              RETURN
           ENDIF
        END IF
        CALL WriteDataC_MPI(IChOutWF_MPI, ind,jnd, &
             CFullWaveFunctions(:), &
             nReflections, 1, IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in WriteDataC_MPI of IChOutWF()"
           RETURN
        ENDIF
        CALL MPI_FILE_CLOSE(IChOutWF_MPI, IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in MPI_FILE_CLOSE of IChOutWF()"
           RETURN
        ENDIF
        
     END IF
     
     !Collection Wave Intensities from all thickness for later writing
     
     IF(IHKLSelectFLAG.EQ.0) THEN
        
        IF(IImageFLAG.LE.2) THEN
           RIndividualReflections(1:IReflectOut,IThicknessIndex,(gnd-ILocalPixelCountMin)+1) = &
                RFullWaveIntensity(1:IReflectOut)
        ELSE
           CAmplitudeandPhase(1:IReflectOut,IThicknessIndex,(gnd-ILocalPixelCountMin)+1) = &
                CFullWavefunctions(1:IReflectOut)
        END IF
     ELSE
        
        IF(IImageFLAG.LE.2) THEN
           DO pnd = 1,IReflectOut
              RIndividualReflections(pnd,IThicknessIndex,(gnd-ILocalPixelCountMin)+1) = &
                   RFullWaveIntensity(IOutputReflections(pnd))
           END DO
        ELSE
           DO pnd = 1,IReflectOut
              CAmplitudeandPhase(pnd,IThicknessIndex,(gnd-ILocalPixelCountMin)+1) = &
                   CFullWavefunctions(IOutputReflections(pnd))
           END DO
        END IF
     END IF
  END DO
  
  
  
  !--------------------------------------------------------------------
  ! OUTPUT EIGENsystem data for given pixel
  !--------------------------------------------------------------------
  
  IMAXCBuffer = 2*13*SIZE(CEigenVectors)+2*13*SIZE(CEigenValues)+3*SIZE(IStrongBeamList)+3*6*ADD_OUT_INFO
  
  IF(IOutputFLAG.GE.2) THEN
     CALL WriteEigenSystem_MPI(IChOutES_MPI, ind,jnd,nReflections,nBeams, &
          CEigenValues,CEigenVectors, IStrongBeamList,nBeams,nBeams, 1, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in WriteEigenSystem_ MPI()"
        RETURN
     ENDIF
  ENDIF
  
  IMAXCBuffer = 2*14*SIZE(CUgMatEffective)+7*6*ADD_OUT_INFO
  
  !--------------------------------------------------------------------
  ! DEALLOCATE eigen problem memory
  !--------------------------------------------------------------------
  
  DEALLOCATE( &
       CUgMatEffective,CPsi0,&
       CInvertedEigenVectors, CAlphaWeightingCoefficients, &
       CEigenValues,CEigenVectors,CEigenValueDependentTerms, &
       CBeamProjectionMatrix, CDummyBeamMatrix,CWavefunctions, &
       RWaveIntensity,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in Deallocation()"
     RETURN
  ENDIF
  
END SUBROUTINE BlochCoefficientCalculation

SUBROUTINE CreateWaveFunctions(rthickness,IErr)

  USE WriteToScreen
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER ind,jnd,knd,hnd,ierr, ifullind, iuniind,gnd,ichnk
  REAL(RKIND) rthickness 
  COMPLEX(CKIND),DIMENSION(:,:),ALLOCATABLE :: &
       CDummyEigenVectors

  IF (my_rank.EQ.0) THEN
     DO WHILE (IMessageCounter .LT.6)
        CALL Message("CreateWaveFunctions",IMust,IErr)
        IMessageCounter = IMessageCounter +1
     END DO
  END IF

  
  !--------------------------------------------------------------------
  ! calculate wavefunctions
  !--------------------------------------------------------------------
  
  CPsi0= CZERO
  IF(nBeams .GE. 0) CPsi0(1) = CONE
  
  ALLOCATE( &
       CDummyEigenVectors(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateWavefunctions(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CDummyEigenVectors"
     RETURN
  ENDIF
  
  ! Invert the EigenVector matrix
  !InvertedEigenVectors= CONJG(TRANSPOSE(CEigenVectors(:,:)))
  CDummyEigenVectors = CEigenVectors

  CALL INVERT(nBeams,CDummyEigenVectors(:,:),CInvertedEigenVectors,IErr)

  !PRINT*,MATMUL(CEigenVectors,CInvertedEigenVectors)
  
  !From EQ 6.32 in Kirkland Advance Computing in EM
  CAlphaWeightingCoefficients = MATMUL(CInvertedEigenVectors(1:nBeams,1:nBeams),CPsi0) 
  
  CEigenValueDependentTerms= CZERO
 
  DO hnd=1,nBeams !IReflectOut 
     
     ! This needs to be a diagonal matrix
     CEigenValueDependentTerms(hnd,hnd) = &
          EXP(CIMAGONE*CMPLX(RThickness,ZERO,CKIND)*CEigenValues(hnd)) 
     
  ENDDO
  
  ! EQ 6.35 in Kirkland Advance Computing in EM
  ! C-1*C*alpha 
  
  CWaveFunctions(:) = &
       MATMUL( &
       MATMUL(CEigenVectors(1:nBeams,1:nBeams),CEigenValueDependentTerms), & 
       CAlphaWeightingCoefficients(:) &
       )
  
  DO hnd=1,nBeams
     RWaveIntensity(hnd)= &
          CONJG(CWaveFunctions(hnd)) * CWaveFunctions(hnd)
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
  
  DEALLOCATE(&
       CDummyEigenVectors, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CreateWavefunctions(", my_rank, ") error ", IErr, &
          " in DEALLOCATE() of DYNAMIC variables CDummyEigenVectors"
     RETURN
  ENDIF
  
END SUBROUTINE CreateWavefunctions


!!$Calculates the x,y,z components of the incident tilted k_vector
SUBROUTINE KVectorsCalculation(Rx0,Ry0,IErr)
  
  USE WriteToScreen
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  REAL(RKIND) Rx0,Ry0
  INTEGER(IKIND) :: IErr

  IF (my_rank.EQ.0) THEN
     DO WHILE (IMessageCounter .LT.2)
        CALL Message("KVectorsCalculation",IMust,IErr)
        IMessageCounter = IMessageCounter +1
     END DO
  END IF

  !!$  k_x - based on crystal orientation
  RTiltedK(1)= Rx0
  !!$  k_y - based on crystal orientation
  RTiltedK(2)= Ry0
  !!$  k_z - taken from: k_z = (k_y)^2 + (k_x)^2 + (k_z)^2 = K^2  
  RTiltedK(3)= SQRT(RBigK**2 - Rx0**2 - Ry0**2)
  
END SUBROUTINE KVectorsCalculation

SUBROUTINE DeviationParameterCalculation(IErr)

  USE WriteToScreen
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  INTEGER(IKIND) knd,IErr
  
  IF (my_rank.EQ.0) THEN
     DO WHILE (IMessageCounter .LT.3)
        CALL Message("DeviationParameterCalculation",IMust,IErr)
        IMessageCounter = IMessageCounter +1
     END DO
  END IF

  DO knd=1,nReflections
     ! DevPara used to be called Sg in the book
     
     RDevPara(knd)= &
          -( RBigK + DOT_PRODUCT(RgVecMatT(knd,:),RTiltedK(:)) /RBigK) + &
          SQRT( &
          ( RBigK**2 + DOT_PRODUCT(RgVecMatT(knd,:),RTiltedK(:)) )**2 /RBigK**2 - &
          (RgVecMag(knd)**2 + &
          2.0D0* DOT_PRODUCT(RgVecMatT(knd,:),RTiltedK(:))) &
          )
  END DO

END SUBROUTINE DeviationParameterCalculation

SUBROUTINE StrongAndWeakBeamsDetermination(IErr)
  
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI
  
  INTEGER(IKIND) ind,knd,IErr,IMinimum,IMaximum,ICheck,jnd,hnd, &
       IAdditionalBmaxStrongBeams,IAdditionalPmaxStrongBeams,&
       IBeamIterationCounter,IFound
  REAL(RKIND) RDummySg(nReflections)
  REAL(RKIND) sumC
  INTEGER(IKIND), DIMENSION(:),ALLOCATABLE  :: &
       IAdditionalBmaxStrongBeamList,IAdditionalPmaxStrongBeamList

    IF (my_rank.EQ.0) THEN
     DO WHILE (IMessageCounter .LT.4)
        CALL Message("StrongAndWeakBeamsDetermination",IMust,IErr)
        IMessageCounter = IMessageCounter +1
     END DO
  END IF

  !----------------------------------------------------------------------------
  ! Determine RBSMaxDeviationPara
  !----------------------------------------------------------------------------

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
  
  RDummySg = ABS(RMeanInnerCrystalPotential/RDevPara)
  
  jnd = 0

  !Determine RBSBethePara

  !IAdditionalStrongBeams = 1

  !IBeamIterationCounter = 0

  !DO WHILE(IAdditionalStrongBeams.NE.0)
  !  IBeamIterationCounter = IBeamIterationCounter +1

  !----------------------------------------------------------------------------
  ! Apply Bmax Criteria 
  !----------------------------------------------------------------------------

  IAdditionalBmaxStrongBeams = 0
  IAdditionalPmaxStrongBeams = 0
  
  IF(IStrongBeamIndex+IMinWeakBeams.GT.nReflections) THEN
     IErr = 1
  END IF
  IF( IErr.NE.0 ) THEN
     PRINT*,"StrongAndWeakBeamDetermination(", my_rank, ") error ", IErr, &
          " Insufficient reflections to accommodate all Strong and Weak Beams"
     RETURN
  ENDIF
    
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
          (ABS(RMeanInnerCrystalPotential/RDevPara(knd)) .GE. RBSBethePara)) THEN
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
  
  ALLOCATE(&
       IAdditionalBmaxStrongBeamList(IWeakBeamIndex),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StrongAndWeakBeamsDetermination(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables IAdditionalBmaxStrongBeamList"
     RETURN
  ENDIF
  
  IAdditionalBmaxStrongBeamList = 0
  
  DO knd = 2,IStrongBeamIndex
     DO hnd = 1,IWeakBeamIndex
        IFound = 0
        sumC = ZERO
        sumC = sumC + &
             REAL(CUgMat(IStrongBeamList(knd),IWeakBeamList(hnd)))* &
             REAL(CUgMat(IWeakBeamList(hnd),1)) / &
             (2*RBigK*RDevPara(IWeakBeamList(hnd)))
        
        sumC = sumC/REAL(CUgMat(IStrongBeamList(knd),1))
        
        !PRINT*,"Perturbation Strength = ",REAL(CUgMat(IStrongBeamList(knd),IWeakBeamList(hnd))), &
        !     RDevPara(IWeakBeamList(hnd)),REAL(CUgMat(IWeakBeamList(hnd),1)),sumC,REAL(sumC) 
        
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
  !DEALLOCATE(&
  !     IAdditionalBmaxStrongBeamList, &
  !     STAT=IErr)
  !IF( IErr.NE.0 ) THEN
  !   PRINT*,"StrongAndWeakBeamDetermination(", my_rank, ") error ", IErr, &
  !        " in DEALLOCATE() of DYNAMIC variables IAdditionalBmaxStrongBeamList"
  !   RETURN
  !ENDIF  
  
  ! END DO

  !----------------------------------------------------------------------------
  ! Apply Pmax Criteria 
  !----------------------------------------------------------------------------
  
  IWeakBeamIndex=0
  IWeakBeamList = 0
  
  DO knd=1,nReflections
     IFound = 0
     IF( (ABS(RDevPara(knd)) .GT. RBSMaxDeviationPara).AND. &
          (ABS(RMeanInnerCrystalPotential/RDevPara(knd)) .GE. RBSBethePara)) THEN
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

  ALLOCATE(&
       IAdditionalPmaxStrongBeamList(IWeakBeamIndex),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StrongAndWeakBeamsDetermination(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables IAdditionalPmaxStrongBeamList"
     RETURN
  ENDIF
  
  IAdditionalPmaxStrongBeamList = 0
  
  DO knd = 1,IAdditionalBmaxStrongBeams
     DO hnd = 1,IWeakBeamIndex
        IFound = 0
        sumC = ZERO
        sumC = sumC + &
             REAL(CUgMat(IAdditionalBmaxStrongBeamList(knd),IWeakBeamList(hnd)))* &
             REAL(CUgMat(IWeakBeamList(hnd),1)) / &
             (2*RBigK*RDevPara(IWeakBeamList(hnd)))
        
        sumC = sumC/REAL(CUgMat(IAdditionalBmaxStrongBeamList(knd),1))
        
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
  
  !PRINT*,IStrongBeamIndex,IAdditionalBmaxStrongBeams,IAdditionalPmaxStrongBeams
!!$  PRINT*,IStrongBeamList

  IWeakBeamIndex=0
  IWeakBeamList = 0
  
  DO knd=1,nReflections
     IFound = 0
     IF( (ABS(RDevPara(knd)) .GT. RBSMaxDeviationPara).AND. &
          (ABS(RMeanInnerCrystalPotential/RDevPara(knd)) .GE. RBSBethePara)) THEN
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

END SUBROUTINE StrongAndWeakBeamsDetermination
