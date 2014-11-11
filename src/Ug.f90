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

SUBROUTINE GMatrixInitialisation (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara

  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,ierr,IUniqueKey,knd,IFound

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"GMatrixInitialisation(",my_rank,")"
  END IF

  DO ind=1,nReflections
     DO jnd=1,nReflections
        RgMatMat(ind,jnd,:)= RgVecMatT(ind,:)-RgVecMatT(jnd,:)
        RgMatMag(ind,jnd)= SQRT(DOT_PRODUCT(RgMatMat(ind,jnd,:),RgMatMat(ind,jnd,:)))
     ENDDO
  ENDDO
   
  RgMatMag = RgMatMag/TWOPI
  
END SUBROUTINE GMatrixInitialisation

SUBROUTINE SymmetryRelatedStructureFactorDetermination (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,ierr,knd,Iuid


  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"SymmetryRelatedStructureFactorDetermination(",my_rank,")"
  END IF

  !Immediately set all the zeros to Relation 1
  
  ISymmetryRelations = 0_IKIND
  
  Iuid = 0
  
  Iuid = Iuid + 1

  WHERE (ABS(CUgMat).LE.RTolerance)
     ISymmetryRelations = Iuid
  END WHERE
  
  DO ind = 1,nReflections
     DO jnd = 1,ind
        IF(ISymmetryRelations(ind,jnd).NE.0) THEN
           CYCLE
        ELSE
           Iuid = Iuid + 1
           WHERE (ABS(ABS(CUgMat)-ABS(CUgMat(ind,jnd))).LE.RTolerance)
              ISymmetryRelations = Iuid
           END WHERE
        END IF
     END DO
  END DO

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Unique Ugs = ",Iuid
  END IF

  ALLOCATE(&
       ISymmetryStrengthKey(Iuid,2),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"UgCalculation(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE(&
       CSymmetryStrengthKey(Iuid),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"UgCalculation(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  
END SUBROUTINE SymmetryRelatedStructureFactorDetermination

!---------------------------------------------------------------------
SUBROUTINE StructureFactorInitialisation (IErr,CZeroMat)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind, jnd, knd, oddindlorentz, evenindlorentz, oddindgauss, &
       evenindgauss,imaxj, IFound, ICount, currentatom,IErr
  INTEGER(IKIND),DIMENSION(2) :: &
       IPos, ILoc
  COMPLEX(CKIND) CVgij
  REAL(RKIND) :: &
       RMeanInnerPotentialVolts,RAtomicFormFactor, Lorentzian,Gaussian
 ! REAL(RKIND),DIMENSION(3) :: RAtomicFormFactorSum
  COMPLEX(CKIND),DIMENSION(:,:), ALLOCATABLE :: &
       CZeroMat

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"StructureFactorInitialisation(",my_rank,")"
  END IF
  

  DO ind=1,nReflections
     imaxj = ind
     DO jnd=1,imaxj 
        
        CVgij= 0.0D0
        
        DO iAtom=1, INAtomsUnitCell
           ICurrentAtom = IAtoms(IAtom)
           ! calculate f_e(q) as in Eq. (C.15) of Kirkland, "Advanced Computing in EM"
           
           SELECT CASE (IScatterFactorMethodFLAG)
              
           CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians
 
              RAtomicFormFactor = ZERO
              DO knd = 1, 3
              
                 !odd and even indicies for Lorentzian function
                 evenindlorentz = knd*2
                 oddindlorentz = knd*2 -1
                 
                 !odd and even indicies for Gaussian function
                 evenindgauss = evenindlorentz + 6
                 oddindgauss = oddindlorentz + 6
                 
                 !Kirkland Method uses summation of 3 Gaussians and 3 Lorentzians (summed in loop)
                 RAtomicFormFactor = RAtomicFormFactor + &
                      !3 Lorentzians
                      LORENTZIAN(RScattFactors(ICurrentAtom,oddindlorentz), RGMatMag(ind,jnd),ZERO,&
                      RScattFactors(ICurrentAtom,evenindlorentz))+ &
                      !3 Gaussians
                      GAUSSIAN(RScattFactors(ICurrentAtom,oddindgauss),RgMatMag(ind,jnd),ZERO, & 
                      1/(SQRT(2*RScattFactors(ICurrentAtom,evenindgauss))),ZERO)

              END DO
              
!----------------------------------------------------------------------------------------------
              !Old Form 
              ! RAtomicFormFactor = &
                      ! 3 Lorentzians
                 !     LORENTZIAN(RScattFactors(ICurrentAtom,1), RGMatMag(ind,jnd),ZERO,&
                  !    RScattFactors(ICurrentAtom,evenind)) + &
                   !   LORENTZIAN(RScattFactors(ICurrentAtom,3), RGMatMag(ind,jnd),ZERO,&
                   !   RScattFactors(ICurrentAtom,4)) + &
                  !     LORENTZIAN(RScattFactors(ICurrentAtom,5), RGMatMag(ind,jnd),ZERO,&
!                       RScattFactors(ICurrentAtom,6)) + &
!                       ! 3 Gaussians
!                       Gaussian(RScattFactors(ICurrentAtom,7),ZERO, & 
!                       1/(SQRT(2*RScattFactors(ICurrentAtom,8))),RgMatMag(ind,jnd),ZERO) + &
!                       Gaussian(RScattFactors(ICurrentAtom,9),ZERO, & 
!                       1/(SQRT(2*RScattFactors(ICurrentAtom,10))),RgMatMag(ind,jnd),ZERO) + &
!                       Gaussian(RScattFactors(ICurrentAtom,11),ZERO, & 
!                       1/(SQRT(2*RScattFactors(ICurrentAtom,12))),RgMatMag(ind,jnd),ZERO)
                 

             
                 !  RScattFactors(ICurrentAtom,7) * &
                 !  EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(ICurrentAtom,8)) + &
                 !  RScattFactors(ICurrentAtom,9) * &
                 !  EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(ICurrentAtom,10)) + &
                 !  RScattFactors(ICurrentAtom,11) * &
                 !  EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(ICurrentAtom,12))
!----------------------------------------------------------------------------------------------
              
           CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
              
               RAtomicFormFactor = ZERO
               DO knd = 1, 4
                 

                  !Peng Method uses summation of 4 Gaussians
                  RAtomicFormFactor = RAtomicFormFactor + &
                       !4 Gaussians
                       GAUSSIAN(RScattFactors(ICurrentAtom,knd),RgMatMag(ind,jnd),ZERO, & 
                       SQRT(2/RScattFactors(ICurrentAtom,knd+4)),ZERO)

               END DO
!----------------------------------------------------------------------------------------------------
              !Old Form
           !    RAtomicFormFactor = &
!                    RScattFactors(ICurrentAtom,1) * &
!                    EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,5)) + &
!                    RScattFactors(ICurrentAtom,2) * &
!                    EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!                    RScattFactors(ICurrentAtom,3) * &
!                    EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,7)) + &
!                    RScattFactors(ICurrentAtom,4) * &
!                    EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,8))
!----------------------------------------------------------------------------------------------------

           CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)

               RAtomicFormFactor = ZERO
               DO knd = 1, 4
                 
                  evenindgauss = knd*2
                  oddindgauss = knd*2 -1
                  !Doyle &Turner uses summation of 4 Gaussians
                  RAtomicFormFactor = RAtomicFormFactor + &
                       !4 Gaussians
                       GAUSSIAN(RScattFactors(ICurrentAtom,oddindgauss),RgMatMag(ind,jnd),ZERO, & 
                       SQRT(2/RScattFactors(ICurrentAtom,evenindgauss)),ZERO)

               END DO

!This doesn't produce a correct result?
             ! RAtomicFormFactor = &
             !      RScattFactors(ICurrentAtom,1) * &
             !      EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,2)) + &
             !      RScattFactors(ICurrentAtom,3) * &
             !      EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,4)) + &
             !      RScattFactors(ICurrentAtom,5) * &
             !      EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
             !      RScattFactors(ICurrentAtom,7) * &
             !      EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,8))

           END SELECT
              
           ! initialize potential as in Eq. (6.10) of Kirkland

           RAtomicFormFactor = RAtomicFormFactor*ROcc(iAtom)

           IF (IAnisoDebyeWallerFactorFlag.EQ.0) THEN
              
              IF(RDWF(iAtom).GT.10.OR.RDWF(iAtom).LT.0) THEN
                 RDWF(iAtom) = RDebyeWallerConstant
              END IF
              
              SELECT CASE (IScatterFactorMethodFLAG)

              CASE (0)

              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-((RgMatMag(ind,jnd)/2.D0)**2)*RDWF(iAtom))
                               
              CASE(1)

              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-((RgMatMag(ind,jnd)/2.D0)**2)*RDWF(iAtom))
                          
                               
              CASE(2)

              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-((RgMatMag(ind,jnd)/2.D0)**2)*RDWF(iAtom))
                               
              END SELECT
              
              !RAtomicFormFactor = RAtomicFormFactor * &
              !     EXP(-RgMatMag(ind,jnd)**2*RDWF(iAtom))
              
           ELSE
              
              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-TWOPI*DOT_PRODUCT(RgMatMat(ind,jnd,:), &
                   MATMUL( RAnisotropicDebyeWallerFactorTensor( &
                   IAnisoDWFT(iAtom),:,:), &
                   RgMatMat(ind,jnd,:))))
              
           END IF
           
           CVgij = CVgij + &
                RAtomicFormFactor * &
                EXP(-CIMAGONE* &
                DOT_PRODUCT(RgMatMat(ind,jnd,:), RrVecMat(iAtom,:)) &
                )
        ENDDO

        CUgMat(ind,jnd)=((((TWOPI**2)* RRelativisticCorrection) / &
             (PI * RVolume)) * CVgij)
              
     ENDDO
  ENDDO

  RMeanInnerCrystalPotential= REAL(CUgMat(1,1))
  RMeanInnerPotentialVolts = RMeanInnerCrystalPotential*(((RPlanckConstant**2)/ &
       (TWO*RElectronMass*RElectronCharge*TWOPI**2))*&
       RAngstromConversion*RAngstromConversion)

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"StructureFactorInitialisation(",my_rank,") RMeanInnerCrystalPotential = ", &
          RMeanInnerCrystalPotential,RMeanInnerPotentialVolts
  END IF
  
  DO ind=1,nReflections
     CUgMat(ind,ind)=CUgMat(ind,ind)-RMeanInnerCrystalPotential
  ENDDO

  CUgMat = CUgMat + CONJG(TRANSPOSE(CUgMat))

  !Now initialisation calls the Ug calculation subroutines
  !Used to be in Setup

  !--------------------------------------------------------------------
  ! high-energy approximation (not HOLZ compatible)
  !--------------------------------------------------------------------
  
  RBigK= SQRT(RElectronWaveVectorMagnitude**2 + RMeanInnerCrystalPotential)

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"StructureFactorSetup(", my_rank, ") BigK=", RBigK
  END IF
    
  IF(IAbsorbFLAG.GT.1) THEN
     CZeroMAT = CZERO
     
     DO ind = 1,nReflections
        DO jnd = 1,ind
           CZeroMAT(ind,jnd) = CONE
        END DO
     END DO
     
     ALLOCATE( &  
          ISymmetryRelations(nReflections,nReflections), &
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        !refinemain was here
        PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
        !call error function
        RETURN
     ENDIF
     
     ISymmetryRelations = ISymmetryRelations*CZeroMat 

     !Symmetry related structure factor calculation
     !---------------------------------------------
     CALL SymmetryRelatedStructureFactorDetermination (IErr)
     IF( IErr.NE.0 ) THEN
        !refinemain was here
        PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
             " in DetermineSymmetryRelatedUgs"
        !call error function
        RETURN
     ENDIF


     ALLOCATE( &  
          RUniqueUgPrimeValues((SIZE(ISymmetryStrengthKey,DIM=1))), &
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        !refinemain was here
        PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
        !call error function
        RETURN
     ENDIF
     
     DO ind = 1,(SIZE(ISymmetryStrengthKey,DIM=1))
        ILoc = MINLOC(ABS(ISymmetryRelations-ind))
        ISymmetryStrengthKey(ind,:) = ILoc
     END DO
  END IF
  
  IF(IAbsorbFLAG.GE.1) THEN

     !Structure Factors when taking into account absorption
     !Flag controlled
     !-----------------------------------------------------
     CALL StructureFactorsWithAbsorptionDetermination(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"StructureFactorSetup(", my_rank, ") error ", IErr, &
             " in StructureFactorsWithAbsorptionDetermination"
        !call error function
        RETURN
     ENDIF


     IF(IAbsorbFLAG.GE.2) THEN
        DO ind = 2,(SIZE(ISymmetryStrengthKey,DIM=1))
           WHERE (ISymmetryRelations.EQ.ind)
              CUgMatPrime = RUniqueUgPrimeValues(ind)*CIMAGONE
           END WHERE
        END DO
        DO ind = 1,nReflections
           CUgMatPrime(ind,ind) = RUniqueUgPrimeValues(1)*CIMAGONE
        END DO
     END IF
     CUgMat =  CUgMat+CUgMatPrime
  
  ENDIF
  

END SUBROUTINE StructureFactorInitialisation

SUBROUTINE StructureFactorsWithAbsorptionDetermination(IErr)         


  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 
  
  INTEGER(IKIND) IErr,ind,jnd,knd
  REAL(RKIND) :: &
       RIntegrationParameterGMagPrime,RAtomicFormFactorGMagPrime,&
       RAtomicFormFactorGMagMinusGMagPrime,RAbsorpativeAtomicFormFactor,&
       RAbsorpativeAtomicFormFactorInterval
!!$       RAbsorpativeAtomicFormFactorUpperInterval,&
!!$       RAbsorpativeAtomicFormFactorMiddleInterval
  INTEGER(IKIND),DIMENSION(2) :: &
       IPos
  COMPLEX(CKIND) CVgij

  IF((my_rank.EQ.0.AND.IWriteFLAG.GE.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"StructureFactorsWithAbsorptionDetermination(",my_rank,")"
  END IF

  CUgMatPrime = CZERO
    
  SELECT CASE (IAbsorbFLAG)

  CASE(1)

!!$     THE PROPORTIONAL MODEL OF ABSORPTION
     
     CUgMatPrime = CUgMatPrime+(REAL(CUgMat)*(RAbsorptionPercentage/100_RKIND)*CIMAGONE)
     
     DO ind = 1,SIZE(CUgMat,DIM=1)
        CUgMatPrime(ind,ind) = REAL(RMeanInnerCrystalPotential)*(RAbsorptionPercentage/100_RKIND)*CIMAGONE
     END DO
    
  CASE Default 
     RUniqueUgPrimeValues = ZERO

!!$     The Einstein TDS model 
!!$     Equations come from Bird and King 1990 Acta Cryst A46, 202-208 

!!$     DO ind=1,nReflections
!!$        DO jnd=1,nReflections
!!$     DO knd = 1,(SIZE(ISymmetryStrengthKey,DIM=1))
     DO knd = 1,8
        ind = ISymmetryStrengthKey(knd,1)
        jnd = ISymmetryStrengthKey(knd,2)
           !PRINT*,"knd =",knd
           CVgij= CZERO
           
           RGVector = RgMatMat(ind,jnd,:)
           RGVectorMagnitude = RgMatMag(ind,jnd)/(FOUR*PI)

           DO IAtom=1, INAtomsUnitCell
              ICurrentAtom = IAtoms(IAtom)
              
              RAbsorpativeAtomicFormFactor = ZERO
              RAbsorpativeAtomicFormFactorInterval = ZERO

!!$              The integral has to be split up into two between 0 and |G|/2 and 
!!$              |G|/2 and (a large reciprocal distance). This means the integral has to be 
!!$              split into four seperate ones as its a double integral 
              
              ROuterIntegralLowerBound = 0.0D0
              ROuterIntegralUpperBound = RGVectorMagnitude/2.0D0
              RInnerIntegralLowerBound = 0.0D0
              RInnerIntegralUpperBound = RGVectorMagnitude/2.0D0
              
              CALL RIntegrateForAbsorption(RAbsorpativeAtomicFormFactorInterval,IErr)
              
              RAbsorpativeAtomicFormFactor = RAbsorpativeAtomicFormFactor +&
                   RAbsorpativeAtomicFormFactorInterval
              RAbsorpativeAtomicFormFactorInterval = ZERO
              
              ROuterIntegralLowerBound = RGVectorMagnitude/2.0D0
              ROuterIntegralUpperBound = 30.0D0
              RInnerIntegralLowerBound = RGVectorMagnitude/2.0D0
              RInnerIntegralUpperBound = 30.0D0
              CALL RIntegrateForAbsorption(RAbsorpativeAtomicFormFactorInterval,IErr)

              RAbsorpativeAtomicFormFactor = RAbsorpativeAtomicFormFactor +&
                   RAbsorpativeAtomicFormFactorInterval
              RAbsorpativeAtomicFormFactorInterval = ZERO

              ROuterIntegralLowerBound = 0.0D0
              ROuterIntegralUpperBound = RGVectorMagnitude/2.0D0
              RInnerIntegralLowerBound = RGVectorMagnitude/2.0D0
              RInnerIntegralUpperBound = 30.0D0

              CALL RIntegrateForAbsorption(RAbsorpativeAtomicFormFactorInterval,IErr)

              RAbsorpativeAtomicFormFactor = RAbsorpativeAtomicFormFactor +&
                   RAbsorpativeAtomicFormFactorInterval
              RAbsorpativeAtomicFormFactorInterval = ZERO

              ROuterIntegralLowerBound = RGVectorMagnitude/2.0D0
              ROuterIntegralUpperBound = 30.0D0 
              RInnerIntegralLowerBound = 0.0D0
              RInnerIntegralUpperBound = RGVectorMagnitude/2.0D0
              CALL RIntegrateForAbsorption(RAbsorpativeAtomicFormFactorInterval,IErr)

              RAbsorpativeAtomicFormFactor = RAbsorpativeAtomicFormFactor +&
                   RAbsorpativeAtomicFormFactorInterval
              RAbsorpativeAtomicFormFactorInterval = ZERO
              
              RAbsorpativeAtomicFormFactor=&
                   RAbsorpativeAtomicFormFactor*ROcc(IAtom)

              CVgij = CVgij + &
                   RAbsorpativeAtomicFormFactor * &
                   EXP(-CIMAGONE* &
                   DOT_PRODUCT(RgMatMat(ind,jnd,:), RrVecMat(iAtom,:)) &
                   )*EXP(-RDWF(IAtom)*RGVectorMagnitude/2.0D0)
           ENDDO

           

           PRINT*,CVgij
           
!!$           Calculate Vg'

           RUniqueUgPrimeValues(knd) = REAL((((RPlanckConstant**3)/&
                (PI*RElectronMass**2*((RElectronVelocity)/RSpeedOfLight)*&
                RSpeedofLight * RVolume))*RAngstromConversion**3)*CVgij,RKIND)
           PRINT*,RUniqueUgPrimeValues(knd)

!!$           Convert Vg' into Ug'

           RUniqueUgPrimeValues(knd) = REAL(RUniqueUgPrimeValues(knd)*&
                (TWO*TWOPI**2*RElectronMass*RElectronCharge)/&
                (RPlanckConstant**2))

!!$           Correct for relativistic effects

           RUniqueUgPrimeValues(knd) = REAL(RUniqueUgPrimeValues(knd)*&
                RRelativisticCorrection,RKIND)

        ENDDO
!!$     ENDDO
     
  END SELECT

  IErr = 0
  
END SUBROUTINE StructureFactorsWithAbsorptionDetermination

REAL(RKIND) FUNCTION RAbsorpativeIntegrand(RIntegrationParameterGMagPrime)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara
  USE IChannels
 
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 

  REAL(RKIND) :: &
       RAtomicFormFactorGMagPrime,&
       RAtomicFormFactorGMagMinusGMagPrime,&
!!       RAbsorpativeIntegrand,&
       RIntegrationParameterGMagPrime,&
       OneDIntegral
  INTEGER(IKIND) ierr,currentatom
  INTEGER(IKIND),DIMENSION(2) :: &
       IPos
  COMPLEX(CKIND) CVgij

!!$  RAbsorpativeIntegrand = SIN(ROuterIntegrationParameterGMagPrime)*SIN(RIntegrationParameterGMagPrime)

  RAbsorpativeIntegrand = OneDIntegral(ROuterIntegrationParameterGMagPrime)*OneDIntegral(RIntegrationParameterGMagPrime)

END FUNCTION RAbsorpativeIntegrand

REAL(RKIND) FUNCTION RGXIntegration(RIntegrationParameterGMagPrime)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 

  REAL(RKIND) :: &
       RIntegrationParameterGMagPrime
  REAL(RKIND) :: &
!!       RGXIntegration,&
       RAbsorpativeIntegrand,RAbsoluteError
  INTEGER(IKIND) :: &
       IErr,IIntegrationSteps
   
  EXTERNAL RAbsorpativeIntegrand

  ROuterIntegrationParameterGMagPrime = RIntegrationParameterGMagPrime

  CALL DQNG(RAbsorpativeIntegrand,RInnerIntegralLowerBound,RInnerIntegralUpperBound,&
       0.0D0,1.0D-3,RGXIntegration,RAbsoluteError,IIntegrationSteps,IErr)

END FUNCTION RGXIntegration

SUBROUTINE RIntegrateForAbsorption(RAbsorpativeAtomicFormFactor,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 

  INTEGER(IKIND) :: &
       IErr,IIntegrationSteps
  REAL(RKIND) RAbsorpativeAtomicFormFactor,&
       RGXIntegration,RAbsoluteError
      
  EXTERNAL RGXIntegration

  CALL DQNG(RGXIntegration,ROuterIntegralLowerBound,ROuterIntegralUpperBound,&
       0.0D0,1.0D-3,RAbsorpativeAtomicFormFactor,RAbsoluteError,IIntegrationSteps,IErr)

END SUBROUTINE RIntegrateForAbsorption

REAL(RKIND) FUNCTION  OneDIntegral(X)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels
 
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 
  
  REAL(RKIND) :: &
       X

  REAL(RKIND) :: &
       RAtomicFormFactorGMagPrime,&
       RAtomicFormFactorGMagMinusGMagPrime,&
       RAbsorpativeIntegrand,&
       RIntegrationParameterGMagPrime,&
       Gaussian,Lorentzian
!!    OneDIntegral,
            

  INTEGER knd,evenindlorentz,oddindlorentz,evenindgauss,oddindgauss
 
  SELECT CASE (IScatterFactorMethodFLAG)
     
  CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians 

     RAtomicFormFactorGMagPrime = ZERO
     RAtomicFormFactorGMagMinusGMagPrime = ZERO
     DO knd = 1, 3
              
        !odd and even indicies for Lorentzian function
        evenindlorentz = knd*2
        oddindlorentz = knd*2 -1
                 
        !odd and even indicies for Gaussian function
        evenindgauss = evenindlorentz + 6
        oddindgauss = oddindlorentz + 6
                 
        !Kirkland Method uses summation of 3 Gaussians and 3 Lorentzians (summed in loop)
        !This may not work... have found no way to test it.
        RAtomicFormFactorGMagPrime = RAtomicFormFactorGMagPrime + &
             !3 Lorentzians
             LORENTZIAN(RScattFactors(ICurrentAtom,oddindlorentz), RGVectorMagnitude/2.0D0,X,&
             RScattFactors(ICurrentAtom,evenindlorentz))+ &
             !3 Gaussians
             GAUSSIAN(RScattFactors(ICurrentAtom,oddindgauss),(RgVectorMagnitude/2.0D0)+X,ZERO, & 
             1/(SQRT(2*RScattFactors(ICurrentAtom,evenindgauss))),ZERO)
              

    ! RAtomicFormFactorGMagPrime = &
          ! 3 Lorentzians
    !      RScattFactors(ICurrentAtom,1) / &
     !     (((RGVectorMagnitude/2.0D0)+X)**2 + RScattFactors(ICurrentAtom,2)) + &
     !     RScattFactors(ICurrentAtom,3) / &
     !     (((RGVectorMagnitude/2.0D0)+X)**2 + RScattFactors(ICurrentAtom,4)) + &
     !     RScattFactors(ICurrentAtom,5) / &
     !     (((RGVectorMagnitude/2.0D0)+X)**2 + RScattFactors(ICurrentAtom,6)) + &
          ! 3 Gaussians
     !     RScattFactors(ICurrentAtom,7) * &
     !     EXP(-(((RGVectorMagnitude/2.0D0)+X)**2)* RScattFactors(ICurrentAtom,8)) + &
     !     RScattFactors(ICurrentAtom,9) * &
     !     EXP(-(((RGVectorMagnitude/2.0D0)+X)**2)* RScattFactors(ICurrentAtom,10)) + &
     !     RScattFactors(ICurrentAtom,11) * &
     !     EXP(-(((RGVectorMagnitude/2.0D0)+X)**2)* RScattFactors(ICurrentAtom,12))
       
        !This may not work -  same reason as above, if problems arise see util.f90 with 
        !gaussian and lorentzian functions
        RAtomicFormFactorGMagMinusGMagPrime = RAtomicFormFactorGMagMinusGMagPrime + &
             !3 Lorentzians
             LORENTZIAN(RScattFactors(ICurrentAtom,oddindlorentz), RGVectorMagnitude/2.0D0,-X,&
             RScattFactors(ICurrentAtom,evenindlorentz))+ &
             !3 Gaussians
             GAUSSIAN(RScattFactors(ICurrentAtom,oddindgauss),(RgVectorMagnitude/2.0D0)-X,ZERO, & 
             1/(SQRT(2*RScattFactors(ICurrentAtom,evenindgauss))),ZERO)

     END DO


    ! RAtomicFormFactorGMagMinusGMagPrime = &
          ! 3 Lorentzians
     !     RScattFactors(ICurrentAtom,1) / &
     !     (((RGVectorMagnitude/2.0D0)-X)**2 + RScattFactors(ICurrentAtom,2)) + &
     !     RScattFactors(ICurrentAtom,3) / &
     !     (((RGVectorMagnitude/2.0D0)-X)**2 + RScattFactors(ICurrentAtom,4)) + &
     !     RScattFactors(ICurrentAtom,5) / &
     !     (((RGVectorMagnitude/2.0D0)-X)**2 + RScattFactors(ICurrentAtom,6)) + &
          ! 3 Gaussians
     !     RScattFactors(ICurrentAtom,7) * &
     !     EXP(-((RGVectorMagnitude/2.0D0)-X)**2 * RScattFactors(ICurrentAtom,8)) + &
     !     RScattFactors(ICurrentAtom,9) * &
     !     EXP(-((RGVectorMagnitude/2.0D0)-X)**2 * RScattFactors(ICurrentAtom,10)) + &
     !     RScattFactors(ICurrentAtom,11) * &
     !     EXP(-((RGVectorMagnitude/2.0D0)-X)**2 * RScattFactors(ICurrentAtom,12))
     
  CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
     RAtomicFormFactorGMagPrime = &
          RScattFactors(ICurrentAtom,1) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,5)) + &
          RScattFactors(ICurrentAtom,2) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,6)) + &
          RScattFactors(ICurrentAtom,3) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,7)) + &
          RScattFactors(ICurrentAtom,4) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,8))
     
     RAtomicFormFactorGMagMinusGMagPrime = &
          RScattFactors(ICurrentAtom,1) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,5)) + &
          RScattFactors(ICurrentAtom,2) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
          RScattFactors(ICurrentAtom,3) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,7)) + &
          RScattFactors(ICurrentAtom,4) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,8))
     
  CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)
     
     RAtomicFormFactorGMagPrime = &
          RScattFactors(ICurrentAtom,1) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,2)) + &
          RScattFactors(ICurrentAtom,3) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,4)) + &
          RScattFactors(ICurrentAtom,5) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,6)) + &
          RScattFactors(ICurrentAtom,7) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,8))
     
     RAtomicFormFactorGMagMinusGMagPrime = &
          RScattFactors(ICurrentAtom,1) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,2)) + &
          RScattFactors(ICurrentAtom,3) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,4)) + &
          RScattFactors(ICurrentAtom,5) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
          RScattFactors(ICurrentAtom,7) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,8))
     
  END SELECT
  
  RAbsorpativeIntegrand = RAtomicFormFactorGMagPrime*RAtomicFormFactorGMagMinusGMagPrime
  
!!$  RAbsorpativeIntegrand = &
!!$       RAbsorpativeIntegrand*ROcc(IAtom)
  
  ! initialize potential as in Eq. (6.10) of Kirkland
  
  
  IF (IAnisoDebyeWallerFactorFlag.EQ.0) THEN
     
     IF(RDWF(IAtom).GT.10.OR.RDWF(IAtom).LT.0) THEN
        RDWF(IAtom) = RDebyeWallerConstant
     END IF
     
     RAbsorpativeIntegrand = RAbsorpativeIntegrand * &
          (1.0D0-EXP(-RDWF(IAtom)*(X**2-(RGVectorMagnitude/4.0D0))))
  ELSE
     
     RAbsorpativeIntegrand = RAbsorpativeIntegrand * &
          EXP(-TWOPI*DOT_PRODUCT(RGVector, &
          MATMUL( RAnisotropicDebyeWallerFactorTensor( &
          IAnisoDWFT(IAtom),:,:), &
          RGVector)))
     
     
  END IF
  
  OneDIntegral = RAbsorpativeIntegrand

END FUNCTION OneDIntegral



!!$REAL FUNCTION  OneDIntegral(X)
!!$  
!!$  USE MyNumbers
!!$  
!!$  USE CConst; USE IConst
!!$  USE IPara; USE RPara; USE CPara
!!$  USE BlochPara
!!$  USE IChannels
!!$  USE MPI
!!$  USE MyMPI
!!$  
!!$  IMPLICIT NONE 
!!$  
!!$  REAL(RKIND) :: &
!!$       X
!!$
!!$  REAL(RKIND) :: &
!!$       RAtomicFormFactorGMagPrime,&
!!$       RAtomicFormFactorGMagMinusGMagPrime,&
!!$       RAbsorpativeIntegrand,&
!!$       RIntegrationParameterGMagPrime,&
!!$       OneDIntegral
!!$       
!!$  
!!$  SELECT CASE (IScatterFactorMethodFLAG)
!!$     
!!$  CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians 
!!$     
!!$     RAtomicFormFactorGMagPrime = &
!!$          ! 3 Lorentzians
!!$          RScattFactors(ICurrentAtom,1) / &
!!$          (X**2 + RScattFactors(ICurrentAtom,2)) + &
!!$          RScattFactors(ICurrentAtom,3) / &
!!$          (X**2 + RScattFactors(ICurrentAtom,4)) + &
!!$          RScattFactors(ICurrentAtom,5) / &
!!$          (X**2 + RScattFactors(ICurrentAtom,6)) + &
!!$          ! 3 Gaussians
!!$          RScattFactors(ICurrentAtom,7) * &
!!$          EXP(-(X**2)* RScattFactors(ICurrentAtom,8)) + &
!!$          RScattFactors(ICurrentAtom,9) * &
!!$          EXP(-(X**2)* RScattFactors(ICurrentAtom,10)) + &
!!$          RScattFactors(ICurrentAtom,11) * &
!!$          EXP(-(X**2)* RScattFactors(ICurrentAtom,12))
!!$     
!!$     RAtomicFormFactorGMagMinusGMagPrime = &
!!$          ! 3 Lorentzians
!!$          RScattFactors(ICurrentAtom,1) / &
!!$          ((RGVectorMagnitude-X)**2 + RScattFactors(ICurrentAtom,2)) + &
!!$          RScattFactors(ICurrentAtom,3) / &
!!$          ((RGVectorMagnitude-X)**2 + RScattFactors(ICurrentAtom,4)) + &
!!$          RScattFactors(ICurrentAtom,5) / &
!!$          ((RGVectorMagnitude-X)**2 + RScattFactors(ICurrentAtom,6)) + &
!!$          ! 3 Gaussians
!!$          RScattFactors(ICurrentAtom,7) * &
!!$          EXP(-(RGVectorMagnitude-X)**2 * RScattFactors(ICurrentAtom,8)) + &
!!$          RScattFactors(ICurrentAtom,9) * &
!!$          EXP(-(RGVectorMagnitude-X)**2 * RScattFactors(ICurrentAtom,10)) + &
!!$          RScattFactors(ICurrentAtom,11) * &
!!$          EXP(-(RGVectorMagnitude-X)**2 * RScattFactors(ICurrentAtom,12))
!!$     
!!$  CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
!!$     RAtomicFormFactorGMagPrime = &
!!$          RScattFactors(ICurrentAtom,1) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,5)) + &
!!$          RScattFactors(ICurrentAtom,2) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!!$          RScattFactors(ICurrentAtom,3) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,7)) + &
!!$          RScattFactors(ICurrentAtom,4) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,8))
!!$     
!!$     RAtomicFormFactorGMagMinusGMagPrime = &
!!$          RScattFactors(ICurrentAtom,1) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,5)) + &
!!$          RScattFactors(ICurrentAtom,2) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!!$          RScattFactors(ICurrentAtom,3) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,7)) + &
!!$          RScattFactors(ICurrentAtom,4) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,8))
!!$     
!!$  CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)
!!$     
!!$     RAtomicFormFactorGMagPrime = &
!!$          RScattFactors(ICurrentAtom,1) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,2)) + &
!!$          RScattFactors(ICurrentAtom,3) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,4)) + &
!!$          RScattFactors(ICurrentAtom,5) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!!$          RScattFactors(ICurrentAtom,7) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,8))
!!$     
!!$     RAtomicFormFactorGMagMinusGMagPrime = &
!!$          RScattFactors(ICurrentAtom,1) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,2)) + &
!!$          RScattFactors(ICurrentAtom,3) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,4)) + &
!!$          RScattFactors(ICurrentAtom,5) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!!$          RScattFactors(ICurrentAtom,7) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,8))
!!$     
!!$  END SELECT
!!$  
!!$  RAbsorpativeIntegrand = RAtomicFormFactorGMagPrime*RAtomicFormFactorGMagMinusGMagPrime
!!$  
!!$  RAbsorpativeIntegrand = &
!!$       RAbsorpativeIntegrand*ROcc(IAtom)
!!$  
!!$  ! initialize potential as in Eq. (6.10) of Kirkland
!!$  
!!$  
!!$  IF (IAnisoDebyeWallerFactorFlag.EQ.0) THEN
!!$     
!!$     IF(RDWF(IAtom).GT.10.OR.RDWF(IAtom).LT.0) THEN
!!$        RDWF(IAtom) = RDebyeWallerConstant
!!$     END IF
!!$     
!!$     RAbsorpativeIntegrand = RAbsorpativeIntegrand * &
!!$          ( &
!!$          EXP(-((RGVectorMagnitude/2.D0)**2)*RDWF(IAtom))-&
!!$          EXP(-((X/2.D0)**2)*RDWF(IAtom))*&
!!$          EXP(-(((RGVectorMagnitude-X)/2.D0)**2)*RDWF(IAtom)))
!!$  ELSE
!!$     
!!$     RAbsorpativeIntegrand = RAbsorpativeIntegrand * &
!!$          EXP(-TWOPI*DOT_PRODUCT(RGVector, &
!!$          MATMUL( RAnisotropicDebyeWallerFactorTensor( &
!!$          IAnisoDWFT(IAtom),:,:), &
!!$          RGVector)))
!!$     
!!$     
!!$  END IF
!!$  
!!$  OneDIntegral = RAbsorpativeIntegrand
!!$
!!$END FUNCTION OneDIntegral
