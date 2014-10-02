!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! FelixSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
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
!  This file is part of FelixSim.
!
!  FelixSim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  FelixSim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with FelixSim.  If not, see <http://www.gnu.org/licenses/>.
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
     PRINT*,"GMatrixInitialisation()"
  END IF

  DO ind=1,nReflections
     DO jnd=1,nReflections
        RgMatMat(ind,jnd,:)= RgVecMatT(ind,:)-RgVecMatT(jnd,:)
        RgMatMag(ind,jnd)= SQRT(DOT_PRODUCT(RgMatMat(ind,jnd,:),RgMatMat(ind,jnd,:)))
     ENDDO
  ENDDO
   
  RgMatMag = RgMatMag/TWOPI
  
END SUBROUTINE GMatrixInitialisation

SUBROUTINE DetermineSymmetryRelatedUgs (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI
  USE CPara
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,ierr,knd,Iuid


  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DetermineSymmetryRelatedUgs(",my_rank,")"
  END IF

  !Immediately set all the zeros to Relation 1
  
  ISymmetryRelations = 0
  
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
  
END SUBROUTINE DetermineSymmetryRelatedUgs

!---------------------------------------------------------------------
SUBROUTINE UgCalculation (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,ierr,imaxj,IFound,ICount,currentatom
  INTEGER(IKIND),DIMENSION(2) :: &
       IPos
  COMPLEX(CKIND) CVgij
  REAL(RKIND) :: &
       RMeanInnerPotentialVolts,RAtomicFormFactor  

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"UgCalculation()"
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
              
              RAtomicFormFactor = &
                   ! 3 Lorentzians
                   RScattFactors(ICurrentAtom,1) / &
                   (RgMatMag(ind,jnd)**2 + RScattFactors(ICurrentAtom,2)) + &
                   RScattFactors(ICurrentAtom,3) / &
                   (RgMatMag(ind,jnd)**2 + RScattFactors(ICurrentAtom,4)) + &
                   RScattFactors(ICurrentAtom,5) / &
                   (RgMatMag(ind,jnd)**2 + RScattFactors(ICurrentAtom,6)) + &
                   ! 3 Gaussians
                   RScattFactors(ICurrentAtom,7) * &
                   EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(ICurrentAtom,8)) + &
                   RScattFactors(ICurrentAtom,9) * &
                   EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(ICurrentAtom,10)) + &
                   RScattFactors(ICurrentAtom,11) * &
                   EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(ICurrentAtom,12))
              
           CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
              RAtomicFormFactor = &
                   RScattFactors(ICurrentAtom,1) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,5)) + &
                   RScattFactors(ICurrentAtom,2) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
                   RScattFactors(ICurrentAtom,3) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,7)) + &
                   RScattFactors(ICurrentAtom,4) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,8))

           CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)

              RAtomicFormFactor = &
                   RScattFactors(ICurrentAtom,1) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,2)) + &
                   RScattFactors(ICurrentAtom,3) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,4)) + &
                   RScattFactors(ICurrentAtom,5) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
                   RScattFactors(ICurrentAtom,7) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,8))

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
     PRINT*,"UgCalculation(",my_rank,") RMeanInnerCrystalPotential = ",RMeanInnerCrystalPotential,RMeanInnerPotentialVolts
  END IF
  
  DO ind=1,nReflections
     CUgMat(ind,ind)=CUgMat(ind,ind)-RMeanInnerCrystalPotential
  ENDDO

  CUgMat = CUgMat + CONJG(TRANSPOSE(CUgMat))

END SUBROUTINE UgCalculation

SUBROUTINE UgAddAbsorption(IErr)         


  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 
  
  INTEGER(IKIND) IErr,ind,jnd
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
     PRINT*,"UgAddAbsorption(",my_rank,")"
  END IF

  CUgMatPrime = CZERO
  
  PRINT*,INAtomsUnitCell
  
  SELECT CASE (IAbsorbFLAG)

  CASE(1)

     !THE PROPORTIONAL MODEL OF ABSORPTION
     
     CUgMatPrime = CUgMatPrime+(REAL(CUgMat)*(RAbsorptionPercentage/100_RKIND)*CIMAGONE)
     
     DO ind = 1,SIZE(CUgMat,DIM=1)
        CUgMatPrime(ind,ind) = REAL(RBigK)*(RAbsorptionPercentage/100_RKIND)*CIMAGONE
     END DO
     
  CASE(2)
     
     DO ind=1,nReflections
        DO jnd=1,nReflections

           
           CVgij= 0.0D0
           
           RGVector = RgMatMat(ind,jnd,:)
           RGVectorMagnitude = RgMatMag(ind,jnd)

           DO IAtom=1, INAtomsUnitCell
              ICurrentAtom = IAtoms(IAtom)
              
              RAbsorpativeAtomicFormFactor = ZERO
              RAbsorpativeAtomicFormFactorInterval = ZERO
              ! NOW INTEGRATE RAbsorpativeAtomicFormFactor OVER RIntegrationParameterGMagPrime FROM 0 TO 30ANGSTROMS

              RIntegralLowerBound = 0.0D0
              RIntegralUpperBound = 200.0D0

              CALL RIntegrateForAbsorption(RAbsorpativeAtomicFormFactorInterval,IErr)

              RAbsorpativeAtomicFormFactor = RAbsorpativeAtomicFormFactor +&
                   RAbsorpativeAtomicFormFactorInterval
              RAbsorpativeAtomicFormFactorInterval = ZERO

              CVgij = CVgij + &
                   RAbsorpativeAtomicFormFactor * &
                   EXP(-CIMAGONE* &
                   DOT_PRODUCT(RgMatMat(ind,jnd,:), RrVecMat(iAtom,:)) &
                   )
           ENDDO
           
           CUgMatPrime(ind,jnd) = ((((TWO*RPlanckConstant**3*TWOPI**2)/&
                (RElectronMass**2*(RElectronVelocity/RSpeedOfLight)*&
                RSpeedofLight * RVolume))*RAngstromConversion**3)/(((RPlanckConstant**2)/ &
                (TWO*RElectronMass*RElectronCharge*TWOPI**2))&
                ))*CVgij*CIMAGONE
        ENDDO
     ENDDO
     
  END SELECT
  
END SUBROUTINE UgAddAbsorption

REAL FUNCTION RAbsorpativeIntegrand(RIntegrationParameterGMagPrime)

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
       RAbsorpativeIntegrand,&
       RIntegrationParameterGMagPrime
  INTEGER(IKIND) ierr,currentatom
  INTEGER(IKIND),DIMENSION(2) :: &
       IPos
  COMPLEX(CKIND) CVgij

!!$  SELECT CASE (IScatterFactorMethodFLAG)
!!$     
!!$  CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians 
!!$     
!!$     RAtomicFormFactorGMagPrime = &
!!$          ! 3 Lorentzians
!!$          RScattFactors(ICurrentAtom,1) / &
!!$          (RIntegrationParameterGMagPrime**2 + RScattFactors(ICurrentAtom,2)) + &
!!$          RScattFactors(ICurrentAtom,3) / &
!!$          (RIntegrationParameterGMagPrime**2 + RScattFactors(ICurrentAtom,4)) + &
!!$          RScattFactors(ICurrentAtom,5) / &
!!$          (RIntegrationParameterGMagPrime**2 + RScattFactors(ICurrentAtom,6)) + &
!!$          ! 3 Gaussians
!!$          RScattFactors(ICurrentAtom,7) * &
!!$          EXP(-RIntegrationParameterGMagPrime**2 * RScattFactors(ICurrentAtom,8)) + &
!!$          RScattFactors(ICurrentAtom,9) * &
!!$          EXP(-RIntegrationParameterGMagPrime**2 * RScattFactors(ICurrentAtom,10)) + &
!!$          RScattFactors(ICurrentAtom,11) * &
!!$          EXP(-RIntegrationParameterGMagPrime**2 * RScattFactors(ICurrentAtom,12))
!!$     
!!$     RAtomicFormFactorGMagMinusGMagPrime = &
!!$          ! 3 Lorentzians
!!$          RScattFactors(ICurrentAtom,1) / &
!!$          ((RGVectorMagnitude-RIntegrationParameterGMagPrime)**2 + RScattFactors(ICurrentAtom,2)) + &
!!$          RScattFactors(ICurrentAtom,3) / &
!!$          ((RGVectorMagnitude-RIntegrationParameterGMagPrime)**2 + RScattFactors(ICurrentAtom,4)) + &
!!$          RScattFactors(ICurrentAtom,5) / &
!!$          ((RGVectorMagnitude-RIntegrationParameterGMagPrime)**2 + RScattFactors(ICurrentAtom,6)) + &
!!$          ! 3 Gaussians
!!$          RScattFactors(ICurrentAtom,7) * &
!!$          EXP(-(RGVectorMagnitude-RIntegrationParameterGMagPrime)**2 * RScattFactors(ICurrentAtom,8)) + &
!!$          RScattFactors(ICurrentAtom,9) * &
!!$          EXP(-(RGVectorMagnitude-RIntegrationParameterGMagPrime)**2 * RScattFactors(ICurrentAtom,10)) + &
!!$          RScattFactors(ICurrentAtom,11) * &
!!$          EXP(-(RGVectorMagnitude-RIntegrationParameterGMagPrime)**2 * RScattFactors(ICurrentAtom,12))
!!$     
!!$  CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
!!$     RAtomicFormFactorGMagPrime = &
!!$          RScattFactors(ICurrentAtom,1) * &
!!$          EXP(-(RIntegrationParameterGMagPrime**2)/4 * RScattFactors(ICurrentAtom,5)) + &
!!$          RScattFactors(ICurrentAtom,2) * &
!!$          EXP(-(RIntegrationParameterGMagPrime**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!!$          RScattFactors(ICurrentAtom,3) * &
!!$          EXP(-(RIntegrationParameterGMagPrime**2)/4 * RScattFactors(ICurrentAtom,7)) + &
!!$          RScattFactors(ICurrentAtom,4) * &
!!$          EXP(-(RIntegrationParameterGMagPrime**2)/4 * RScattFactors(ICurrentAtom,8))
!!$     
!!$     RAtomicFormFactorGMagMinusGMagPrime = &
!!$          RScattFactors(ICurrentAtom,1) * &
!!$          EXP(-((RGVectorMagnitude-RIntegrationParameterGMagPrime)**2)/4 * RScattFactors(ICurrentAtom,5)) + &
!!$          RScattFactors(ICurrentAtom,2) * &
!!$          EXP(-((RGVectorMagnitude-RIntegrationParameterGMagPrime)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!!$          RScattFactors(ICurrentAtom,3) * &
!!$          EXP(-((RGVectorMagnitude-RIntegrationParameterGMagPrime)**2)/4 * RScattFactors(ICurrentAtom,7)) + &
!!$          RScattFactors(ICurrentAtom,4) * &
!!$          EXP(-((RGVectorMagnitude-RIntegrationParameterGMagPrime)**2)/4 * RScattFactors(ICurrentAtom,8))
!!$     
!!$  CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)
!!$     
!!$     RAtomicFormFactorGMagPrime = &
!!$          RScattFactors(ICurrentAtom,1) * &
!!$          EXP(-(RIntegrationParameterGMagPrime**2)/4 * RScattFactors(ICurrentAtom,2)) + &
!!$          RScattFactors(ICurrentAtom,3) * &
!!$          EXP(-(RIntegrationParameterGMagPrime**2)/4 * RScattFactors(ICurrentAtom,4)) + &
!!$          RScattFactors(ICurrentAtom,5) * &
!!$          EXP(-(RIntegrationParameterGMagPrime**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!!$          RScattFactors(ICurrentAtom,7) * &
!!$          EXP(-(RIntegrationParameterGMagPrime**2)/4 * RScattFactors(ICurrentAtom,8))
!!$     
!!$     RAtomicFormFactorGMagMinusGMagPrime = &
!!$          RScattFactors(ICurrentAtom,1) * &
!!$          EXP(-((RGVectorMagnitude-RIntegrationParameterGMagPrime)**2)/4 * RScattFactors(ICurrentAtom,2)) + &
!!$          RScattFactors(ICurrentAtom,3) * &
!!$          EXP(-((RGVectorMagnitude-RIntegrationParameterGMagPrime)**2)/4 * RScattFactors(ICurrentAtom,4)) + &
!!$          RScattFactors(ICurrentAtom,5) * &
!!$          EXP(-((RGVectorMagnitude-RIntegrationParameterGMagPrime)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!!$          RScattFactors(ICurrentAtom,7) * &
!!$          EXP(-((RGVectorMagnitude-RIntegrationParameterGMagPrime)**2)/4 * RScattFactors(ICurrentAtom,8))
!!$     
!!$  END SELECT
!!$  
!!$  RAbsorpativeIntegrand = RAtomicFormFactorGMagPrime*RAtomicFormFactorGMagMinusGMagPrime
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
!!$          EXP(-((RIntegrationParameterGMagPrime/2.D0)**2)*RDWF(IAtom))*&
!!$          EXP(-(((RGVectorMagnitude-RIntegrationParameterGMagPrime)/2.D0)**2)*RDWF(IAtom)))
!!$     
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
  
  RAbsorpativeIntegrand = SIN(RIntegrationParameterGMagPrime)

END FUNCTION RAbsorpativeIntegrand

REAL FUNCTION RGXIntegration(RIntegrationParameterGMagPrime)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 

  REAL(RKIND) :: &
       RIntegrationParameterGMagPrime,RGXIntegration,&
       RAbsorpativeIntegrand,RAbsoluteError
  INTEGER(IKIND) :: &
       IErr,IIntegrationSteps

  EXTERNAL RAbsorpativeIntegrand

  CALL DQNG(RAbsorpativeIntegrand,RIntegralLowerBound,RIntegralUpperBound,&
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

  CALL DQNG(RGXIntegration,RIntegralLowerBound,RIntegralUpperBound,&
       0.0D0,1.0D-3,RAbsorpativeAtomicFormFactor,RAbsoluteError,IIntegrationSteps,IErr)

END SUBROUTINE RIntegrateForAbsorption
