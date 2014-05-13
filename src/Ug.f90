!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! FelixSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
!
! (C) 2013/14, all rights reserved
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
  
  INTEGER ind,jnd,ierr

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

!---------------------------------------------------------------------
SUBROUTINE UgCalculation (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,ierr, currentatom, iAtom
  COMPLEX(CKIND) CVgij
  REAL(RKIND) RAtomicFormFactor
  REAL(RKIND) :: &
       RMeanInnerPotentialVolts

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"UgCalculation()"
  END IF
  
  DO ind=1,nReflections
     DO jnd=1,nReflections
        
        CVgij= 0.0D0
        
        DO iAtom=1, INAtomsUnitCell
           currentatom = IAtoms(iAtom)
           ! calculate f_e(q) as in Eq. (C.15) of Kirkland, "Advanced Computing in EM"
           
           SELECT CASE (IScatterFactorMethodFLAG)
              
           CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians 
              
              RAtomicFormFactor = &
                   ! 3 Lorentzians
                   RScattFactors(currentatom,1) / &
                   (RgMatMag(ind,jnd)**2 + RScattFactors(currentatom,2)) + &
                   RScattFactors(currentatom,3) / &
                   (RgMatMag(ind,jnd)**2 + RScattFactors(currentatom,4)) + &
                   RScattFactors(currentatom,5) / &
                   (RgMatMag(ind,jnd)**2 + RScattFactors(currentatom,6)) + &
                   ! 3 Gaussians
                   RScattFactors(currentatom,7) * &
                   EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(currentatom,8)) + &
                   RScattFactors(currentatom,9) * &
                   EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(currentatom,10)) + &
                   RScattFactors(currentatom,11) * &
                   EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(currentatom,12))
              
           CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
              RAtomicFormFactor = &
                   RScattFactors(currentatom,1) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(currentatom,5)) + &
                   RScattFactors(currentatom,2) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(currentatom,6)) + &
                   RScattFactors(currentatom,3) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(currentatom,7)) + &
                   RScattFactors(currentatom,4) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(currentatom,8))

           CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)

              RAtomicFormFactor = &
                   RScattFactors(currentatom,1) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(currentatom,2)) + &
                   RScattFactors(currentatom,3) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(currentatom,4)) + &
                   RScattFactors(currentatom,5) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(currentatom,6)) + &
                   RScattFactors(currentatom,7) * &
                   EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(currentatom,8))

           END SELECT
              
          ! initialize potential as in Eq. (6.10) of Kirkland

           RAtomicFormFactor = RAtomicFormFactor*ROcc(iAtom)

           IF (IAnisoDebyeWallerFactorFlag.EQ.0) THEN
              
              IF(RDWF(iAtom).GT.1.OR.RDWF(iAtom).LT.0) THEN
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
                
        IF (IAbsorbFlag.EQ.1) THEN
           CUgMat(ind,jnd) = &
                CUgMat(ind,jnd) + &
                ABS(CUgMat(ind,jnd))*(RAbsorptionPercentage/100.D0)*CIMAGONE       
        END IF
        
     ENDDO
  ENDDO

  RMeanInnerCrystalPotential= REAL(CUgMat(1,1))
  RMeanInnerPotentialVolts = ((RMeanInnerCrystalPotential*RPlanckConstant**2)/ &
       (TWO*RElectronMass*RElectronCharge*TWOPI**2))*&
       RAngstromConversion*RAngstromConversion

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"UgCalculation(",my_rank,") RMeanInnerCrystalPotential = ",RMeanInnerCrystalPotential,RMeanInnerPotentialVolts
  END IF
  
  DO ind=1,nReflections
     CUgMat(ind,ind)=CUgMat(ind,ind)-RMeanInnerCrystalPotential
  ENDDO


END SUBROUTINE UgCalculation
