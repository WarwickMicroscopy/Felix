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
  
  INTEGER ind,jnd,ierr,IUniqueKey,knd,IFound
!!$  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
!!$       RUniqueKeyComplete

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"GMatrixInitialisation()"
  END IF

!!$  DO ind = 1,10
!!$     PRINT*
!!$  END DO

!!$  IUniqueKey = 0
  
  DO ind=1,nReflections
     DO jnd=1,nReflections
!!$        IFound = 0
        
        RgMatMat(ind,jnd,:)= RgVecMatT(ind,:)-RgVecMatT(jnd,:)
        RgMatMag(ind,jnd)= SQRT(DOT_PRODUCT(RgMatMat(ind,jnd,:),RgMatMat(ind,jnd,:)))
!!$        IF(IUniqueKey.EQ.0) THEN
!!$           IUniqueKey = 1
!!$           RUniqueKey(IUniqueKey,2:4) = RgMatMat(ind,jnd,:)
!!$           RUniqueKey(IUniqueKey,1) = IUniqueKey
!!$        ELSE
!!$           DO knd = 1,IUniqueKey
!!$              IF(RgMatMat(ind,jnd,1).EQ.RUniqueKey(knd,2).AND.&
!!$                   RgMatMat(ind,jnd,2).EQ.RUniqueKey(knd,3).AND.&
!!$                   RgMatMat(ind,jnd,3).EQ.RUniqueKey(knd,4)) THEN
!!$                 IFound = knd
!!$                 EXIT
!!$              ELSE
!!$                 CYCLE      
!!$              END IF
!!$           END DO
!!$           IF(IFound.EQ.0) THEN
!!$              IUniqueKey = IUniqueKey + 1
!!$              RUniqueKey(IUniqueKey,1) = IUniqueKey
!!$              RUniqueKey(IUniqueKey,2:4) = RgMatMat(ind,jnd,:)
!!$              ISymmetryRelations(ind,jnd) = IUniqueKey
!!$           ELSE
!!$              ISymmetryRelations(ind,jnd) = knd              
!!$           END IF
!!$        END IF
     ENDDO
  ENDDO
 
!!$  ALLOCAte(&
!!$      RUniqueKeyComplete(IUniqueKey,5),&
!!$      STAT = IErr)
!!$  IF( IErr.NE.0 ) THEN
!!$     PRINT*,"GMatrixInitialisation(", my_rank, ") error ", IErr, &
!!$          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
!!$     RETURN
!!$  ENDIF
!!$
!!$  RUniqueKeyComplete = RUniqueKey(:IUniqueKey,:)
!!$
!!$  DEALLOCATE(&
!!$       RUniqueKey,STAT=IErr)
!!$  IF( IErr.NE.0 ) THEN
!!$     PRINT*,"GMatrixInitialisation(", my_rank, ") error ", IErr, &
!!$          " in DEALLOCATE() of DYNAMIC variables Reflection Matrix"
!!$     RETURN
!!$  ENDIF
!!$ 
!!$  ALLOCAte(&
!!$      RUniqueKey(IUniqueKey,5),&
!!$      STAT = IErr)
!!$  IF( IErr.NE.0 ) THEN
!!$     PRINT*,"GMatrixInitialisation(", my_rank, ") error ", IErr, &
!!$          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
!!$     RETURN
!!$  ENDIF
!!$
!!$  RUniqueKey  = RUniqueKeyComplete
!!$
!!$  
!!$
!!$  DEALLOCATE(&
!!$       RUniqueKeyComplete,STAT=IErr)
!!$  IF( IErr.NE.0 ) THEN
!!$     PRINT*,"GMatrixInitialisation(", my_rank, ") error ", IErr, &
!!$          " in DEALLOCATE() of DYNAMIC variables Reflection Matrix"
!!$     RETURN
!!$  ENDIF
!!$  PRINT*,"IUniqueKey = ",IUniqueKey
  
  RgMatMag = RgMatMag/TWOPI
  
END SUBROUTINE GMatrixInitialisation

SUBROUTINE DetermineSymmetryRelatedUgs (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,ierr,knd,IArrayLoc
  INTEGER(IKIND),DIMENSION(:,:),ALLOCATABLE :: &
       IUgStrength
  INTEGER(IKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       ISymmetryhklRelations
  INTEGER(IKIND),DIMENSION(THREEDIM) :: &
       Ihklval
  INTEGER(IKIND),DIMENSION(1) :: &
       Ihklmax,h,k,l

  ALLOCATE(&
       IUgStrength(NINT(((FOUR*IHKLMAXValue)+ONE)**THREE),2),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DetermineSymmetryRelatedUgs(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     RETURN
  ENDIF

  ALLOCATE( &  
       ISymmetryhklRelations(nReflections,nReflections,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     RETURN
  ENDIF

  IUgStrength(:,2) = 0
  DO ind = 1,SIZE(IUgStrength,DIM=1)
     IUgStrength(ind,1) = ind
  END DO

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"GMatrixInitialisation()"
  END IF

  DO ind=1,nReflections
     DO jnd=1,nReflections
        
        
        Ihklval = RHKL(ind,:)-RHKL(jnd,:)
!!$        Ihklval = ihklval + (2*IHKLMAXValue)
        IArrayLoc = (&
             (ihklval(1)+(2*IHKLMAXValue))*(2*IHKLMAXValue)**2+&
             (ihklval(2)+(2*IHKLMAXValue))*(2*IHKLMAXValue)+&
             (ihklval(3)+(2*IHKLMAXValue))+&
             1)
        IUgStrength(IArrayLoc,2)= IUgStrength(IArrayLoc,2)+1
        ISymmetryRelations(ind,jnd) = IArrayLoc
        ISymmetryhklRelations(ind,jnd,:) = Ihklval     
     ENDDO
  ENDDO
  Ihklmax = MAXLOC(IUgStrength(:,2))-1
!!$  Ihklmax-1
  PRINT*,IHKLMAXValue
  h = Ihklmax/((2*IHKLMAXValue)**2)-2*IHKLMAXValue
  PRINT*,"h = ",h
  k = MOD(ihklmax-((h+(2*IHKLMAXValue))*(2*IHKLMAXValue)),(2*IHKLMAXValue))
  PRINT*,"k = ",k
  l = ihklmax-&
       ((((Ihklmax-1)/(2*IHKLMAXValue)**2)-2*IHKLMAXValue)+&
       MOD(ihklmax,((Ihklmax-1)/(2*IHKLMAXValue)**2)-2*IHKLMAXValue))
  PRINT*,Ihklmax,MAXVAL(IUgStrength(:,2))
  !RgMatMag = RgMatMag/TWOPI

  DO ind=1,6
     PRINT*,ISymmetryRelations(ind,:6)
  END DO
  
  DO ind=1,6
     DO jnd=1,6
        PRINT*,ISymmetryhklRelations(ind,jnd,:)
     END DO
  END DO
  
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
  
  INTEGER(IKIND) ind,jnd,ierr, currentatom, iAtom,imaxj,IFound,ICount
  INTEGER(IKIND),DIMENSION(2) :: &
       IPos
  COMPLEX(CKIND) CVgij
  COMPLEX(CKIND), DIMENSION(:,:),ALLOCATABLE :: &
       CUgMatUnique
  REAL(RKIND) RAtomicFormFactor
  REAL(RKIND) :: &
       RMeanInnerPotentialVolts

  ALLOCATE(&
       CUgMatUnique(nReflections,nReflections),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"UgCalculation(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"UgCalculation()"
  END IF
  
  DO ind=1,nReflections
     imaxj = ind
     DO jnd=1,imaxj
        
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
  RMeanInnerPotentialVolts = ((RMeanInnerCrystalPotential*RPlanckConstant**2)/ &
       (TWO*RElectronMass*RElectronCharge*TWOPI**2))*&
       RAngstromConversion*RAngstromConversion

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"UgCalculation(",my_rank,") RMeanInnerCrystalPotential = ",RMeanInnerCrystalPotential,RMeanInnerPotentialVolts
  END IF
  
  DO ind=1,nReflections
     CUgMat(ind,ind)=CUgMat(ind,ind)-RMeanInnerCrystalPotential
  ENDDO

  CUgMat = CUgMat + CONJG(TRANSPOSE(CUgMat))

  WHERE(REAL(REAL(CUgMat)).LT.TINY)
     CUgmat = ZERO + REAL(AIMAG(CUgMat))*CIMAGONE
  END WHERE

  WHERE(REAL(AIMAG(CUgMat)).LT.TINY)
     CUgmat = REAL(REAL(CUgMat)) + CZERO
  END WHERE

!!$  WHERE(CUgMat.EQ.CUgMat(1,2))
!!$     CUgMatUnique = CONE
!!$  ELSEWHERE
!!$     CUgMatUnique = CZERO
!!$  END WHERE
!!$
!!$  IF(my_rank.EQ.0) THEN
!!$     IFound = 1
!!$     ICount = 0
!!$     DO WHILE (IFound.EQ.1)
!!$        IF(MAXVAL(REAL(CUgMatUnique)).EQ.ONE) THEN
!!$           IPos = MAXLOC(REAL(CUgMatUnique,RKIND))
!!$           IFound = 1
!!$           ICount = ICount + 1
!!$           !PRINT*,IPos,CUgMat(IPos(1),IPos(2))
!!$           CUgMatUnique(IPos(1),IPos(2))=CZERO
!!$        ELSE
!!$           IFound = 0
!!$        END IF
!!$     END DO
!!$     PRINT*,"Total Symmetry Related Ugs",ICount,nReflections**2-ICount
!!$  END IF
  
  DO ind=1,10
     PRINT*,RHKL(ind,:),RgVecMatT(ind,:),CUgmat(ind,:1)
  END DO
  

END SUBROUTINE UgCalculation

SUBROUTINE UgAddAbsorption(IErr)         


  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 
  
  INTEGER(IKIND) IErr

  CUgMatPrime = CZERO

  CUgMatPrime = CUgMatPrime+ABS(CUgMat)*(RAbsorptionPercentage/100.D0)*CIMAGONE

!!$  IF (IAbsorbFlag.EQ.1) THEN
!!$     CUgMat(ind,jnd) = &
!!$          CUgMat(ind,jnd) + &
!!$          ABS(CUgMat(ind,jnd))*(RAbsorptionPercentage/100.D0)*CIMAGONE       
!!$  END IF
  
END SUBROUTINE UgAddAbsorption
