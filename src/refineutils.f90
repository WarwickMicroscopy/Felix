!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all right reserved
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

SUBROUTINE DetermineFluxSteps(IErr)
 
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  !--------------------------------------------------------------------
  ! local variable definitions
  !--------------------------------------------------------------------
  
  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr
  
  SELECT CASE (IRefineModeFLAG)
     
  CASE(0)

     IFluxIterationSteps = ((RFinalDebyeWallerFactor-RInitialDebyeWallerFactor)/&
          RDeltaDebyeWallerFactor +1)**SIZE(IElementList)
     
  CASE(1)
     
     IFluxIterationSteps = (RUpperBoundUgChange-RLowerBoundUgChange)/RDeltaUgChange + 1

  CASE(2)

     IFluxIterationSteps = 1
     
     IThicknessCount= (RFinalThickness- RInitialThickness)/RDeltaThickness + 1
     
  CASE DEFAULT
     
     IErr = 1
     IF( IErr.NE.0 ) THEN
        PRINT*,"DetermineFluxSteps(", my_rank, ") error ", IErr, &
             " Refinement Mode Not Defined"
        RETURN
     ENDIF
     
  END SELECT
END SUBROUTINE DetermineFluxSteps

SUBROUTINE AngularOffset (RImageSim,RImageExpi,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER ind,jnd, ierr
  REAL(RKIND),DIMENSION(IImageSizeXY(1),IImageSizeXY(2)) :: &
       RImageSim,RImageExpi  

END SUBROUTINE AngularOffset

SUBROUTINE CorrectCentreOffset (RUncorrImage,RCorrImage,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER ind,jnd, ierr
  REAL(RKIND),DIMENSION(IImageSizeXY(1),IImageSizeXY(2)) :: &
       RUnCorrImage
  REAL(RKIND),DIMENSION(IImageSizeXY(1)-IOffset(1),IImageSizeXY(1)-IOffset(1)),INTENT(OUT) :: &
       RCorrImage

  RCorrImage = RUncorrImage((IOffset(1)+1):,(IOffset(2)+1):) 
  
END SUBROUTINE CorrectCentreOffset

SUBROUTINE PhaseCorrelate(RImageSim,RImageExpi,IErr,IXsizeIn,IYSizeIn)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  
  USE MPI
  USE MyMPI
  USE MyFFTW

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,IXsizeIn,IYSizeIn

  REAL(RKIND),DIMENSION(IXSizeIn,IYSizeIn) :: &
       RImageExpi,RImageSim
  
  type(C_PTR) :: &
       Iplan

  real(C_DOUBLE), pointer :: &
       RImageSimDummy(:,:)
  type(C_PTR) :: & 
       p1,p2,p3,p4
  complex(C_DOUBLE_COMPLEX), pointer :: &
       CDummy1(:,:),CDummy2(:,:),CCorrelatedImage(:,:)
  integer(C_INT) :: IX,IY
  
  IX = IXSizeIn
  IY = IYSizein

  !PRINT*,"IX, IY =",IX,IY

  p1 = fftw_alloc_real(INT(IXSizeIn*IYSizeIn, C_SIZE_T))
  p2 = fftw_alloc_complex(INT(IXSizeIn*IYSizeIn, C_SIZE_T))
  p3 = fftw_alloc_complex(INT(IXSizeIn*IYSizeIn, C_SIZE_T))
  p4 = fftw_alloc_complex(INT(IXSizeIn*IYSizeIn, C_SIZE_T))
  call c_f_pointer(p1, RImageSimDummy, [IXSizeIn,IYSizeIn])
  call c_f_pointer(p2, CDummy1, [IXSizeIn,IYSizeIn])
  call c_f_pointer(p3, CDummy2, [IXSizeIn,IYSizeIn])
  call c_f_pointer(p4, CCorrelatedImage, [IXSizeIn,IYSizeIn])
  !...use arr and arr(i,j) as usual...
  
  ! Set the dummy array to the input simulated data

  RImageSimDummy = RImageSim


  !PRINT*,RImageSimDummy(:2,:2)

  ! Plan and Execute the fft of the Simulated Data 
  
  Iplan = FFTW_PLAN_DFT_r2c_2D(IX,IY,RImageSimDummy,CDummy1,FFTW_ESTIMATE)
  CALL FFTW_EXECUTE_DFT_R2C(Iplan,RImageSimDummy,CDummy1)
  CALL FFTW_DESTROY_PLAN(Iplan)

  ! Set the dummy array to the input experimental data

  RImageSimDummy = RImageExpi


  !PRINT*,RImageSimDummy(:2,:2)
 
  ! Plan and Execute the fft of the Experimental Data 

  Iplan = FFTW_PLAN_DFT_R2C_2D(IX,IY,RImageSimDummy,CDummy2,FFTW_ESTIMATE)

  CALL FFTW_EXECUTE_DFT_R2C(Iplan,RImageSimDummy,CDummy2)
  CALL FFTW_DESTROY_PLAN(Iplan)

  !Calculate the Phase Correlation

  CCorrelatedImage = (CDummy1*CONJG(CDummy2))/&
       (&
       ABS(CDummy1*CONJG(CDummy2)+&
       CMPLX(TINY,TINY,C_DOUBLE_COMPLEX)))

  
  ! Plan and Execute the inverse fft of the phase correlation

  Iplan = FFTW_PLAN_DFT_C2R_2D(IX,IY,CCorrelatedImage,RImageSimDummy,FFTW_ESTIMATE)

  CALL FFTW_EXECUTE_DFT_C2R(Iplan,CCorrelatedImage,RImageSimDummy)

  CALL FFTW_DESTROY_PLAN(Iplan)

  
  RCrossCorrelation = MAXVAL(RImageSimDummy)
  IOffset = MAXLOC(RImageSimDummy)
  
  !PRINT*,RImageSimDummy(:2,:2)

  call fftw_free(p1)
  call fftw_free(p2)
  call fftw_free(p3)
  call fftw_free(p4)
  
END SUBROUTINE PhaseCorrelate

!!$SUBROUTINE BiCubicResampling(RImin,IImSize,RImout,RScaleFactor,IErr)
!!$  
!!$USE MyNumbers
!!$  
!!$  USE CConst; USE IConst
!!$  USE IPara; USE RPara
!!$  USE IChannels
!!$  
!!$  USE MPI
!!$  USE MyMPI
!!$
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER(IKIND) :: &
!!$       ind,jnd, ierr
!!$  INTEGER(IKIND),DIMENSION(2),Intent(IN) :: &
!!$       IImSize
!!$  INTEGER(IKIND),DIMENSION(2):: &
!!$       IImNewSize
!!$  REAL(RKIND),DIMENSION(IImSize(2),IImSize(2)),INTENT(IN) :: &
!!$       RImin
!!$  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT) :: &
!!$       RImout
!!$  REAL(RKIND) :: &
!!$       RScaleFactor,RrowFunction,RTotalInterpolation
!!$  REAL(RKIND),DIMENSION(:) :: &
!!$       Ralist,Rblist
!!$
!!$  IF(RScaleFactor.LT.ONE) THEN
!!$     IImNewSize(1) = FLOOR(IImSize(1)*RScaleFactor)
!!$     IImNewSize(2) = FLOOR(IImSize(2)*RScaleFactor)
!!$  ELSE
!!$     IImNewSize(1) = CEILING(IImSize(1)*RScaleFactor)
!!$     IImNewSize(2) = CEILING(IImSize(2)*RScaleFactor)
!!$  END IF
!!$
!!$  ALLOCATE(&
!!$       Ralist(IImNewSize(1)),&
!!$       STAT=IErr)
!!$  IF( IErr.NE.0 ) THEN
!!$     PRINT*,"BiCubicResampling(", my_rank, ") error ", IErr, &
!!$          " in allocation"
!!$     RETURN
!!$  ENDIF
!!$  ALLOCATE(&
!!$       RImout(IImNewSize(1),IImNewSize(2)),&
!!$       STAT=IErr)
!!$  IF( IErr.NE.0 ) THEN
!!$     PRINT*,"BiCubicResampling(", my_rank, ") error ", IErr, &
!!$          " in allocation"
!!$     RETURN
!!$  ENDIF
!!$
!!$  ALLOCATE(&
!!$       Rblist(IImNewSize(2)),&
!!$       STAT=IErr)
!!$  IF( IErr.NE.0 ) THEN
!!$     PRINT*,"BiCubicResampling(", my_rank, ") error ", IErr, &
!!$          " in allocation"
!!$     RETURN
!!$  ENDIF
!!$
!!$  Ralist = (1:IImNewSize(1))*(1/RScaleFactor)
!!$  Rblist = (1:IImNewSize(2))*(1/RScaleFactor)
!!$
!!$  DO ind=1,IImNewSize(1)
!!$     DO jnd=1,IImNewSize(2)
!!$     RrowFunction = -Rblist(ind)*(1-Rblist(ind))**2*RImin() + &
!!$          (1-2*Rblist(ind)**2+Rblist(ind)**3)*RImin() + &
!!$          Rblist(ind)*(1+Rblist(ind)-Rblist(ind)**3)*RImin() - &
!!$          (Rblist(ind)**2)*(1-Rblist(ind))*RImin()
!!$     RImout(ind,jnd) = 
!!$     END DO
!!$  END DO
!!$  
!!$ 
!!$END SUBROUTINE BiCubicResampling

SUBROUTINE ReSortUgs( ISymmetryIntegers,CUgs, N )
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER (IKIND) N,IDummy,ISymmetryIntegers(N,2)
  REAL(RKIND) RhklarraySearch(THREEDIM), RhklarrayCompare(THREEDIM)
  COMPLEX(CKIND) CUgSearch,CUgCompare,CUgs(N)
  REAL(KIND=RKIND) ALN2I, LocalTINY
  PARAMETER (ALN2I=1.4426950D0, LocalTINY=1.D-5)
  
  INTEGER (IKIND) NN,M,L,K,J,I,LOGNB2, index
  COMPLEX(CKIND) Cdummy

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ReSort()"
  END IF
  
  LOGNB2=INT(LOG(REAL(N))*ALN2I+LocalTINY)
  M=N
  DO 12 NN=1,LOGNB2
     M=M/2
     K=N-M
     DO 11 J=1,K
        I=J
3       CONTINUE
        L=I+M
        
        CUgSearch = CUgs(L)
        CUgCompare = CUgs(I)
        IF( &
             (REAL(CUgSearch**2)) .GT. &
             (REAL(CUgCompare**2))) THEN
 !          DO 100
              !IF(my_rank.eq.0) THEN
              !   PRINT*,I
              !END IF
              Cdummy = CUgs(I)
              CUgs(I)= CUgs(L)
              Cugs(L)= Cdummy
              Idummy = ISymmetryIntegers(I,2)
              ISymmetryIntegers(I,2)= ISymmetryIntegers(L,2)
              ISymmetryIntegers(L,2)= Idummy
!100        ENDDO
           
           I=I-M
           IF(I.GE.1) GOTO 3
        ENDIF
11   ENDDO
12 ENDDO
  
  !PRINT*,"Finishing ResortHKL"

  !	PRINT*,"array0(1),array0(N)",array0(1),array0(N)
  RETURN

END SUBROUTINE ReSortUgs
