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

! All procedures conatained in this file/module:
! PhaseCorrelate()
! ReSortUgs()
! ResidualSumofSquares()
! Normalised2DCrossCorrelation()
! MaskedCorrelation()
! SortHKL()
! GreatestCommonDivisor()
! Lorentzian()
! Gaussian()

!>
!! Module-description: 
!!
MODULE utilities_mod

  IMPLICIT NONE

  CONTAINS

  !>
  !! Procedure-description: Plan and Execute the fft of the Simulated Data. 
  !! Plan and Execute the fft of the Experimental Data. Plan and Execute the
  !! inverse fft of the phase correlation.
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  REAL(RKIND) FUNCTION PhaseCorrelate(RImageSim,RImageExpiDummy,IErr,IXsizeIn,IYSizeIn)
    
    USE MyNumbers
    
    USE SConst; USE IConst
    USE IPara; USE RPara
    USE IChannels
    
    USE MPI
    USE MyMPI
    USE MyFFTW

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr,IXsizeIn,IYSizeIn
    REAL(RKIND),DIMENSION(IXSizeIn,IYSizeIn) :: RImageExpiDummy,RImageSim
    
    type(C_PTR) :: Iplan
    real(C_DOUBLE), pointer :: RImageSimDummy(:,:)
    type(C_PTR) :: p1,p2,p3,p4
    complex(C_DOUBLE_COMPLEX), pointer :: CDummy1(:,:),CDummy2(:,:),CCorrelatedImage(:,:)
    integer(C_INT) :: IX,IY
    
    IX = IXSizeIn
    IY = IYSizein

    !CALL message ( "IX, IY =", (/ IX,IY /) )

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

    ! CALL message ( "RImageSimDummy", RImageSimDummy(:2,:2) )

    ! Plan and Execute the fft of the Simulated Data 
    Iplan = FFTW_PLAN_DFT_r2c_2D(IX,IY,RImageSimDummy,CDummy1,FFTW_ESTIMATE)
    CALL FFTW_EXECUTE_DFT_R2C(Iplan,RImageSimDummy,CDummy1)
    CALL FFTW_DESTROY_PLAN(Iplan)

    ! Set the dummy array to the input experimental data
    RImageSimDummy = RImageExpiDummy 
    ! CALL message ( "RImageSimDummy", RImageSimDummy(:2,:2) )
   
    ! Plan and Execute the fft of the Experimental Data 
    Iplan = FFTW_PLAN_DFT_R2C_2D(IX,IY,RImageSimDummy,CDummy2,FFTW_ESTIMATE)
    CALL FFTW_EXECUTE_DFT_R2C(Iplan,RImageSimDummy,CDummy2)
    CALL FFTW_DESTROY_PLAN(Iplan)

    !Calculate the Phase Correlation
    WHERE(ABS(CDummy1*CONJG(CDummy2)).NE.ZERO)
       CCorrelatedImage = (CDummy1*CONJG(CDummy2))/&
            (ABS(CDummy1*CONJG(CDummy2)))
    ELSEWHERE 
       CCorrelatedImage = CZERO
    END WHERE
    
    ! Plan and Execute the inverse fft of the phase correlation
    Iplan = FFTW_PLAN_DFT_C2R_2D(IX,IY,CCorrelatedImage,RImageSimDummy,FFTW_ESTIMATE)
    CALL FFTW_EXECUTE_DFT_C2R(Iplan,CCorrelatedImage,RImageSimDummy)
    CALL FFTW_DESTROY_PLAN(Iplan)

    ! RCrossCorrelation = MAXVAL(RImageSimDummy)/(IX*IY)
    PhaseCorrelate = MAXVAL(RImageSimDummy)/(IX*IY)
    IOffset = MAXLOC(RImageSimDummy)
    
    ! CALL message ( "RImageSimDummy", RImageSimDummy(:2,:2) )

    call fftw_free(p1)
    call fftw_free(p2)
    call fftw_free(p3)
    call fftw_free(p4)
    
  END FUNCTION  PhaseCorrelate

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !>
  !! Procedure-description: 
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE ReSortUgs( ISymmetryIntegers,CUgs, N )
    
    USE MyNumbers
    
    USE SConst; USE IConst
    USE IPara; USE RPara
    USE IChannels
    USE message_mod
    USE MPI
    USE MyMPI

    IMPLICIT NONE

    INTEGER(IKIND) :: N,IDummy,ISymmetryIntegers(N)
    INTEGER(IKIND) :: NN,M,L,K,J,I,LOGNB2
    REAL(RKIND) :: RhklarraySearch(ITHREE),RhklarrayCompare(ITHREE)
    REAL(RKIND) :: ALN2I,LocalTINY
    COMPLEX(CKIND) :: CUgSearch,CUgCompare,Cdummy,CUgs(N)
    PARAMETER (ALN2I=1.4426950D0, LocalTINY=1.D-5)

    CALL message ( LL, "Sorting Ugs" )
    
    LOGNB2=INT(LOG(REAL(N))*ALN2I+LocalTINY)
    M=N
    DO 12 NN=1,LOGNB2
       M=M/2
       K=N-M
       DO 11 J=1,K
          I=J
3         CONTINUE
          L=I+M     
          CUgSearch = CUgs(L)
          CUgCompare = CUgs(I)
          IF( (ABS(CUgSearch)).GT.(ABS(CUgCompare)) ) THEN ! RB sort on modulus ABS
                Cdummy = CUgs(I)
                CUgs(I)= CUgs(L)
                Cugs(L)= Cdummy
                Idummy = ISymmetryIntegers(I)
                ISymmetryIntegers(I)= ISymmetryIntegers(L)
                ISymmetryIntegers(L)= Idummy
             I=I-M
             IF(I.GE.1) GOTO 3
          END IF
11     ENDDO
12   ENDDO
    
    RETURN

  END SUBROUTINE ReSortUgs

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !>
  !! Procedure-description: 
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  REAL(RKIND) FUNCTION ResidualSumofSquares(RImage1,RImage2,IErr)
    !todo - this is excessively small underused subroutine & the message should be on debug mode
    USE MyNumbers
    
    USE SConst; USE IConst
    USE IPara; USE RPara
    USE IChannels
    USE message_mod
    USE MPI
    USE MyMPI

    IMPLICIT NONE
    
    INTEGER(IKIND) :: &
         IErr
    REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: &
         RImage1,RImage2

    RImage2 =  RImage2/(2.0**16.0)

    CALL message ( LM, "Residual Sum of Squares Before",&
          (/ ResidualSumofSquares,MAXVAL(RImage1),MAXVAL(RImage2) /) )

    ResidualSumofSquares = SUM((RImage2-RImage1)**2)

    CALL message ( LM, "Residual Sum of Squares After",ResidualSumofSquares )

  END FUNCTION ResidualSumofSquares

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !>
  !! Procedure-description: 
  !!
  !! Major-Authors: 'kidwhizz' (2015), Richard Beanland (2016)
  !!
  REAL(RKIND) FUNCTION Normalised2DCrossCorrelation(Rimg1,Rimg2,IErr)

    USE MyNumbers
    
    USE SConst; USE IConst
    USE IPara; USE RPara
    USE IChannels
    USE MPI
    USE MyMPI

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr
    REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: Rimg1,Rimg2
    REAL(RKIND) :: Rimg1Mean,Rimg2Mean,Rimg1StDev,Rimg2StDev,RPixelTotal,&
                   Rimg1Min,Rimg2Min,Rimg1Max,Rimg2Max
    
    RPixelTotal = REAL(2*IPixelCount*2*IPixelCount,RKIND)  
    Rimg1Mean = SUM(Rimg1)/RPixelTotal
    Rimg2Mean = SUM(Rimg2)/RPixelTotal
    Rimg1StDev=SQRT(SUM(((Rimg1-Rimg1Mean)**2)/RPixelTotal))
    Rimg2StDev=SQRT(SUM(((Rimg2-Rimg2Mean)**2)/RPixelTotal))
    Normalised2DCrossCorrelation=SUM(((Rimg1-Rimg1Mean)*(Rimg2-Rimg2Mean)))/&
         (Rimg1StDev*Rimg2StDev*RPixelTotal)

  END FUNCTION Normalised2DCrossCorrelation

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: 
  !!
  !! Major-Authors: Richard Beanland (2016)
  !!
  REAL(RKIND) FUNCTION MaskedCorrelation(Rimg1,Rimg2,RBinaryMask,IErr)

    USE MyNumbers
    
    USE SConst; USE IConst
    USE IPara; USE RPara
    USE IChannels
    USE message_mod
    USE MPI
    USE MyMPI

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr
    REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: Rimg1,Rimg2,RBinaryMask
    REAL(RKIND) :: Rimg1Mean,Rimg2Mean,Rimg1StDev,Rimg2StDev,RPixelTotal,&
                   Rimg1Min,Rimg2Min,Rimg1Max,Rimg2Max
    CHARACTER*200 :: path

    IF (SUM(RBinaryMask).GT.ZERO) THEN
      RPixelTotal = SUM(RBinaryMask)!Mask has value 1 for pixels of interest and zero elsewhere
      CALL message ( LL, dbg6, "Mask pixels = ", RPixelTotal )
    ELSE
      RPixelTotal = REAL(2*IPixelCount*2*IPixelCount,RKIND)
      CALL message ( LL, dbg6, "Empty, mask" )
    END IF  
    Rimg1=Rimg1*RBinaryMask
    Rimg2=Rimg2*RBinaryMask
    Rimg1Mean = SUM(Rimg1)/RPixelTotal
    Rimg2Mean = SUM(Rimg2)/RPixelTotal
    Rimg1StDev=SQRT(SUM(((Rimg1-Rimg1Mean)**2)/RPixelTotal))
    Rimg2StDev=SQRT(SUM(((Rimg2-Rimg2Mean)**2)/RPixelTotal))
    MaskedCorrelation=SUM(((Rimg1-Rimg1Mean)*(Rimg2-Rimg2Mean)))/&
         (Rimg1StDev*Rimg2StDev*RPixelTotal)

  END FUNCTION MaskedCorrelation

  !>
  !! Procedure-description: Sorts array into descending order
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE SortHKL( Rhklarray,N,IErr )

    !--------------------------------------------------------------------
    !	Sort:
    !	sort s.t. the largest comes first. RESORT()
    !	is based on ShellSort from "Numerical Recipes", routine SHELL().
    !---------------------------------------------------------------------  

    USE MyNumbers
      
    USE SConst; USE IConst
    USE IPara; USE RPara

    USE IChannels

    USE MPI
    USE MyMPI

    IMPLICIT NONE

    INTEGER (IKIND) :: IErr,NN,M,L,K,J,I,LOGNB2,ind
    INTEGER (IKIND),INTENT(IN) :: N
    REAL(RKIND),INTENT(INOUT) :: Rhklarray(N,ITHREE)
    REAL(RKIND) :: RhklarraySearch(ITHREE), RhklarrayCompare(ITHREE)
    REAL(RKIND) :: ALN2I, LocalTINY
    PARAMETER (ALN2I=1.4426950D0, LocalTINY=1.D-5)
    
    REAL(RKIND) :: dummy
   
    NN = 0
    M = 0
    L = 0
    K = 0
    J = 0
    I = 0
    LOGNB2 = 0
    ind = 0
    RhklarraySearch = 0.0D0
    RhklarrayCompare = 0.0D0
    dummy = 0.0D0

    LOGNB2=INT(LOG(REAL(N))*ALN2I+LocalTINY)
    M=N
    DO 12 NN=1,LOGNB2
       M=M/2
       K=N-M
       DO 11 J=1,K
          I=J
3         CONTINUE
          L=I+M
          RhklarraySearch = Rhklarray(L,1)*RarVecO + &
               Rhklarray(L,2)*RbrVecO + &
               Rhklarray(L,3)*RcrVecO    
          RhklarrayCompare = Rhklarray(I,1)*RarVecO + &
               Rhklarray(I,2)*RbrVecO + &
               Rhklarray(I,3)*RcrVecO
          IF( DOT_PRODUCT(RhklarraySearch(:),RhklarraySearch(:)) .LT. &
              DOT_PRODUCT(RhklarrayCompare(:),RhklarrayCompare(:))) THEN
             DO 100 ind=1,ITHREE
                dummy        = Rhklarray(I,ind)
                Rhklarray(I,ind)= Rhklarray(L,ind)
                Rhklarray(L,ind)= dummy
100          ENDDO
             I=I-M
             IF(I.GE.1) GOTO 3
          ENDIF
11     ENDDO
12   ENDDO
    
    RETURN

  END SUBROUTINE SortHKL

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !>
  !! Procedure-description: Defines a Lorentzian Distribution for any parameter input
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  FUNCTION Lorentzian(FWHM,x,x_0,offset)

    USE MyNumbers
    
    USE RPara; 

    IMPLICIT NONE

    REAL(RKIND):: FWHM,x,x_0,offset,LORENTZIAN

    LORENTZIAN = FWHM/(((x+x_0)**2)+offset)
    
  END FUNCTION Lorentzian

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !>
  !! Procedure-description: Defines a Gaussian distribution for any parameter input 
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  FUNCTION Gaussian(height,x,peakcentre,standarddeviation,intercept)

    USE MyNumbers

    IMPLICIT NONE

    REAL(RKIND):: height,x,peakcentre,standarddeviation,intercept,gaussian

    Gaussian = height*exp(-(((x-peakcentre)**2)/(2*(standarddeviation**2))))+ intercept
    
  END FUNCTION Gaussian

END MODULE utilities_mod
