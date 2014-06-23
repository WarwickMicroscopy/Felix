!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! FelixSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
!
! (C) 2013/14, all rights reserved
!
! Version: VERSION
! Date:    DATE
! Time:    TIME
! Status:  STATUS
! Build:   BUILD
! Author:  AUTHOR
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: util.f90,v 1.95 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!--------------------------------------------------------------------
!	ReSort:
!
!	sort the Lyapunov eigenvalues s.t. the largest comes first. RESORT()
!	is based on ShellSort from "Numerical Recipes", routine SHELL().
!---------------------------------------------------------------------

SUBROUTINE ReSortHKL( RHKLarray, N )

  USE MyNumbers
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER (IKIND) N
  REAL(RKIND) RHKLarray(N,THREEDIM)
  REAL(RKIND) RhklarraySearch(THREEDIM), RhklarrayCompare(THREEDIM)
  
  REAL(KIND=RKIND) ALN2I, LocalTINY
  PARAMETER (ALN2I=1.4426950D0, LocalTINY=1.D-5)
  
  INTEGER (IKIND) NN,M,L,K,J,I,LOGNB2, index
  REAL(KIND=RKIND) dummy

  IF((IWriteFLAG.EQ.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
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
        RhklarraySearch = Rhklarray(L,1)*RarVecO + &
             Rhklarray(L,2)*RbrVecO + &
             Rhklarray(L,3)*RcrVecO    
        RhklarrayCompare = Rhklarray(I,1)*RarVecO + &
             Rhklarray(I,2)*RbrVecO + &
             Rhklarray(I,3)*RcrVecO
        IF( &
             DOT_PRODUCT(RHKLarraySearch(:),RHKLarraySearch(:)) .LT. &
             DOT_PRODUCT(RHKLarrayCompare(:),RHKLarrayCompare(:))) THEN
           DO 100 index=1,THREEDIM
              dummy        = RHKLarray(I,index)
              RHKLarray(I,index)= RHKLarray(L,index)
              RHKLarray(L,index)= dummy
100        ENDDO
           
           I=I-M
           IF(I.GE.1) GOTO 3
        ENDIF
11   ENDDO
12 ENDDO
  
  !PRINT*,"Finishing ResortHKL"

  !	PRINT*,"array0(1),array0(N)",array0(1),array0(N)
  RETURN

END SUBROUTINE ReSortHKL

!---------------------------------------------------------------------
SUBROUTINE CONVERTAtomName2Number(name, number, IErr)

  USE IPara
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER IErr, ind, number
  CHARACTER*2 name

  INTEGER, PARAMETER :: NElements=103

  CHARACTER*2 A(NElements)

  DATA A/" H", "He", "Li", "Be", " B", " C", " N", "O", "F", "Ne", &
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", &
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", &
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", &
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",& 
        "Md","No","Lr"/

  DO ind=1,NElements
     IF(TRIM(name)==TRIM(A(ind))) THEN
        number= ind
        IF((IWriteFLAG.EQ.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"DBG: name, number ", name, number
        END IF
           RETURN
     ENDIF
  ENDDO

  PRINT*,"CONVERTAtomName2Number(): could not find index for atom of name ", name
  IErr=1
  RETURN

  PRINT*,"DBG: name, number ", name, number
END SUBROUTINE CONVERTAtomName2Number

!---------------------------------------------------------------------
SUBROUTINE ImageMask (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER ind,jnd, ierr
  REAL(RKIND) :: Rradius, RImageRadius
  
  IF((IWriteFLAG.EQ.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: ImageMask()"
  END IF

  IPixelTotal =0
  SELECT CASE (IMaskFLAG)
  CASE(0) ! circle
     DO ind=1,2*IPixelCount
        DO jnd=1,2*IPixelCount
           Rradius= (ind-(REAL(IPixelCount,RKIND)+0.5))**2 + &
                (jnd-(REAL(IPixelCount,RKIND)+0.5))**2
           Rradius=SQRT(DBLE(Rradius))
           RImageRadius = IPixelCount+0.5
           IF(Rradius.LE.RImageRadius) THEN
              RMask(jnd,ind) = 1
              IPixelTotal= IPixelTotal + 1
              
           ELSE
              RMask(jnd,ind) = 0
           ENDIF
        ENDDO
     ENDDO
  CASE(1) ! square
     RMask = 1
     IPixelTotal = (2*IPixelCount)**2
  END SELECT
  
  ALLOCATE( &
       IPixelLocations(IPixelTotal,2), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Imagemask(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  IPixelTotal = 0
 
  SELECT CASE (IMaskFLAG)
  CASE(0) ! circle
     DO ind=1,2*IPixelCount
        DO jnd=1,2*IPixelCount
           Rradius= (ind-(REAL(IPixelCount,RKIND)+0.5))**2 + &
                (jnd-(REAL(IPixelCount,RKIND)+0.5))**2
           Rradius=SQRT(DBLE(Rradius))
           RImageRadius = IPixelCount+0.5
           IF(Rradius.LE.RImageRadius) THEN
              !RMask(jnd,ind) = 1
              IPixelTotal= IPixelTotal + 1
              IPixelLocations(IPixelTotal,1) = ind
              IPixelLocations(IPixelTotal,2) = jnd
           ELSE
              !RMask(jnd,ind) = 0
           ENDIF
        ENDDO
     ENDDO
  CASE(1) ! square
     IPixelTotal = 0
     DO ind = 1,2*IPixelCount
        DO jnd = 1,2*IPixelCount
           IPixelTotal = IPixelTotal+1
           IPixelLocations(IPixelTotal,1) = ind
           IPixelLocations(IPixelTotal,2) = jnd
        END DO
     END DO
  END SELECT
  
END SUBROUTINE ImageMask

!---------------------------------------------------------------------
SUBROUTINE CountTotalAtoms(IErr)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER ind,jnd,knd,hnd,ierr, ifullind, iuniind
  LOGICAL Lunique

  
  IF((IWriteFLAG.EQ.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"CountTotalAtoms()"
  END IF

  ALLOCATE( &
       RFullAtomicFracCoordVec( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1),&
       THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE( &
       SFullAtomicNameVec( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  
  RFullAtomicFracCoordVec = ZERO
  
  IF((IWriteFLAG.GE.10.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"SIZE OF RFULLATOMICFRACCOORDVEC = ",SIZE(RFullAtomicFracCoordVec,1)
  END IF

  DO ind=1, SIZE(RSymVec,DIM=1)
     
     DO jnd=1, SIZE(RAtomSiteFracCoordVec,DIM=1)
       
        Ifullind= SIZE(RSymVec,1)*(jnd-1) + ind
        
        RFullAtomicFracCoordVec(Ifullind,:)= &
             MATMUL(RSymMat(ind,:,:),RAtomSiteFracCoordVec(jnd,:)) &
             + RSymVec(ind,:)
        SFullAtomicNameVec(Ifullind) = SAtomName(jnd)
        
        ! renormalize such that all values are non-negative
        DO knd=1,THREEDIM
           IF( RFullAtomicFracCoordVec(Ifullind,knd) .LT. ZERO) THEN
              RFullAtomicFracCoordVec(Ifullind,knd)= &
                   RFullAtomicFracCoordVec(Ifullind,knd)+1.D0
           ENDIF
        ENDDO
        
     ENDDO
     
  ENDDO

  DO ind = 1,SIZE(RFullAtomicFracCoordVec,DIM=1)
     DO jnd = 1,SIZE(RFullAtomicFracCoordVec,DIM=2)
        IF (RFullAtomicFracCoordVec(ind,jnd).LT.ZERO) THEN
           RFullAtomicFracCoordVec(ind,jnd) = &
                RFullAtomicFracCoordVec(ind,jnd) + ONE
        END IF
        IF (RFullAtomicFracCoordVec(ind,jnd).GE.ONE) THEN
           RFullAtomicFracCoordVec(ind,jnd) = &
                RFullAtomicFracCoordVec(ind,jnd) - ONE
        END IF
     END DO
  END DO

  
  ! Calculate the set of unique fractional atomic positions
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"CountTotalAtoms(",my_rank,") ITotalAtoms = ",ITotalAtoms
  END IF

  ALLOCATE( &
       MNP(1000,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       SMNP(1000), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  
  MNP = ZERO
  
  MNP(1,:)= RFullAtomicFracCoordVec(1,:)
  SMNP(1)= SFullAtomicNameVec(1)
  
  Iuniind = 1
  
  IF(ITotalAtoms.EQ.0)THEN
     
     DO ind=2,SIZE(RFullAtomicFracCoordVec,1)
        
        DO jnd=1,Iuniind
           
           IF ( RFullAtomicFracCoordVec(ind,1) .EQ. MNP(jnd,1) .AND. &
                RFullAtomicFracCoordVec(ind,2) .EQ. MNP(jnd,2) .AND. &
                RFullAtomicFracCoordVec(ind,3) .EQ. MNP(jnd,3) .AND. &
                SFullAtomicNameVec(ind) .EQ. SMNP(jnd) ) THEN
              !this seems NOT a unique coordinate
              Lunique=.FALSE.
              EXIT
           ENDIF
           Lunique=.TRUE.
        ENDDO
        
        IF(Lunique .EQV. .TRUE.) THEN
           Iuniind=Iuniind+1
           MNP(Iuniind,:)= RFullAtomicFracCoordVec(ind,:)
           SMNP(Iuniind)= SFullAtomicNameVec(ind)
        ENDIF
        
     ENDDO
     
     ITotalAtoms = Iuniind
     
  END IF
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"CountTotalAtoms(",my_rank,") ITotalAtoms = ",ITotalAtoms
  END IF

  DEALLOCATE( &
       MNP,SMNP, &
       RFullAtomicFracCoordVec, &
       SFullAtomicNameVec,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
  
END SUBROUTINE CountTotalAtoms

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
