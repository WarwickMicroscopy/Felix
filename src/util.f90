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

SUBROUTINE GreatestCommonDivisor(ITotalProcesses,INooDWFs,ISubgroups)

USE MyNumbers

INTEGER(IKIND) :: &
     a,b,c
INTEGER(IKIND), INTENT(IN) :: &
     ITotalProcesses,INooDWFs
INTEGER(IKIND), INTENT(OUT) :: &
     ISubgroups

a = ITotalProcesses
b = INooDWFs
c = 0


  DO                    ! now we have a <= b
     c = MOD(a, b)      !    compute c, the reminder
     IF (c == 0) EXIT   !    if c is zero, we are done.  GCD = b
     a = b              !    otherwise, b becomes a
     b = c              !    and c becomes b
  END DO                !    go back

ISubgroups = b

END SUBROUTINE GreatestCommonDivisor
