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

!--------------------------------------------------------------------
!	Sort:
!	sort s.t. the largest comes first. RESORT()
!	is based on ShellSort from "Numerical Recipes", routine SHELL().
!---------------------------------------------------------------------

SUBROUTINE SortHKL( Rhklarray,N,IErr )

  USE MyNumbers
  USE WriteToScreen
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE WriteToScreen

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER (IKIND) :: IErr,NN,M,L,K,J,I,LOGNB2,index
  INTEGER (IKIND),INTENT(IN) :: N
  REAL(RKIND),INTENT(INOUT) :: Rhklarray(N,ITHREE)
  REAL(RKIND) :: RhklarraySearch(ITHREE), RhklarrayCompare(ITHREE)
  REAL(RKIND) :: ALN2I, LocalTINY
  PARAMETER (ALN2I=1.4426950D0, LocalTINY=1.D-5)
  
  REAL(RKIND) :: dummy

  CALL Message("SortHkl",IMust,IErr)
  
  NN = 0
  M = 0
  L = 0
  K = 0
  J = 0
  I = 0
  LOGNB2 = 0
  index = 0
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
3       CONTINUE
        L=I+M
        RhklarraySearch = Rhklarray(L,1)*RarVecO + &
             Rhklarray(L,2)*RbrVecO + &
             Rhklarray(L,3)*RcrVecO    
        RhklarrayCompare = Rhklarray(I,1)*RarVecO + &
             Rhklarray(I,2)*RbrVecO + &
             Rhklarray(I,3)*RcrVecO
        IF( &
             DOT_PRODUCT(RhklarraySearch(:),RhklarraySearch(:)) .LT. &
             DOT_PRODUCT(RhklarrayCompare(:),RhklarrayCompare(:))) THEN
           DO 100 index=1,ITHREE
              dummy        = Rhklarray(I,index)
              Rhklarray(I,index)= Rhklarray(L,index)
              Rhklarray(L,index)= dummy
100        ENDDO
           
           I=I-M
           IF(I.GE.1) GOTO 3
        ENDIF
11   ENDDO
12 ENDDO
  
  RETURN

END SUBROUTINE SortHKL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CONVERTAtomName2Number(name, number, IErr)
!!$  %    Converts atomic symbols to atomic numbers, used to read cif file
  
  USE WriteToScreen
  USE IPara
  USE MPI
  USE MyMPI
  USE IConst
  USE CConst

  IMPLICIT NONE
  
  INTEGER :: &
       IErr, ind, number
  CHARACTER*2 :: &
       name

!!$  Subroutine within loop, therefore only want to print this message once
  DO WHILE (IMessageCounter.LT.1)
     CALL Message("CONVERTAtomName2Number",IMust,IErr)
     CALL Message("CONVERTAtomName2Number",IMust+IDEBUG,IErr,MessageString = "Is looping")
     IMessageCounter = IMessageCounter +1
  END DO
 

  DO ind=1,NElements
     IF(TRIM(name)==TRIM(SElementSymbolMatrix(ind))) THEN
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

END SUBROUTINE CONVERTAtomName2Number

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CountTotalAtoms(IErr)!RB now redundant

  USE MyNumbers
  USE WriteToScreen

  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER (IKIND) :: ind,jnd,knd,hnd,ierr, ifullind, iuniind
  LOGICAL :: Lunique

  CALL Message("CountTotalAtoms",IMust,IErr)
     
  IMaxPossibleNAtomsUnitCell = 0

  ALLOCATE(RAtomPosition(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1),ITHREE),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE(SAtomName(SIZE(RSymVec,1)*SIZE(RBasisAtomPosition,1)),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  RAtomPosition = ZERO

  CALL Message("CountTotalAtoms",IAllInfo,IErr, &
       MessageVariable = "Size of RAtomPosition", &
       IVariable = SIZE(RAtomPosition,1))

  !Apply symmetry elements to generate equivalent positions	   
  DO ind=1, SIZE(RSymVec,DIM=1)
    DO jnd=1, SIZE(RBasisAtomPosition,DIM=1)     
      Ifullind= SIZE(RSymVec,1)*(jnd-1) + ind
      RAtomPosition(Ifullind,:)= &
             MATMUL(RSymMat(ind,:,:),RBasisAtomPosition(jnd,:))+RSymVec(ind,:)
      SAtomName(Ifullind) = SBasisAtomName(jnd)
     ENDDO
  ENDDO
  RAtomPosition=MOD(RAtomPosition,ONE)

  ! Calculate the set of unique fractional atomic positions
  CALL Message("CountTotalAtoms",IMoreInfo,IErr, MessageVariable = "IMaxPossibleNAtomsUnitCell", &
       IVariable = IMaxPossibleNAtomsUnitCell)

  ALLOCATE(MNP(1000,ITHREE),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE(SMNP(1000),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  
  MNP = ZERO
  
  MNP(1,:)= RAtomPosition(1,:)
  SMNP(1)= SAtomName(1)
  
  Iuniind = 1
  
  IF(IMaxPossibleNAtomsUnitCell.EQ.0)THEN
     
     DO ind=2,SIZE(RAtomPosition,1)
        DO jnd=1,Iuniind
           
           IF ( RAtomPosition(ind,1) .EQ. MNP(jnd,1) .AND. &
                RAtomPosition(ind,2) .EQ. MNP(jnd,2) .AND. &
                RAtomPosition(ind,3) .EQ. MNP(jnd,3) .AND. &
                SAtomName(ind) .EQ. SMNP(jnd) ) THEN
              !this seems NOT a unique coordinate
              Lunique=.FALSE.
              EXIT
           ENDIF
           Lunique=.TRUE.
        ENDDO
        IF(Lunique .EQV. .TRUE.) THEN
           Iuniind=Iuniind+1
           MNP(Iuniind,:)= RAtomPosition(ind,:)
           SMNP(Iuniind)= SAtomName(ind)
        ENDIF
     ENDDO
     
     IMaxPossibleNAtomsUnitCell = Iuniind
     
  END IF
  
  CALL Message("CountTotalAtoms",IMoreInfo,IErr,MessageVariable = "IMaxPossibleNAtomsUnitCell",IVariable=IMaxPossibleNAtomsUnitCell)

  DEALLOCATE(MNP,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(",my_rank,")error",IErr,"deallocating MNP"
     RETURN
  ENDIF
  DEALLOCATE(SMNP,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(",my_rank,")error",IErr,"deallocating SMNP"
     RETURN
  END IF
  DEALLOCATE(RAtomPosition,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(",my_rank,")error",IErr,"deallocating RAtomPosition"
     RETURN
  END IF
  DEALLOCATE(SAtomName,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(",my_rank,")error",IErr,"deallocating SAtomName"
     RETURN
  ENDIF
  
END SUBROUTINE CountTotalAtoms

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE GreatestCommonDivisor(ITotalProcesses,INooDWFs,ISubgroups)

USE MyNumbers

INTEGER(IKIND) :: a,b,c
INTEGER(IKIND), INTENT(IN) :: ITotalProcesses,INooDWFs
INTEGER(IKIND), INTENT(OUT) :: ISubgroups

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Defines a Lorentzian Distribution for any parameter input
FUNCTION Lorentzian(FWHM,x,x_0,offset)

  USE MyNumbers
  
  USE RPara; 

  IMPLICIT NONE

  REAL(RKIND):: FWHM,x,x_0,offset,LORENTZIAN

  LORENTZIAN = FWHM/(((x+x_0)**2)+offset)
  
END FUNCTION Lorentzian

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Defines a Gaussian distribution for any parameter input 
FUNCTION Gaussian(height,x,peakcentre,standarddeviation,intercept)

  USE MyNumbers

  IMPLICIT NONE

  REAL(RKIND):: height,x,peakcentre,standarddeviation,intercept,gaussian

  Gaussian = height*exp(-(((x-peakcentre)**2)/(2*(standarddeviation**2))))+ intercept
  
END FUNCTION Gaussian
