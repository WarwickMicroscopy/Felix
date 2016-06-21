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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Defines a Kirkland scattering factor 
FUNCTION Kirkland(Iz,Rq)
!From Appendix C of Kirkland, "Advanced Computing in Electron Microscopy", 2nd ed.
!z is atomic number, q is magnitude of scattering vector in Angstroms (NB exp(2*pi*i*q.r), optical convention)
  USE MyNumbers
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: Iz,ind
  REAL(RKIND):: Kirkland,Ra,Rb,Rc,Rd,Rq

  Kirkland=ZERO;
  !Equation C.15
  DO ind = 1,3
    Ra=RScattFactors(Iz,ind*2-1);
    Rb=RScattFactors(Iz,ind*2);
    Rc=RScattFactors(Iz,ind*2+5);
    Rd=RScattFactors(Iz,ind*2+6);
    Kirkland = Kirkland + Ra/((Rq**2)+Rb)+Rc*EXP(-(Rd*Rq**2));
  END DO
  
END FUNCTION Kirkland
