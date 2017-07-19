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

! All procedures conatained in this file:
! SortHKL()
! GreatestCommonDivisor()
! Lorentzian()
! Gaussian()


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
  USE WriteToScreen
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE WriteToScreen

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
3       CONTINUE
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
100        ENDDO
           I=I-M
           IF(I.GE.1) GOTO 3
        ENDIF
11   ENDDO
12 ENDDO
  
  RETURN

END SUBROUTINE SortHKL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!>
!! Procedure-description:
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
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
