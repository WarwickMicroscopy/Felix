!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-17, all rights reserved
!
! Version: 2.0
! Date: 31-08-2022
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build: cRED
! Author:  r.beanland@warwick.ac.uk
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

!>
!! Module-description: 
!!
MODULE utilities_mod

  IMPLICIT NONE

  CONTAINS
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

