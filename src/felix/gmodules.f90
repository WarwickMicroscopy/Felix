!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: 1.3
! Date: 13-05-2024
! Time:
! Status:  
! Build: g-vector limit   
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
!! Module-description: Some numeric computational and mathematical paramter
!! constants as well as maths functions - INIT_NUMBERS, ARG, CROSS, DOT
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
MODULE MyNumbers     
  IMPLICIT NONE

  INTEGER, PARAMETER :: IKIND = KIND(1)
  INTEGER, PARAMETER :: RKIND = KIND(1.0D0)
  INTEGER, PARAMETER :: CKIND = RKIND 
  INTEGER(IKIND), PARAMETER :: INElements=110
  INTEGER(IKIND), PARAMETER :: ITHREE=3_IKIND
  
  REAL(RKIND) :: PI,TWOPI,FOURPI,ONEPLS,ONEMNS,SQRTHALF,SQRTTWO,DEG2RADIAN,TWODEG2RADIAN,RADIAN2DEG
  REAL(RKIND), PARAMETER :: ZERO=0.0_RKIND,ONE=1.0_RKIND,TWO=2.0_RKIND,THREE=3.0_RKIND, &
        FOUR=4.0_RKIND,SIX=6.0_RKIND,TEN=10.0_RKIND,HUNDRED=100.0_RKIND,THOUSAND=1000.0_RKIND
  REAL(RKIND), PARAMETER :: HALF = 0.5_RKIND, QUARTER = 0.25_RKIND, EIGHTH = 0.125_RKIND, &
       THIRD=0.3333333333333333_RKIND, TWOTHIRD=0.6666666666666_RKIND,  &
       NEGTHIRD=-0.3333333333333333_RKIND, NEGTWOTHIRD=-0.6666666666666_RKIND
  REAL(RKIND) :: TINY= 1.0D-9,NEGHUGE=-1.0D9
  COMPLEX(RKIND), PARAMETER :: CZERO = (0.0_RKIND,0.0_RKIND), CONE = (1.0_RKIND,0.0_RKIND), &
       CIMAGONE= (0.0_RKIND,1.0_RKIND)            
  
  REAL(RKIND) :: RKiloByte,RMegaByte,RGigaByte,RTeraByte 
  
CONTAINS
  SUBROUTINE INIT_NUMBERS
    PI       = 4.0_RKIND* ATAN(1.0_RKIND)
    TWOPI    = 8.0_RKIND* ATAN(1.0_RKIND)
    FOURPI   = 16.0_RKIND* ATAN(1.0_RKIND)
    ONEMNS   = SQRT(EPSILON(ONEMNS))
    ONEPLS   = ONE + ONEMNS
    ONEMNS   = ONE - ONEMNS
    SQRTHALF = DSQRT(0.5_RKIND)
    SQRTTWO  = DSQRT(2.0_RKIND)
    RKiloByte = 2.0_RKIND**10.0_RKIND
    RMegaByte = 2.0_RKIND**20.0_RKIND
    RGigaByte = 2.0_RKIND**30.0_RKIND
    RTeraByte = 2.0_RKIND**40.0_RKIND
    DEG2RADIAN= PI/180.D0
    TWODEG2RADIAN=TWOPI/180.D0
    RADIAN2DEG=180.D0/PI
  END SUBROUTINE INIT_NUMBERS

  FUNCTION ARG(X,Y)
    
    REAL(KIND=RKIND) ARG, X, Y
    
    IF( X > 0. ) THEN 
       ARG= ATAN(Y/X)
    ELSE IF ( (X == 0.) .and. (Y > 0. )) THEN 
       ARG = PI/2.0_RKIND
    ELSE IF ( (X == 0.) .and. (Y < 0. )) THEN 
       ARG = -PI/2.0_RKIND
    ELSE IF ( (X < 0. ) .and. (Y >= 0.)) THEN 
       ARG = PI + ATAN(Y/X)
    ELSE IF ( (X < 0. ) .and. (Y < 0. )) THEN 
       ARG = -PI + ATAN(Y/X)
    ELSE IF ( (X == 0.0) .and. (Y == 0.)) THEN
       PRINT*, "ARG: both X and Y ==0, undefined --- using ARG=0"
       ARG=0.0_RKIND
    ENDIF
    
    RETURN
  END FUNCTION ARG

  FUNCTION CROSS(a, b)
    REAL(RKIND), DIMENSION(3) :: CROSS
    REAL(RKIND), DIMENSION(3), INTENT(IN) :: a, b
    
    CROSS(1) = a(2) * b(3) - a(3) * b(2)
    CROSS(2) = a(3) * b(1) - a(1) * b(3)
    CROSS(3) = a(1) * b(2) - a(2) * b(1)
  END FUNCTION CROSS

  FUNCTION DOT(a, b)
    REAL(RKIND) :: DOT
    REAL(RKIND), DIMENSION(3), INTENT(IN) :: a, b
 
    DOT= a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
  END FUNCTION DOT
  
END MODULE MyNumbers


!>
!! Module-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
MODULE MyMPI !?? grants MPI, so use of both in the code is unnecessary

  USE MPI
  USE MyNumbers, ONLY : IKIND

  INTEGER(IKIND) :: my_rank, p, srce, dest
  INTEGER, PARAMETER :: root = 0
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status_info
  
END MODULE MyMPI

!>
!! Module-description:
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
MODULE MyFFTW   

  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  INTEGER, PARAMETER :: C_FFTW_R2R_KIND = C_INT32_T
  INCLUDE  'fftw3.f03'
  
END MODULE MyFFTW
