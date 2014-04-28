!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! BlochSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
!
! (C) 2013/14, all right reserved
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This file is part of BlochSim.
!
!  BlochSim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  BlochSim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with BlochSim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: gmodules.f90,v 1.11 2014/03/25 15:37:30 phsht Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!$Log: gmodules.f90,v $
!Revision 1.11  2014/03/25 15:37:30  phsht
!more on GPL
!
!Revision 1.10  2014/03/25 15:35:34  phsht
!included consistent start of each source file and GPL statements
!
!Revision 1.9  2014/03/21 15:55:36  phslaz
!New Lacbed code Working
!
!Revision 1.7  2013/12/19 16:30:27  phsht
!new version of HKLMake(), using .cif information
!
!Revision 1.6  2013/11/26 18:27:18  phsht
!MPI file IO now working for WI (intensities)
!
!Revision 1.5  2013/11/25 18:26:33  phsht
!progress on the MPI version
!
!Revision 1.4  2013/06/11 14:53:08  phsht
!more work in the Diffraction defs part
!
!Revision 1.3  2013/06/10 08:20:28  phsht
!more work before realizing that MATLAB file is not quite the latest version
!
!Revision 1.2  2013/06/07 07:11:27  phsht
!more ground work
!
!Revision 1.1  2013/04/03 19:35:45  phsht
!first installation of basic Fortran routines/structure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE MyNumbers     
  IMPLICIT NONE

  INTEGER, PARAMETER :: IKIND = KIND(1)
  INTEGER, PARAMETER :: RKIND = KIND(1.0D0)
  INTEGER, PARAMETER :: CKIND = RKIND

  REAL(KIND=RKIND) :: PI, TWOPI, ONEPLS, ONEMNS, &
       SQRTHALF, SQRTTWO

  REAL(KIND=RKIND), PARAMETER :: ZERO = 0.0, ONE = 1.0 ,TWO = 2.0, &
       THREE = 3.0, FOUR = 4.0
  COMPLEX(KIND=RKIND), PARAMETER :: CZERO = (0.0d0,0.0d0), CONE = (1.0d0,0.0d0), &
       CIMAGONE= (0.0d0,1.0d0)            

  REAL(KIND=RKIND), PARAMETER :: HALF = 0.5D0, QUARTER = 0.25D0, EIGHTH = 0.125D0, &
       THIRD=0.3333333333333333D0, TWOTHIRD=0.6666666666666D0, &
       NEGTHIRD=-0.3333333333333333D0, NEGTWOTHIRD=-0.6666666666666D0

  REAL(KIND=RKIND) :: TINY= 1.0D-9

  REAL(RKIND) :: RKiloByte,RMegaByte,RGigaByte,RTeraByte 
  
CONTAINS
  SUBROUTINE INIT_NUMBERS
    PI       = 4.0D0* ATAN(1.0D0)
    TWOPI    = 8.0D0* ATAN(1.0D0)
    ONEMNS   = SQRT(EPSILON(ONEMNS))
    ONEPLS   = ONE + ONEMNS
    ONEMNS   = ONE - ONEMNS
    SQRTHALF = DSQRT(0.5D0)
    SQRTTWO  = DSQRT(2.0D0)
    RKiloByte = 2.0D0**10.0D0
    RMegaByte = 2.0D0**20.0D0
    RGigaByte = 2.0D0**30.0D0
    RTeraByte = 2.0D0**40.0D0
    
  END SUBROUTINE INIT_NUMBERS

  FUNCTION ARG(X,Y)
    
    REAL(KIND=RKIND) ARG, X, Y
    
    IF( X > 0. ) THEN 
       ARG= ATAN(Y/X)
    ELSE IF ( (X == 0.) .and. (Y > 0. )) THEN 
       ARG = PI/2.0D0
    ELSE IF ( (X == 0.) .and. (Y < 0. )) THEN 
       ARG = -PI/2.0D0
    ELSE IF ( (X < 0. ) .and. (Y >= 0.)) THEN 
       ARG = PI + ATAN(Y/X)
    ELSE IF ( (X < 0. ) .and. (Y < 0. )) THEN 
       ARG = -PI + ATAN(Y/X)
    ELSE IF ( (X == 0.0) .and. (Y == 0.)) THEN
       PRINT*, "ARG(): both X and Y ==0, undefined --- using ARG=0"
       ARG=0.0D0
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

MODULE MyMPI

  USE MPI
  USE MyNumbers

  INTEGER(IKIND) :: my_rank, p, srce, dest
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status_info
  
END MODULE MyMPI



