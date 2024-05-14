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
! Time:    :TIME:
! Status:  :RLSTATUS:
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
!! Module-description: 
!!
MODULE image_initialisation_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ImageMaskInitialisation

  CONTAINS

  !>
  !! Procedure-description: Creates a circular or square image mask depending on
  !! the value of IMaskFLAG and assigns pixel locations for each one for MPI load
  !! balancing
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE ImageMaskInitialisation (IErr)
    
    ! this procedure is called once in felixrefine during initial setup
    USE MyNumbers
    USE message_mod

    ! global output
    USE IPARA, ONLY : IPixelTotal, IPixelLocations

    ! global input
    USE IPARA, ONLY : IPixelCount
    
    IMPLICIT NONE
    
    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND) :: ind, jnd, InnerRadiusFLAG
    REAL(RKIND) :: Rradius, RImageRadius

    IPixelTotal =0
    ! SELECT CASE (IMaskFLAG)
    ! CASE(0) ! circle
      ! DO ind=1,2*IPixelCount
        ! DO jnd=1,2*IPixelCount
            ! ! IMask(jnd,ind) = 1
            ! IPixelTotal= IPixelTotal + 1			  
        ! ENDDO
      ! ENDDO
    ! CASE(1) ! square
      ! IMask = 1
      IPixelTotal = (2*IPixelCount)**2

    ! Removed InnerConvergenceAngle here

    ALLOCATE(IPixelLocations(IPixelTotal,2),STAT=IErr)
    IF(l_alert(IErr,"ImageMaskInitialisation","allocate IMask")) RETURN

    IPixelTotal = 0

    ! SELECT CASE (IMaskFLAG)
    ! CASE(0) ! circle
      ! DO ind=1,2*IPixelCount
        ! DO jnd=1,2*IPixelCount
            ! IPixelTotal= IPixelTotal + 1
            ! IPixelLocations(IPixelTotal,1) = ind
            ! IPixelLocations(IPixelTotal,2) = jnd
        ! ENDDO
      ! ENDDO
    ! CASE(1) ! square
      DO ind = 1,2*IPixelCount
        DO jnd = 1,2*IPixelCount
            IPixelTotal = IPixelTotal+1
            IPixelLocations(IPixelTotal,1) = ind
            IPixelLocations(IPixelTotal,2) = jnd
        END DO
      END DO
   
  END SUBROUTINE ImageMaskInitialisation

END MODULE image_initialisation_mod
