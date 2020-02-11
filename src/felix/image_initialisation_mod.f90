!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: :VERSION: RB_coord / 1.14 /
! Date:    :DATE: 15-01-2019
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build:   :BUILD: Mode F: test different lattice types" 
! Author:  :AUTHOR: r.beanland
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
    USE IPARA, ONLY : IMask, IPixelTotal, IPixelLocations

    ! global input
    USE IPARA, ONLY : IPixelCount, IMaskFLAG
    
    IMPLICIT NONE
    
    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND) :: ind, jnd, InnerRadiusFLAG
    REAL(RKIND) :: Rradius, RImageRadius

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
            IMask(jnd,ind) = 1
            IPixelTotal= IPixelTotal + 1			  
          ELSE
            IMask(jnd,ind) = 0
          END IF
        ENDDO
      ENDDO
    CASE(1) ! square
      IMask = 1
      IPixelTotal = (2*IPixelCount)**2
    END SELECT

    ! Removed InnerConvergenceAngle here

    ALLOCATE(IPixelLocations(IPixelTotal,2),STAT=IErr)
    IF(l_alert(IErr,"ImageMaskInitialisation","allocate IMask")) RETURN

    IPixelTotal = 0

    SELECT CASE (IMaskFLAG)
    CASE(0) ! circle
      DO ind=1,2*IPixelCount
        DO jnd=1,2*IPixelCount
          IF(IMask(ind,jnd).GT.ZERO) THEN
            IPixelTotal= IPixelTotal + 1
            IPixelLocations(IPixelTotal,1) = ind
            IPixelLocations(IPixelTotal,2) = jnd
          ENDIF
        ENDDO
      ENDDO
    CASE(1) ! square
      DO ind = 1,2*IPixelCount
        DO jnd = 1,2*IPixelCount
          IF(IMask(ind,jnd).GT.ZERO) THEN
            IPixelTotal = IPixelTotal+1
            IPixelLocations(IPixelTotal,1) = ind
            IPixelLocations(IPixelTotal,2) = jnd
          END IF
        END DO
      END DO
    END SELECT
   
  END SUBROUTINE ImageMaskInitialisation

END MODULE image_initialisation_mod
