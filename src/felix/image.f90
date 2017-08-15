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

!>
!! Module-description: 
!!
MODULE image

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ImageInitialisation, ImageMaskInitialisation

  CONTAINS

  !>
  !! Procedure-description: Sets the positions of the centres of the disks
  !! and calculate the size of the final image
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE ImageInitialisation( IErr )

    !?? called once in felixrefine setup

    USE MyNumbers
    USE message_mod

    ! global outputs
    USE RPARA, ONLY : RhklPositions
    USE IPARA, ONLY : IImageSizeXY
    
    ! global inputs
    USE RPARA, ONLY : RgPool, RMinimumGMag, RConvergenceAngle
    USE IPARA, ONLY : IPixelCount, INoOfLacbedPatterns, nReflections, IHKLSelectFLAG, &
                      IOutputReflections
    
    IMPLICIT NONE

    INTEGER(IKIND), INTENT(OUT) :: IErr
    REAL(RKIND) :: DummyConvergenceAngle
    INTEGER(IKIND) :: ind, jnd

    ! positions of the centres of the disks
    DO ind=1,nReflections
      RhklPositions(ind,1) = RgPool(ind,1)/RMinimumGMag
      RhklPositions(ind,2) = RgPool(ind,2)/RMinimumGMag
    ENDDO

    ! size of final image
    IF(RConvergenceAngle .LT. ONE) THEN
      DummyConvergenceAngle=RConvergenceAngle
    ELSE
      DummyConvergenceAngle=0.95_RKIND
    END IF
    IF(IHKLSelectFLAG.EQ.0) THEN
      DO ind=1,2
        IImageSizeXY(ind)= CEILING(&
             FOUR*REAL(IPixelCount,RKIND)/DummyConvergenceAngle * &
            (MAXVAL(ABS(RhklPositions(1:INoOfLacbedPatterns,ind)))+ONE) )
      ENDDO
    ELSE
      DO ind=1,2
        DO jnd = 1,INoOfLacbedPatterns
          IImageSizeXY(ind)= CEILING(&
             FOUR*REAL(IPixelCount,RKIND)/DummyConvergenceAngle * &
            (MAXVAL(ABS(RhklPositions(IOutputReflections(1:INoOfLacbedPatterns),ind)))+ONE) )
        END DO
      ENDDO
    END IF
    
    RETURN

  END SUBROUTINE ImageInitialisation

  !!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  !>
  !! Procedure-description: Creates a circular or square image mask depending on
  !! the value of IMaskFLAG and assigns pixel locations for each one for MPI load
  !! balancing
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE ImageMaskInitialisation (IErr)
    
    !?? called once in felixrefine setup

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

END MODULE image
