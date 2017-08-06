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
! ImageSetup( )
! ImageInitialisation( )
! MontageInitialisation( )
! ImageMaskInitialisation( )
! CountPixels( )

!>
!! Procedure-description: Image setup, including defining image masks
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE ImageSetup (IErr) 

  USE MyNumbers

  USE IPara; USE RPara

  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr

  CALL ImageInitialisation( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"Error:ImageSetup(",my_rank,")error in ImageInitialistion"
     RETURN
  END IF
 
  !--------------------------------------------------------------------
  ! define image masks
  CALL ImageMaskInitialisation(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Error:ImageSetup(",my_rank,")error in ImageMaskInitialisation"
     RETURN
  END IF

END SUBROUTINE ImageSetup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description: Sets the positions of the centres of the disks
!! and calculate the size of the final image
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE ImageInitialisation( IErr )

  USE MyNumbers

  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  REAL(RKIND) :: DummyConvergenceAngle
  INTEGER(IKIND) :: IErr,ind,jnd

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
!! Procedure-description: Places Calculated pixels into montage 1 pixel, per reflection per call
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE MontageInitialisation(IPixelHorizontalPosition,IPixelVerticalPosition,&
     IThicknessindex,RMontageImage,RIntensityValues,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ::IThicknessindex, IMontagePixelVerticalPosition, IMontagePixelHorizontalPosition,&
       hnd,IPixelHorizontalPosition,IPixelVerticalPosition,Ierr
  REAL(RKIND), DIMENSION(MAXVAL(IImageSizeXY),MAXVAL(IImageSizeXY),IThicknessCount),INTENT(OUT) :: RMontageImage
  REAL(RKIND), DIMENSION(INoOfLacbedPatterns) :: RIntensityValues
 
 DO hnd = 1,INoOfLacbedPatterns
   IF(IHKLSelectFLAG.EQ.0) THEN
     IF (RConvergenceAngle.LT.ONE) THEN
       IMontagePixelVerticalPosition = &
           NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(hnd,2))
       IMontagePixelHorizontalPosition = &
           NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(hnd,1))
     ELSE
       ! If the Convergence angle is > 1 causing disk overlap in experimental pattern, 
       ! then plot as if convergence angle was 0.95 (non-physical but makes a pretty picture)
       IMontagePixelVerticalPosition = &
                NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/0.95D0)*RhklPositions(hnd,2))
       IMontagePixelHorizontalPosition = &
                NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/0.95D0)*RhklPositions(hnd,1))
     END IF
     RMontageImage(IMontagePixelVerticalPosition-IPixelCount+IPixelVerticalPosition,&
             IMontagePixelHorizontalPosition-IPixelCount+IPixelHorizontalPosition,&
             IThicknessIndex) = &
             RMontageImage(IMontagePixelVerticalPosition-IPixelCount+IPixelVerticalPosition,&
             IMontagePixelHorizontalPosition-IPixelCount+IPixelHorizontalPosition,&
             IThicknessIndex) + RIntensityValues(hnd)
     ELSE  
       IF (RConvergenceAngle.LT.ONE) THEN
         IMontagePixelVerticalPosition = &
                NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(IOutputReflections(hnd),2))
         IMontagePixelHorizontalPosition = &
                NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(IOutputReflections(hnd),1))
       ELSE
         ! If the Convergence angle is > 1 causing disk overlap in experimental pattern, 
         ! then plot as if convergence angle was 0.95 (non-physical but makes a pretty picture)
         IMontagePixelVerticalPosition = &
                NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/0.95D0)*RhklPositions(IOutputReflections(hnd),2))
         IMontagePixelHorizontalPosition = &
                NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/0.95D0)*RhklPositions(IOutputReflections(hnd),1))
       END IF
       RMontageImage(&
             IMontagePixelVerticalPosition-IPixelCount+IPixelVerticalPosition,&
             IMontagePixelHorizontalPosition-IPixelCount+IPixelHorizontalPosition,&
             IThicknessIndex) = &
             RMontageImage(&
             IMontagePixelVerticalPosition-IPixelCount+IPixelVerticalPosition,&
             IMontagePixelHorizontalPosition-IPixelCount+IPixelHorizontalPosition,&
             IThicknessIndex) + &
             RIntensityValues(hnd)
    END IF
  END DO

END SUBROUTINE MontageInitialisation

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!>
!! Procedure-description: Creates a circular or square image mask depending on
!! the value of IMaskFLAG and assigns pixel locations for each one for MPI load
!! balancing
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE ImageMaskInitialisation (IErr)
  
  USE MyNumbers

  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,jnd, ierr,InnerRadiusFLAG
  REAL(RKIND) :: Rradius, RImageRadius

  ALLOCATE(RMask(2*IPixelCount,2*IPixelCount),STAT=IErr)
  IF( IErr.NE.0 ) THEN
    PRINT*,"Error:ImageMaskInitialisation(",my_rank,")error allocating RMask"
    RETURN
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
           END IF
        ENDDO
     ENDDO
  CASE(1) ! square
     RMask = 1
     IPixelTotal = (2*IPixelCount)**2
  END SELECT

  !Removed InnerConvergenceAngle here

  ALLOCATE(IPixelLocations(IPixelTotal,2),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Error:ImagemaskInitialization(",my_rank,")error allocating IPixelLocations"
     RETURN
  ENDIF

  IPixelTotal = 0

  SELECT CASE (IMaskFLAG)
  CASE(0) ! circle
     DO ind=1,2*IPixelCount
        DO jnd=1,2*IPixelCount
           IF(RMask(ind,jnd).GT.ZERO) THEN
              IPixelTotal= IPixelTotal + 1
              IPixelLocations(IPixelTotal,1) = ind
              IPixelLocations(IPixelTotal,2) = jnd
           ENDIF
        ENDDO
     ENDDO
  CASE(1) ! square
     IPixelTotal = 0
     DO ind = 1,2*IPixelCount
        DO jnd = 1,2*IPixelCount
           IF(RMask(ind,jnd).GT.ZERO) THEN
              IPixelTotal = IPixelTotal+1
              IPixelLocations(IPixelTotal,1) = ind
              IPixelLocations(IPixelTotal,2) = jnd
           END IF
        END DO
     END DO
  END SELECT
 
END SUBROUTINE ImageMaskInitialisation

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description: Counts pixels in requested image for memory allocation
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
INTEGER(IKIND) FUNCTION CountPixels(IErr)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE terminal_output
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,jnd, IErr
  REAL(RKIND) :: Rradius, RImageRadius
  
  CALL message ( LL, "Counting Pixels")
  
  CountPixels =0

  SELECT CASE (IMaskFLAG)

  CASE(0) ! circle

     DO ind=1,2*IPixelCount
        DO jnd=1,2*IPixelCount
           Rradius= (ind-(REAL(IPixelCount,RKIND)+0.5))**2 + &
                (jnd-(REAL(IPixelCount,RKIND)+0.5))**2
           Rradius=SQRT(DBLE(Rradius))
           RImageRadius = IPixelCount+0.5
           IF(Rradius.LE.RImageRadius) THEN
              CountPixels =  CountPixels + 1
           ENDIF
        ENDDO
     ENDDO

  CASE(1) ! square

     CountPixels = (2*IPixelCount)**2

  END SELECT
END FUNCTION CountPixels
