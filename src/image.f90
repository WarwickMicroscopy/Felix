!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14/15/16, all rights reserved
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

SUBROUTINE ImageInitialisation( IErr )

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$%
!!$%     Determines Montage size as twice the distance to the furtherest pixel + 1
!!$%
!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  USE WriteToScreen

  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  REAL(RKIND) :: DummyConvergenceAngle
  INTEGER(IKIND) :: IErr,ind,jnd

  CALL Message("ImageInitialisation",IMust,IErr)
  
  !Determine Positions of reflections in final image (may not need to be here)
  CALL Message("ImageInitialisation",IInfo,IErr, MessageVariable = "nReflections", IVariable = nReflections)
  CALL Message("ImageInitialisation",IInfo,IErr, MessageVariable = "RMinimumGMag", RVariable = RMinimumGMag)
  PRINT*,"RB boo ImageInitialisation1"

  ! positions of the centres of the disks
  DO ind=1,nReflections
     RhklPositions(ind,1) = RgPoolT(ind,1)/RMinimumGMag
     RhklPositions(ind,2) = RgPoolT(ind,2)/RMinimumGMag
  ENDDO
  PRINT*,"RB ImageInitialisation"
  PRINT*,"RB ImageInitialisation IOutputReflections",IOutputReflections
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

  DO ind=1,2
     CALL Message("ImageInitialisation",IInfo,IErr, &
          MessageVariable = "IImageSizeXY(ind)", IVariable = IImageSizeXY(ind))
  END DO
  
  RETURN

END SUBROUTINE ImageInitialisation

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE MontageInitialisation(IPixelHorizontalPosition,IPixelVerticalPosition,&
     IThicknessindex,RMontageImage,RIntensityValues,IErr)

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$%
!!$%      Places Calculated pixels into montage 1 pixel, per reflection per call
!!$%
!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE WriteToScreen
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ::&
       IThicknessindex, IMontagePixelVerticalPosition, IMontagePixelHorizontalPosition,&
       hnd,IPixelHorizontalPosition,IPixelVerticalPosition,Ierr
  REAL(RKIND), DIMENSION(MAXVAL(IImageSizeXY),MAXVAL(IImageSizeXY),IThicknessCount),INTENT(OUT) :: &
       RMontageImage
  REAL(RKIND), DIMENSION(INoOfLacbedPatterns) :: RIntensityValues

!!$  Only print out once when first entered - use message counter

  IF (my_rank.EQ.0) THEN
     DO WHILE (IMessageCounter .LT.1)
        CALL Message("MontageInitialisation",IMust,IErr)
        CALL Message("MontageInitialisation",IMust+IDebug,IErr,MessageString="Is looping (called from MontageSetup)")
        IMessageCounter = IMessageCounter +1
     END DO
  END IF
 
 DO hnd = 1,INoOfLacbedPatterns

     IF(IHKLSelectFLAG.EQ.0) THEN
        
        IF (RConvergenceAngle.LT.ONE) THEN
           IMontagePixelVerticalPosition = &
                NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(hnd,2))
           IMontagePixelHorizontalPosition = &
                NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(hnd,1))
        ELSE
        
!!$           If the Convergence angle is > 1 causing disk overlap in experimental pattern, 
!!$           then plot as if convergence angle was 0.95 (non-physical but makes a pretty picture)
           IMontagePixelVerticalPosition = &
                NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/0.95D0)*RhklPositions(hnd,2))
           IMontagePixelHorizontalPosition = &
                NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/0.95D0)*RhklPositions(hnd,1))
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
        
     ELSE  
      
        IF (RConvergenceAngle.LT.ONE) THEN
           IMontagePixelVerticalPosition = &
                NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(IOutputReflections(hnd),2))
           IMontagePixelHorizontalPosition = &
                NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(IOutputReflections(hnd),1))
        ELSE
!!$           If the Convergence angle is > 1 causing disk overlap in experimental pattern, 
!!$           then plot as if convergence angle was 0.95 (non-physical but makes a pretty picture)
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

!---------------------------------------------------------------------
!
SUBROUTINE ImageMaskInitialisation (IErr)

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$%
!!$%    Creates a circular or square image mask depending on the value of IMaskFLAG
!!$%       and assigns pixel locations for each one for MPI load balancing
!!$%
!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  USE MyNumbers
  USE WriteToScreen

  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,jnd, ierr,InnerRadiusFLAG
  REAL(RKIND) :: Rradius, RImageRadius
  
  CALL Message("ImageMaskInitialisation",IMust,IErr)
  PRINT*,"RB boo2"

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

  IF (RInnerConvergenceAngle.GT.ZERO) THEN
     DO ind=1,2*IPixelCount
        DO jnd=1,2*IPixelCount
           Rradius= (ind-(REAL(IPixelCount,RKIND)+0.5))**2 + &
                (jnd-(REAL(IPixelCount,RKIND)+0.5))**2
           Rradius=SQRT(DBLE(Rradius))
           RImageRadius = (IPixelCount+0.5)*(RInnerConvergenceAngle/RConvergenceAngle)
           IF(Rradius.LE.RImageRadius) THEN
              RMask(jnd,ind) = 0
              IPixelTotal= IPixelTotal - 1
           END IF
        ENDDO
     ENDDO
  END IF

  ALLOCATE(IPixelLocations(IPixelTotal,2), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ImagemaskInitialization(", my_rank, ") error ", IErr, " in ALLOCATE of IPixelLocations"
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

INTEGER(IKIND) FUNCTION CountPixels(IErr)
  
!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$%
!!$%     Counts pixels in requested image for memory allocation
!!$%
!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       ind,jnd, IErr
  REAL(RKIND) :: &
       Rradius, RImageRadius
  
  IF((IWriteFLAG.EQ.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"CountPixels (",my_rank,")"
  END IF
  
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
