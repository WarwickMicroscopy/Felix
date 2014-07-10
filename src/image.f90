!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! FelixSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
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
!  This file is part of FelixSim.
!
!  FelixSim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  FelixSim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with FelixSim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ImageInitialization( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  REAL(RKIND) dummyCA

  INTEGER IErr, ind
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ImageInitialization()"
  END IF
  
!Determine Positions of reflections in final image (may not need to be here)

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ImageInitialization(",my_rank,") nReflections,MinGMag =", nReflections, RMinimumGMag
  END IF

  ! positions of the centres of the disks

  DO ind=1,nReflections
     Rhklpositions(ind,1) = &
          RgVecMatT(ind,1)/RMinimumGMag
     Rhklpositions(ind,2) = &
          RgVecMatT(ind,2)/RMinimumGMag
  ENDDO
  
  ! size of final image
  
  IF(RConvergenceAngle .LT. ONE) THEN
     dummyCA=RConvergenceAngle
  ELSE
     dummyCA=0.95D0
  ENDIF
  
DO ind=1,SIZE(Rhklpositions,DIM=2)
     IImageSizeXY(ind)= CEILING(&
          4.0D0*REAL(IPixelCount,RKIND)/dummyCA * &
          (MAXVAL(ABS(Rhklpositions(1:IReflectOut,ind)))+1.0D0) )
  ENDDO
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ImageInitialization(",my_rank,") IImageSizeXY=", IImageSizeXY
  END IF
    
  RETURN

END SUBROUTINE ImageInitialization

SUBROUTINE MakeMontagePixel(ind,jnd,ithicknessindex,RMontageImage,RIntensity,Ierr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) IThicknessindex, IXpos, IYpos,hnd,ind,jnd,knd,Ierr
  REAL(RKIND), DIMENSION(MAXVAL(IImageSizeXY),MAXVAL(IImageSizeXY),IThicknessCount),INTENT(OUT) :: &
       RMontageImage
  REAL(RKIND), DIMENSION(IReflectOut) :: RIntensity
 
  DO hnd = 1,IReflectOut
     
     IF(IHKLSelectFLAG.EQ.0) THEN
        
        IF (RConvergenceAngle.LT.ONE) THEN
           IXpos = NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(hnd,2))
           IYpos = NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(hnd,1))
        ELSE
           
           !If the Convergence angle is > 1 causing disk overlap in experimental pattern, then plot as if convergence angle was 0.95 (non-physical but makes a pretty picture)
           IXpos = NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/0.95D0)*RhklPositions(hnd,2))
           IYpos = NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/0.95D0)*RhklPositions(hnd,1))
        ENDIF
        !PRINT*,"IXpos,IYpos,IPixelCount,jnd,ind = ",IXpos,IYpos,IPixelCount,jnd,ind
        RMontageImage(IXpos-IPixelCount+jnd,IYpos-IPixelCount+ind,IThicknessIndex) = &
             RMontageImage(IXpos-IPixelCount+jnd,IYpos-IPixelCount+ind,IThicknessIndex) + &
             RIntensity(hnd)
        
     ELSE
        
        IF (RConvergenceAngle.LT.ONE) THEN
           IXpos = NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(IOutputReflections(hnd),2))
           IYpos = NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(IOutputReflections(hnd),1))
        ELSE
           
           !If the Convergence angle is > 1 causing disk overlap in experimental pattern, then plot as if convergence angle was 0.95 (non-physical but makes a pretty picture)
           IXpos = NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/0.95D0)*RhklPositions(IOutputReflections(hnd),2))
           IYpos = NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/0.95D0)*RhklPositions(IOutputReflections(hnd),1))
        ENDIF
        !PRINT*,"IXpos,IYpos,IPixelCount,jnd,ind = ",IXpos,IYpos,IPixelCount,jnd,ind
        RMontageImage(IXpos-IPixelCount+jnd,IYpos-IPixelCount+ind,IThicknessIndex) = &
             RMontageImage(IXpos-IPixelCount+jnd,IYpos-IPixelCount+ind,IThicknessIndex) + &
             RIntensity(hnd)
     END IF
  END DO

END SUBROUTINE MakeMontagePixel
