!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
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

SUBROUTINE ImageInitialization( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  REAL(RKIND) dummyCA

  !changed - was missing (IKIND)
  INTEGER(IKIND) IErr, ind,jnd
  
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
  IF(IHKLSelectFLAG.EQ.0) THEN
     DO ind=1,SIZE(Rhklpositions,DIM=2)
        IImageSizeXY(ind)= CEILING(&
             4.0D0*REAL(IPixelCount,RKIND)/dummyCA * &
             (MAXVAL(ABS(Rhklpositions(1:IReflectOut,ind)))+1.0D0) )
     ENDDO
  ELSE
     DO ind=1,SIZE(Rhklpositions,DIM=2)
        DO jnd = 1,IReflectOut
           IImageSizeXY(ind)= CEILING(&
                4.0D0*REAL(IPixelCount,RKIND)/dummyCA * &
                (MAXVAL(ABS(Rhklpositions(IOutputReflections(1:IReflectOut),ind)))+1.0D0) )
        END DO
     ENDDO
  END IF
     
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

!---------------------------------------------------------------------
!
SUBROUTINE ImageMaskInitialization (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER ind,jnd, ierr
  REAL(RKIND) :: Rradius, RImageRadius
  
  IF((IWriteFLAG.EQ.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: ImageMaskInitialization()"
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
           ENDIF
        ENDDO
     ENDDO
  CASE(1) ! square
     RMask = 1
     IPixelTotal = (2*IPixelCount)**2
  END SELECT
  
  ALLOCATE( &
       IPixelLocations(IPixelTotal,2), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ImagemaskInitialization(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  IPixelTotal = 0
 
  SELECT CASE (IMaskFLAG)
  CASE(0) ! circle
     DO ind=1,2*IPixelCount
        DO jnd=1,2*IPixelCount
           Rradius= (ind-(REAL(IPixelCount,RKIND)+0.5))**2 + &
                (jnd-(REAL(IPixelCount,RKIND)+0.5))**2
           Rradius=SQRT(DBLE(Rradius))
           RImageRadius = IPixelCount+0.5
           IF(Rradius.LE.RImageRadius) THEN
              !RMask(jnd,ind) = 1
              IPixelTotal= IPixelTotal + 1
              IPixelLocations(IPixelTotal,1) = ind
              IPixelLocations(IPixelTotal,2) = jnd
           ELSE
              !RMask(jnd,ind) = 0
           ENDIF
        ENDDO
     ENDDO
  CASE(1) ! square
     IPixelTotal = 0
     DO ind = 1,2*IPixelCount
        DO jnd = 1,2*IPixelCount
           IPixelTotal = IPixelTotal+1
           IPixelLocations(IPixelTotal,1) = ind
           IPixelLocations(IPixelTotal,2) = jnd
        END DO
     END DO
  END SELECT
  
END SUBROUTINE ImageMaskInitialization
