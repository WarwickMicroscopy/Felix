!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all right reserved
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: specimentsetup.f90,v 1.59 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE WriteOutput( CAmplitudeandPhaseRoot,RIndividualReflectionsRoot,RFinalMontageImageRoot,IErr)

  USE MyNumbers
  USE CPara; USE IPara; USE SPara
  USE RPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,knd,gnd,hnd
  INTEGER(IKIND) IThickness, IErr

  REAL(RKIND) RThickness
  REAL(RKIND),DIMENSION(IReflectOut,IThicknessCount,IPixelTotal):: &
       RIndividualReflectionsRoot
  REAL(RKIND),DIMENSION(MAXVAL(IImageSizeXY),&
          MAXVAL(IImageSizeXY),IThicknessCount):: RFinalMontageImageRoot
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RImage

  COMPLEX(CKIND),DIMENSION(IReflectOut,IThicknessCount,IPixelTotal):: &
       CAmplitudeandPhaseRoot

  CHARACTER*40 surname, path
  CHARACTER*25 CThickness, CThicknessLength
  
  !IF (my_rank.EQ.0) THEN

  ALLOCATE( &
       RImage(2*IPixelCount,2*IPixelCount), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteOutput(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RImage"
     RETURN
  ENDIF
  
  
  IF(IWriteFLAG.GE.10) THEN
     PRINT*,"Writing Images" 
  END IF
  
  
  DO knd = 1,IThicknessCount
     
     !--------------------------------------------------------
     ! Write Montage
     !--------------------------------------------------------
     
     RThickness = RInitialThickness + (knd-1)*RDeltaThickness 
     IThickness = RInitialThickness + (knd-1)*RDeltaThickness 
     
     WRITE(CThickness,*) IThickness
     WRITE(CThicknessLength,*) SCAN(CThickness,'0123456789',.TRUE.)-SCAN(CThickness,'0123456789')+1
     
     
     IF(IImageFLAG.EQ.0.OR.IImageFLAG.EQ.2.OR.IImageFLAG.EQ.4.OR.IImageFLAG.EQ.6) THEN
        WRITE(surname,"(A2,A1,I5.5,A2,I5.5)") &
             "M-","T",IThickness,"-P",MAXVAL(IImageSizeXY)
        
        CALL OpenReflectionImage(MontageOut,surname,IErr,0,MAXVAL(IImageSizeXY)) 
        IF( IErr.NE.0 ) THEN
           PRINT*,"WriteOutput(", my_rank, ") error in OpenData()"
           RETURN
        ENDIF
        
        IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN         
           PRINT*,"WriteOutput(", my_rank, ") working on RThickness=", RThickness
        END IF
        
        CALL WriteReflectionImage(MontageOut,RFinalMontageImageRoot(:,:,knd), &
             IErr,MAXVAL(IImageSizeXY),MAXVAL(IImageSizeXY))
        CLOSE(MontageOut,ERR=9999)
        
     END IF
     
     !--------------------------------------------------------
     ! Write Reflections
     !--------------------------------------------------------
     
     
     IF(IImageFLAG.EQ.1.OR.IImageFLAG.EQ.2.OR.IImageFLAG.EQ.5.OR.IImageFLAG.EQ.6) THEN
        
        WRITE(path,"(A2,A1,I1.1,A2,I1.1,A2,I1.1,A2,I4.4,A2,I5.5)") &
             "F-",&
             "S", IScatterFactorMethodFLAG, &
             "_B", ICentralBeamFLAG, &
             "_M", IMaskFLAG, &
             "_P", IPixelCount, &
             "_T", IThickness
        
        call system('mkdir ' // path)
        
        DO ind = 1,IReflectOut
           CALL OpenReflectionImage(IChOutWIImage,path,IErr,ind,2*IPixelCount)
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error in OpenReflectionImage()"
              RETURN
           ENDIF
           
           RImage = ZERO
           DO jnd = 1,IPixelTotal
              gnd = IPixelLocations(jnd,1)
              hnd = IPixelLocations(jnd,2)
              RImage(gnd,hnd) = RIndividualReflectionsRoot(ind,knd,jnd)
              PRINT*, RImage
           END DO
           
           CALL WriteReflectionImage(IChOutWIImage,&
                RImage,IErr,2*IPixelCount,2*IPixelCount)       
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error in WriteReflectionImage()"
              RETURN
           ENDIF
           CLOSE(IChOutWIImage,ERR=9999)
        END DO
     END IF
     
     IF(IImageFLAG.GE.3) THEN
        
        WRITE(path,"(A2,A1,I1.1,A2,I1.1,A2,I1.1,A2,I4.4,A2,I5.5)") &
             "F-",&
             "S", IScatterFactorMethodFLAG, &
             "_B", ICentralBeamFLAG, &
             "_M", IMaskFLAG, &
             "_P", IPixelCount, &
             "_T", IThickness
        
        call system('mkdir ' // path)
        
        DO ind = 1,IReflectOut
           CALL OpenReflectionImage(IChOutWFImageReal,path,IErr,ind,2*IPixelCount)
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error in OpenAmplitudeImage()"
              RETURN
           ENDIF
           
           CALL OpenReflectionImage(IChOutWFImagePhase,path,IErr,ind,2*IPixelCount)
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error in OpenPhaseImage()"
              RETURN
           ENDIF
           
           !-----------------------------------------------------------------------------
           ! Create An Image
           !-----------------------------------------------------------------------------
           
           RImage = ZERO
           DO jnd = 1,IPixelTotal
              gnd = IPixelLocations(jnd,1)
              hnd = IPixelLocations(jnd,2)
              RImage(gnd,hnd) = REAL(CAmplitudeandPhaseRoot(ind,knd,jnd))
              PRINT*, RImage
           END DO
              
           CALL WriteReflectionImage(IChOutWFImageReal,&
                RImage,IErr,2*IPixelCount,2*IPixelCount)       
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error in WriteReflectionImage()"
              RETURN
           ENDIF
           
           RImage = ZERO
           DO jnd = 1,IPixelTotal
              gnd = IPixelLocations(jnd,1)
              hnd = IPixelLocations(jnd,2)
              RImage(gnd,hnd) = AIMAG(CAmplitudeandPhaseRoot(ind,knd,jnd))
           END DO
           
           CALL WriteReflectionImage(IChOutWFImagePhase,&
                RImage,IErr,2*IPixelCount,2*IPixelCount)       
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error in WriteReflectionImage()"
              RETURN
           ENDIF
           
           
           CLOSE(IChOutWFImageReal,ERR=9999)
           CLOSE(IChOutWFImagePhase,ERR = 9999)
        END DO
     END IF
  END DO
  
  DEALLOCATE( &
       RImage,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteOutput(", my_rank, ") error in Deallocation of RImage"
     RETURN
  ENDIF
  
  ! END IF
  
     
  !--------------------------------------------------------------------
  ! Shut down MPI
  !--------------------------------------------------------------------
9999 &
     CALL MPI_Finalize(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteOutput(", my_rank, ") error ", IErr, " in MPI_Finalize()"
        STOP
     ENDIF
     
     ! clean shutdown
     STOP
  

END SUBROUTINE WriteOutput
