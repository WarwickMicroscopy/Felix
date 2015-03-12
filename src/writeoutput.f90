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
  USE WriteToScreen

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
    
  CALL Message("WriteOutput",IMust,IErr)

  ALLOCATE( &
       RImage(2*IPixelCount,2*IPixelCount), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteOutput(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RImage"
     RETURN
  ENDIF

  CALL Message("WriteOutput",IAllInfo,IErr,MessageString = "Writing Images")
  
  DO knd = 1,IThicknessCount
     
     !--------------------------------------------------------
     ! Write Montage
     !--------------------------------------------------------
     
     RThickness = RInitialThickness + (knd-1)*RDeltaThickness 
     IThickness = RInitialThickness + (knd-1)*RDeltaThickness 
     
     WRITE(CThickness,*) IThickness
     WRITE(CThicknessLength,*) SCAN(CThickness,'0123456789',.TRUE.)-SCAN(CThickness,'0123456789')+1
     
     
     IF(IImageFLAG.EQ.0.OR.IImageFLAG.EQ.2.OR.IImageFLAG.EQ.4.OR.IImageFLAG.EQ.6) THEN

        WRITE(surname,"(A2,I1.1,I1.1,I1.1,I1.1,A2,I5.5,A2,I5.5,A2,I5.5)") &
             "f-",&
             IScatterFactorMethodFLAG, &
             IZolzFLAG, &
             IAbsorbFLAG, &
             IAnisoDebyeWallerFactorFlag,&
             "-T",IThickness,&
             "-P",MAXVAL(IImageSizeXY),&
             "-P",MAXVAL(IImageSizeXY)
!!$        WRITE(surname,"(A2,A1,I5.5,A2,I5.5)") &
!!$             "M-","T",IThickness,"-P",MAXVAL(IImageSizeXY)
        
        CALL OpenReflectionImage(MontageOut,surname,IErr,0,MAXVAL(IImageSizeXY),knd)
        IF( IErr.NE.0 ) THEN
           PRINT*,"WriteOutput(", my_rank, ") error in OpenData()"
           RETURN
        ENDIF
           

        CALL Message("WriteOutput",IMoreInfo,IErr,MessageVariable = "working on RThickness",RVariable = RThickness )

        CALL WriteReflectionImage(MontageOut,RFinalMontageImageRoot(:,:,knd), &
             IErr,MAXVAL(IImageSizeXY),MAXVAL(IImageSizeXY),knd)

        CLOSE(MontageOut,IOSTAT =IErr)
        CALL ErrorChecks("WriteOutput","WriteOutput",ICritError,IErr)
        
        
     END IF
     
     !--------------------------------------------------------
     ! Write Reflections
     !--------------------------------------------------------
     
     IF(IImageFLAG.EQ.1.OR.IImageFLAG.EQ.2.OR.IImageFLAG.EQ.5.OR.IImageFLAG.EQ.6) THEN
        
        WRITE(path,"(A2,I1.1,I1.1,I1.1,I1.1,A2,I5.5,A2,I5.5,A2,I5.5)") &
             "f-",&
             IScatterFactorMethodFLAG, &
             IZolzFLAG, &
             IAbsorbFLAG, &
             IAnisoDebyeWallerFactorFlag,&
             "-T",IThickness,&
             "-P",2*IPixelcount,&
             "-P",2*IPixelcount

        call system('mkdir ' // path)
   
        DO ind = 1,IReflectOut
           CALL OpenReflectionImage(IChOutWIImage,path,IErr,ind,2*IPixelCount,knd)
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error in OpenReflectionImage()"
              RETURN
           ENDIF
          
           RImage = ZERO
           DO jnd = 1,IPixelTotal
              gnd = IPixelLocations(jnd,1)
              hnd = IPixelLocations(jnd,2)
              RImage(gnd,hnd) = RIndividualReflectionsRoot(ind,knd,jnd) 
              
           END DO
           
!!$           RImage = TRANSPOSE(RImage);

           CALL WriteReflectionImage(IChOutWIImage,&
                RImage,IErr,2*IPixelCount,2*IPixelCount,knd)       
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error in WriteReflectionImage()"
              RETURN
           ENDIF    
           CLOSE(IChOutWIImage,IOSTAT = IErr)
           CALL ErrorChecks("WriteOutput","WriteOutput",ICritError,IErr)
        END DO
     END IF
     
     IF(IImageFLAG.GE.3) THEN
        
!!$        WRITE(path,"(A2,A1,I1.1,A2,I1.1,A2,I1.1,A2,I4.4,A2,I5.5)") &
!!$             "F-",&
!!$             "S", IScatterFactorMethodFLAG, &
!!$             "_B", ICentralBeamFLAG, &
!!$             "_M", IMaskFLAG, &
!!$             "_P", IPixelCount, &
!!$             "_T", IThickness
        
        WRITE(path,"(A2,I1.1,I1.1,I1.1,I1.1,A2,I5.5,A2,I5.5,A2,I5.5)") &
             "f-",&
             IScatterFactorMethodFLAG, &
             IZolzFLAG, &
             IAbsorbFLAG, &
             IAnisoDebyeWallerFactorFlag,&
             "-T",IThickness,&
             "-P",2*IPixelcount,&
             "-P",2*IPixelcount

        call system('mkdir ' // path)
        
        DO ind = 1,IReflectOut
           CALL OpenReflectionImage(IChOutWFImageReal,path,IErr,ind,2*IPixelCount,knd)
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error in OpenAmplitudeImage()"
              RETURN
           ENDIF
           
           CALL OpenReflectionImage(IChOutWFImagePhase,path,IErr,ind,2*IPixelCount,knd)
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
           END DO
              
           CALL WriteReflectionImage(IChOutWFImageReal,&
                RImage,IErr,2*IPixelCount,2*IPixelCount,knd)       
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
                RImage,IErr,2*IPixelCount,2*IPixelCount,knd)       
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error in WriteReflectionImage()"
              RETURN
           ENDIF
           
           
           CLOSE(IChOutWFImageReal,IOSTAT = IErr)
           CALL ErrorChecks("WriteOutput","WriteOutput",ICritError,IErr)

           CLOSE(IChOutWFImagePhase,IOSTAT = IErr)
           CALL ErrorChecks("WriteOutput","WriteOutput",ICritError,IErr)
         
        END DO
     END IF
  END DO

!!$  Resets the Message Counter (For future entering subroutine messages)
  IMessageCounter = 0  
  
 !IF (my_rank ==0) THEN
 !END IF

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

END SUBROUTINE WriteOutput
