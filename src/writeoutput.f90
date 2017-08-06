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
! WriteOutput()
! OpenReflectionImage()
! WriteReflectionImage()

!>
!! Procedure-description: Makes output directory, writes montage, writes 
!! reflections, and creates the image
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE WriteOutput( CAmplitudeandPhaseImages,RReflectionImages,RMontageImages,IErr)

  USE MyNumbers

  USE CPara; USE IPara; USE SPara
  USE RPara
  USE BlochPara
    
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,jnd,knd,gnd,hnd,IThickness,IErr
  REAL(RKIND) :: RThickness
  REAL(RKIND),DIMENSION(INoOfLacbedPatterns,IThicknessCount,IPixelTotal):: RReflectionImages
  REAL(RKIND),DIMENSION(MAXVAL(IImageSizeXY),MAXVAL(IImageSizeXY),IThicknessCount):: RMontageImages
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: RImage
  COMPLEX(CKIND),DIMENSION(INoOfLacbedPatterns,IThicknessCount,IPixelTotal):: CAmplitudeandPhaseImages
  CHARACTER*40 :: surname, path
  CHARACTER*25 :: SThickness, SThicknessLength

  ALLOCATE(RImage(2*IPixelCount,2*IPixelCount),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Error:WriteOutput(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RImage"
     RETURN
  ENDIF

  !--------------------------------------------------------
  ! Make an output directory
  call system('mkdir felixsim_output/')!RB
  
  DO knd = 1,IThicknessCount
     
     !--------------------------------------------------------
     ! Write Montage
     RThickness = RInitialThickness + (knd-1)*RDeltaThickness 
     IThickness = RInitialThickness + (knd-1)*RDeltaThickness 
     
     WRITE(SThickness,*) IThickness
     WRITE(SThicknessLength,*) SCAN(SThickness,'0123456789',.TRUE.)-SCAN(SThickness,'0123456789')+1
     
     
     IF(IImageFLAG.EQ.0.OR.IImageFLAG.EQ.2.OR.IImageFLAG.EQ.4.OR.IImageFLAG.EQ.6) THEN

        WRITE(surname,"(A8,I4.4,A4,I5.5,A1,I5.5)") &
            "Montage_",IThickness/10,"_nm_",MAXVAL(IImageSizeXY),"x",MAXVAL(IImageSizeXY)
        
        CALL OpenReflectionImage(MontageOut,surname,IErr,0,MAXVAL(IImageSizeXY),knd)
        IF( IErr.NE.0 ) THEN
           PRINT*,"Error:WriteOutput(", my_rank, ") error in OpenData()"
           RETURN
        ENDIF

        CALL WriteReflectionImage(MontageOut,RMontageImages(:,:,knd), &
             IErr,MAXVAL(IImageSizeXY),MAXVAL(IImageSizeXY),knd)

        CLOSE(MontageOut,IOSTAT =IErr)
        CALL ErrorChecks("WriteOutput","WriteOutput",ICritError,IErr)
        
        
     END IF
     
     !--------------------------------------------------------
     ! Write Reflections
     !--------------------------------------------------------
     
     IF(IImageFLAG.EQ.1.OR.IImageFLAG.EQ.2.OR.IImageFLAG.EQ.5.OR.IImageFLAG.EQ.6) THEN

        WRITE(path,"(A9,I3.3,A3,I4.4,A1,I4.4,A1)") &
            "felixsim_",IThickness/10,"nm_",2*IPixelcount,"x",2*IPixelcount,"/"
        call system('mkdir -p ' // path)
   
        DO ind = 1,INoOfLacbedPatterns
           CALL OpenReflectionImage(IChOutWIImage,path,IErr,ind,2*IPixelCount,knd)
           IF( IErr.NE.0 ) THEN
              PRINT*,"Error:WriteOutput(", my_rank, ") error in OpenReflectionImage()"
              RETURN
           ENDIF
          
           RImage = ZERO
           DO jnd = 1,IPixelTotal
              gnd = IPixelLocations(jnd,1)
              hnd = IPixelLocations(jnd,2)
              RImage(gnd,hnd) = RReflectionImages(ind,knd,jnd) 
              
           END DO
           
           CALL WriteReflectionImage(IChOutWIImage,&
                RImage,IErr,2*IPixelCount,2*IPixelCount,knd)       
           IF( IErr.NE.0 ) THEN
              PRINT*,"Error:WriteOutput(", my_rank, ") error in WriteReflectionImage()"
              RETURN
           ENDIF    
           CLOSE(IChOutWIImage,IOSTAT = IErr)
           CALL ErrorChecks("WriteOutput","WriteOutput",ICritError,IErr)
        END DO
     END IF
     
     IF(IImageFLAG.GE.3) THEN
        
        WRITE(path,"(A2,I1.1,I1.1,I1.1,I1.1,A2,I5.5,A2,I5.5,A2,I5.5)") &
             "f-",&
             IScatterFactorMethodFLAG, &
             IHolzFLAG, &
             IAbsorbFLAG, &
             IAnisoDebyeWallerFactorFlag,&
             "-T",IThickness,&
             "-P",2*IPixelcount,&
             "-P",2*IPixelcount

        call system('mkdir ' // path)
        
        DO ind = 1,INoOfLacbedPatterns
           CALL OpenReflectionImage(IChOutWFImageReal,path,IErr,ind,2*IPixelCount,knd)
           IF( IErr.NE.0 ) THEN
              PRINT*,"WriteOutput(", my_rank, ") error in OpenAmplitudeImage()"
              RETURN
           ENDIF
           
           CALL OpenReflectionImage(IChOutWFImagePhase,path,IErr,ind,2*IPixelCount,knd)
           IF( IErr.NE.0 ) THEN
              PRINT*,"Error:WriteOutput(", my_rank, ") error in OpenPhaseImage()"
              RETURN
           ENDIF
           
           !-----------------------------------------------------------------------------
           ! Create An Image
           !-----------------------------------------------------------------------------
           
           RImage = ZERO
           DO jnd = 1,IPixelTotal
              gnd = IPixelLocations(jnd,1)
              hnd = IPixelLocations(jnd,2)
              RImage(gnd,hnd) = REAL(CAmplitudeandPhaseImages(ind,knd,jnd))
           END DO
              
           CALL WriteReflectionImage(IChOutWFImageReal,&
                RImage,IErr,2*IPixelCount,2*IPixelCount,knd)       
           IF( IErr.NE.0 ) THEN
              PRINT*,"Error:WriteOutput(", my_rank, ") error in WriteReflectionImage()"
              RETURN
           ENDIF
           
           RImage = ZERO
           DO jnd = 1,IPixelTotal
              gnd = IPixelLocations(jnd,1)
              hnd = IPixelLocations(jnd,2)
              RImage(gnd,hnd) = AIMAG(CAmplitudeandPhaseImages(ind,knd,jnd))
           END DO
           
           CALL WriteReflectionImage(IChOutWFImagePhase,&
                RImage,IErr,2*IPixelCount,2*IPixelCount,knd)       
           IF( IErr.NE.0 ) THEN
              PRINT*,"Error:WriteOutput(", my_rank, ") error in WriteReflectionImage()"
              RETURN
           ENDIF
           
           CLOSE(IChOutWFImageReal,IOSTAT = IErr)
           CALL ErrorChecks("WriteOutput","WriteOutput",ICritError,IErr)
           CLOSE(IChOutWFImagePhase,IOSTAT = IErr)
           CALL ErrorChecks("WriteOutput","WriteOutput",ICritError,IErr)
         
        END DO
     END IF
  END DO

  DEALLOCATE(RImage,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"Error:WriteOutput(", my_rank, ") error in Deallocation of RImage"
     RETURN
  ENDIF
   
END SUBROUTINE WriteOutput

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description: Open Reflection Image
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE OpenReflectionImage(IChOutWrite, surname, IErr,IReflectWriting,IImageSizeX,ind)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara
  USE terminal_output
  USE MPI
  USE MyMPI

  USE IChannels
  IMPLICIT NONE

  CHARACTER(*) :: surname
  CHARACTER*20 :: prefix,postfix,h,k,l
  INTEGER(IKIND) :: IChOutWrite, IErr,IReflectWriting,IImageSizeX
  CHARACTER*250 filename
  CHARACTER*40 fileext
  CHARACTER*60 Simagesize
  INTEGER index,ind

  SELECT CASE(IChOutWrite)
  CASE(MontageOut)
  CASE DEFAULT
     IF(IHKLSelectFLAG.EQ.0) THEn
        WRITE(h,*)  NINT(Rhkl(IReflectWriting,1))
        WRITE(k,*)  NINT(Rhkl(IReflectWriting,2))
        WRITE(l,*)  NINT(Rhkl(IReflectWriting,3))
     ELSE
        WRITE(h,*)  NINT(Rhkl(IOutPutReflections(IReflectWriting),1))
        WRITE(k,*)  NINT(Rhkl(IOutPutReflections(IReflectWriting),2))
        WRITE(l,*)  NINT(Rhkl(IOutPutReflections(IReflectWriting),3))
     END IF
  END SELECT

  WRITE(Simagesize,"(A2,I3.3,A2,I3.3)") &
       "_",IImageSizeX,&
       "x",IImageSizeX

  WRITE(fileext,*) TRIM(ADJUSTL(".bin")) 
  
  SELECT CASE(IChOutWrite)
  
    CASE(IChOutWFImageReal)        
      CALL message ( LL, "OpenImage: opening image for WAVE FUNCTION REAL PART (WR*.txt)" )
      WRITE(filename,*) TRIM(ADJUSTL(surname)),"/Real-",&
      TRIM(ADJUSTL(h)),&
      TRIM(ADJUSTL(k)),&
      TRIM(ADJUSTL(l)),&
      TRIM(ADJUSTL(fileext))
	   
    CASE(IChOutWFImagePhase)        
      CALL message ( LL, "OpenImage: opening image for WAVE FUNCTION PHASE PART (WP*.txt)")
      WRITE(filename,*) TRIM(ADJUSTL(surname)),"/Imag-",&
      TRIM(ADJUSTL(h)),&
      TRIM(ADJUSTL(k)),&
      TRIM(ADJUSTL(l)),&
      TRIM(ADJUSTL(fileext))
	   
    CASE(IChOutWIImage) 
      CALL message ( LL, "OpenImage: opening image for INTENSITIES" )
      WRITE(filename,*) TRIM(ADJUSTL(surname)),"/",&
      TRIM(ADJUSTL(h)),&
      TRIM(ADJUSTL(k)),&
      TRIM(ADJUSTL(l)),&
      TRIM(ADJUSTL(fileext))
	   
    CASE(MontageOut)        
      WRITE(filename,*) TRIM(ADJUSTL(surname)),"Montage_",&
            TRIM(ADJUSTL(fileext))
      CALL message ( LL, "OpenImage: opening image for WAVE INTENSITIES" )

    CASE DEFAULT
      CALL message ( LL, "OpenImage: opening UNKNOWN channel ", IChOutWrite )
	 
  END SELECT

  OPEN(UNIT=IChOutWrite, ERR=10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),FORM='UNFORMATTED',&
       ACCESS='DIRECT',IOSTAT=Ierr,RECL=IImageSizeX*8)
  RETURN
   
  ! error in OPEN detected
10 PRINT*,"Error:OpenReflectionImage(): ERR in OPEN()",IErr
  IErr= 1
  RETURN
  
END SUBROUTINE OpenReflectionImage

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!>
!! Procedure-description: Write reflection images
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE WriteReflectionImage( IChOutWrite, data, IErr,IImageSizeX,IImageSizeY,knd)
  !this is now redundant
  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(KIND=IKIND) IErr,IImageSizeY,IImageSizeX
  REAL(KIND=RKIND),DIMENSION(IImageSizeX,IImageSizeY) :: data
  CHARACTER*100 CSizeofData
  INTEGER ind,knd, IChOutWrite
  CHARACTER*100 SFormatString
     
  DO ind = 1,(IImageSizeY)
     WRITE(IChOutWrite,rec=ind) data(ind,:)
  END DO

  RETURN
  ! error in WRITE detected
20 PRINT*,"Error:WriteReflectionImage(): ERR in WRITE()",Ierr
  IErr= 1
  RETURN

  ! buffer too short
30 PRINT*,"Error:WriteImageR_MPI(): EOF in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteReflectionImage
