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
!!  Module-Description: This is a generic fortran module for reading in data from dm3 files.
!!  It is generally not felix-specific but may specify types and formatting to match that are
!!  expected for felix code.
!!
!!  Major-Authors: Jacob Richardson (2017)
!!
MODULE read_dm3_mod

  ! SUBROUTINE (INPUTS..., IErr, OUTPUTS... ) is the convention followed here
  ! For test programs using this module, e.g. use felix make to compile read_dm3_mod, 
  ! then in src/felix/, gfortran -o use_dm3_module use_dm3_module.f90 dread_dm3_mod.o

  !?? aside from byte reading should standardise real and int kind

  !?? despite the pixel sizes RImageMatrixDM3,RImageArray are outputted as square using max side
  !?? felix convert scripts currently requires square output images, but to read-in .dm3 need rectangular

  !?? felix KIND(1.0D0) is 64bit, corresponds to 8bytes for RECORD LENGTH

  !?? standardise convention here of x, y axis of image from top left corner...

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ReadDM3TagsAndImage, ReadInDM3directory, WriteOutImageArrayDirectory 

  CONTAINS
    
  !>
  !! Procedure-description: This reads in bytes from a .dm3 file, recognises the '%%%%' 
  !! delimiters, reads the tag labels and then the image data. Image data is assumed to be 
  !! 32bit float type and in the file in an array under the 2nd 'Data' tag.
  !!
  !! Major-Authors: Jacob Richardson (2017)
  !!
  SUBROUTINE ReadDM3TagsAndImage( SFilePath, IXPixels, IYPixels, IErr, RImageMatrixDM3 )

    !?? despite the pixel sizes RImageMatrixDM3 is outputted as square using max side, refer to module notes above

    CHARACTER(*),INTENT(IN) :: SFilePath
    INTEGER(4), INTENT(IN) :: IXPixels, IYPixels !?? despite these RImageArray is outputted as square
    INTEGER(4),INTENT(OUT) :: IErr
    REAL(4),INTENT(OUT) :: RImageMatrixDM3(MAX(IXPixels,IYPixels),MAX(IXPixels,IYPixels))
    ! RImageMatrixDM3( x_pixel, y_pixel )

    ! parameters (some of these could be inputs but currently unncessary)
    LOGICAL,PARAMETER  ::  &
      LGetImageData = .TRUE.,&
      LPrintTags = .FALSE.,&
      LPrintingAllowed = .FALSE.       ! if false, excluding errors printing to terminal is suppressed
    INTEGER(1),PARAMETER :: &
      ICorrespondingImageDataTag = 2  ! specifies which 'Data' tag corresponds to image data 

    ! local variables
    INTEGER(4) :: i,j, INoOfBytes, INoOfTags, INoOfDataTags, INoOfIgnoredBytes, Ix, Iy
    LOGICAL :: LFindingTags, LReadingImageData, LLookingForDataTags
    CHARACTER(36) :: STagLabel
    INTEGER(1) :: IPrevious4Bytes(4), IPreviousBytes(40)
    INTEGER(4) :: IIgnore4Bytes
    REAL(4) :: RDataBytes ! image data is assumed to be 32bit float type
    ! intialise output variables
    RImageMatrixDM3=0 
    IErr=0
    ! intialise local variables
    LFindingTags=.TRUE.; LReadingImageData=.FALSE.; LLookingForDataTags=.TRUE.
    INoOfBytes=0; INoOfTags=0; INoOfIgnoredBytes=0; INoOfDataTags=0; Ix=0; Iy=1;
    IPrevious4Bytes=0; IPreviousBytes=0; IIgnore4Bytes=0;
    RDataBytes=0;

    OPEN(UNIT=1,FILE=SFilePath,STATUS='OLD',ACCESS='STREAM',ACTION='READ',IOSTAT=IErr)
    IF(IErr.NE.0) THEN
      WRITE(*,'(A)') 'Error opening .dm3 file, file path =', SFilePath
      RETURN
    ENDIF

    DO i = 1,4000000 ! should exit due to end of file before this is reached
      ! each iterationa a set of bytes are read in from the file
      ! depending upon what is being read-in, this can be 1-4 bytes at a time

      IF( LFindingTags ) THEN ! read another byte to search and find tags

        ! remove oldest (leftmost) byte and move list of bytes left by 1
        IPreviousBytes = CSHIFT(IPreviousBytes,1)
        IPrevious4Bytes = CSHIFT(IPrevious4Bytes,1)
        
        ! read byte from file
        READ(1,IOSTAT=IErr) IPreviousBytes(SIZE(IPreviousBytes)) ! read byte into end of IPreviousBytes array
        IF(IErr.LT.0) THEN ! End of file reached
          IF(LPrintTags.AND.LPrintingAllowed) WRITE(*,'(A43,I15)') 'End Of File Reached, Total bytes =         ',INoOfBytes
          IErr=0
          EXIT
        ENDIF
        INoOfBytes = INoOfBytes + 1
        IPrevious4Bytes(4) = IPreviousBytes(SIZE(IPreviousBytes))

        ! check if at a DM3 tag delimiter '%%%%' ('%' is 37 in ASCII)
        ! and read tag label
        IF( ALL(IPrevious4Bytes.EQ.[37,37,37,37]) ) THEN 
          INoOfTags = INoOfTags + 1 
          ! convert bytes before delimiter into corresponding tag label
          STagLabel = ''
          DO j = SIZE(IPreviousBytes)-4,1,-1
            ! only convert ASCII printable bytes into characters for labels
            IF( 32.LE.IPreviousBytes(j) .AND. IPreviousBytes(j).LE.126 ) THEN
              STagLabel = CHAR(IPreviousBytes(j))//TRIM(STagLabel)
            ELSE
              EXIT ! non ASCII character found found so no longer reading label
            ENDIF
          ENDDO
          IF(LPrintTags.AND.LPrintingAllowed) WRITE(*,'(A36,2x,I5,I15)') STagLabel,INoOfTags,INoOfBytes
          IF(LLookingForDataTags) THEN
            IF(LGetImageData.AND.STagLabel.EQ.'Data') INoOfDataTags = INoOfDataTags + 1
            IF(INoOfDataTags.EQ.ICorrespondingImageDataTag) THEN
              LLookingForDataTags=.FALSE.
              IF(LPrintingAllowed) WRITE(*,'(A)') 'Now reading image data'
              LFindingTags = .FALSE.
              LReadingImageData = .TRUE.
            ENDIF
          ENDIF
        ENDIF
  
      ELSEIF ( LReadingImageData ) THEN ! read another set of bytes of image data
        
        ! at the beggining of data array there are 12 bytes to ignore (three 4-byte integers)
        IF(INoOfIgnoredBytes.LT.12) THEN
          READ(1,IOSTAT=IErr) IIgnore4Bytes
          IF(IErr.NE.0) THEN ! error reading dm3
            WRITE(*,'(A,I0)') 'Error while reading intial data bytes to-be-ignored, Byte number = ',INoOfBytes
            RETURN
          ENDIF
          INoOfBytes = INoOfBytes + 4
          INoOfIgnoredBytes = INoOfIgnoredBytes + 4
        ELSE ! read image data bytes
          READ(1,IOSTAT=IErr) RDataBytes
          IF(IErr.NE.0) THEN ! error reading dm3
            WRITE(*,'(A,I0)') 'Error while reading data, Byte number = ',INoOfBytes
            RETURN
          ENDIF
          INoOfBytes = INoOfBytes + 4
          Ix = Ix + 1
          IF(Ix .EQ. IXPixels + 1) THEN
            Ix = 1
            Iy = Iy + 1
          ENDIF
          IF(Iy .LT. IYPixels + 1) THEN
            RImageMatrixDM3(Ix,Iy) = RDataBytes
          ELSE ! reached end of image data
            LFindingTags = .TRUE.
            LReadingImageData = .FALSE.
          ENDIF                   
        ENDIF

      ELSE ! error, this should never happen
        IErr=1; WRITE(*,'(A)') 'Error reading DM3, cannot search for tags nor read image data'; RETURN

      ENDIF

    ENDDO    
    CLOSE(1)
  
  END SUBROUTINE ReadDM3TagsAndImage
  
  
  !>
  !! Procedure-description: Read all .dm3 files in a directory into an array of 2D image matrices.
  !! Manually specify expected filenames in this subroutine. 
  !! Importantly ImageArray becomes square using max IXPixels, IYPixels
  !!
  !! Major-Authors: Jacob Richardson (2017)
  !!
  SUBROUTINE ReadInDM3directory( SDirectoryPath, IXPixels, IYPixels, IErr, RImageArray )
    
    !?? despite the pixel sizes RImageArray is outputted as square using max side, refer to module notes above

    CHARACTER(*),INTENT(IN) :: SDirectoryPath ! should include '/'
    INTEGER(4),INTENT(IN) :: IXPixels,IYPixels
    INTEGER(4),INTENT(OUT) :: IErr
    REAL(KIND(1.0D0)),ALLOCATABLE,INTENT(OUT) :: RImageArray(:,:,:) 
    ! RImageArray( x_pixel, y_pixel, image_number )

    ! manually specify these in code
    CHARACTER(200) :: SFilePath
    INTEGER(4) :: INoOfImages
    
    INTEGER(4) :: i, ImageNumber
    REAL(4),ALLOCATABLE :: RImageMatrixDM3(:,:) ! RImageMatrixDM3( x_pixel, y_pixel )
    ! this real kind matches felix format
    
    ! manually specify how many images expected
    INoOfImages=168

    ALLOCATE(RImageMatrixDM3(MAX(IXPixels,IYPixels),MAX(IXPixels,IYPixels)))
    ALLOCATE(RImageArray(MAX(IXPixels,IYPixels),MAX(IXPixels,IYPixels),INoOfImages))

    DO i = 1,INoOfImages

      ! manually specify expected .dm3 filenames inside directory
      ImageNumber=i
      WRITE(SFilePath,'(A,A,I3.3,A)') TRIM(SDirectoryPath),'006_D-LACBED_',ImageNumber,'.dm3'
      !WRITE(*,'(A,A)') 'Reading .dm3 with file path = ',TRIM(SFilePath)

      ! read in single image into image array   
      CALL ReadDM3TagsAndImage(SFilePath, IXPixels, IYPixels, IErr, RImageMatrixDM3)
      IF(IErr.NE.0) THEN
        WRITE(*,'(A)') 'Error found in ReadDM3TagsAndImage'
        RETURN
      ENDIF
      RImageArray(:,:,i) = REAL(RImageMatrixDM3,KIND(1.0D0))   
      !WRITE(*,'(A,ES8.1,ES8.1)') 'Mix/max value respectively ',MINVAL(RImageMatrixDM3),MAXVAL(RImageMatrixDM3)

    ENDDO

  END SUBROUTINE ReadInDM3directory


  !>
  !! Procedure-description: Write out any real image array into binary .bin files in a directory.
  !! Manually specify output image filenames in this subroutine.
  !!
  !! Major-Authors: Jacob Richardson (2017)
  !!
  SUBROUTINE WriteOutImageArrayDirectory( SDirectoryPath, RImageArray )
    CHARACTER(*),INTENT(IN) :: SDirectoryPath ! should include '/'
    REAL(KIND(1.0D0)),INTENT(IN),DIMENSION(:,:,:) :: RImageArray ! RImageArray( x_pixel, y_pixel, image_number )
    REAL(KIND(1.0D0)),ALLOCATABLE,DIMENSION(:,:,:) :: RImageArrayDummy
    ! this real kind matches felix format
    INTEGER(4) :: i,j,ImageNumber,IXPixels,IYPixels,IErr
    ! manually specify these in code
    CHARACTER(200) :: SFilePath,SSystemCommand

    WRITE(SSystemCommand,'(A,A)') 'mkdir ',TRIM(SDirectoryPath)
    CALL SYSTEM(SSystemCommand)

    DO i = 1,SIZE(RImageArray,3)
      WRITE(*,'(A,ES8.1,ES8.1)') 'Mix/max value respectively ',MINVAL(RImageArray(:,:,i)),MAXVAL(RImageArray(:,:,i))
    ENDDO
  
    ! process image for output
    RImageArrayDummy=RImageArray
    WHERE ( RImageArrayDummy > 15000 )
      RImageArrayDummy = 15000
    END WHERE
    RImageArrayDummy = RImageArrayDummy/MAXVAL(RImageArrayDummy)

    ! iteratively write .bin file
    DO i = 1,SIZE(RImageArray,3)
      ! manually specify .bin output filename
      ImageNumber=i
      IXPixels=SIZE(RImageArray,1)
      IYPixels=SIZE(RImageArray,2)
      WRITE(SFilePath,'(A,A,I3.3,A,I3.3,A,I3.3,A)')&
            TRIM(SDirectoryPath),'006_D-LACBED_',IXPixels,'x',IYPixels,'_',ImageNumber,'.bin'
      ! NB need pixel size in filenames in format 'bla_bla_NxN_bla.bin' for felix convert image script
      WRITE(*,'(A,A)') 'Writing .bin with file path = ',TRIM(SFilePath)

      ! write .bin file
      OPEN(UNIT=2,STATUS='UNKNOWN',FILE=SFilePath,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=672*8,IOSTAT=IErr)
      IF(IErr.NE.0) THEN
        WRITE(*,'(A)') 'Error in OPEN, making .bin file to write to'
        RETURN
      ENDIF
      DO j = 1,SIZE(RImageArray,2)
        WRITE(2,REC=j) RImageArrayDummy(:,j,i)
      END DO
      CLOSE(2) 
    ENDDO

  END SUBROUTINE WriteOutImageArrayDirectory
    
END MODULE
