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
!!  The following describes DM3 file format: https://imagej.nih.gov/ij/plugins/DM3Format.gj.html
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

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ReadDM3TagsAndImage 

  CONTAINS
    
  !>
  !! Procedure-description: This reads in bytes from a .dm3 file, recognises the '%%%%' 
  !! delimiters, reads the TagEntry labels and then the image data. Image data is assumed to be 
  !! 32bit float big endian and under the 2nd 'Data' tag.
  !!
  !! Major-Authors: Jacob Richardson (2017)
  !!
  SUBROUTINE ReadDM3TagsAndImage( SFilePath, IYPixels, IXPixels, IErr, RImageMatrixDM3 )

    ! NB The file is read-in with big endian to match the image data format, however
    ! the tag data is little endian format. The tag characters are read-in byte by byte and
    ! the pixel size is read-in incorrectly and then converted to the little endian.

    !?? despite the pixel sizes RImageMatrixDM3 is outputted as square using max side, refer to module notes above

    CHARACTER(*),INTENT(IN) :: SFilePath
    INTEGER(4), INTENT(IN) :: IYPixels, IXPixels !?? despite these, RImageArray is outputted as square
    INTEGER(4),INTENT(OUT) :: IErr
    REAL(4),INTENT(OUT) :: RImageMatrixDM3(MAX(IYPixels,IXPixels),MAX(IYPixels,IXPixels))
    ! RImageMatrixDM3( x_pixel, y_pixel )

    ! parameters (some of these could be inputs but currently unncessary)
    LOGICAL,PARAMETER  ::  &
      LGetImageData = .TRUE.,&
      LPrintTags = .FALSE.,&
      LPrintingAllowed = .FALSE.       ! if false, excluding errors, printing to terminal is suppressed
    INTEGER(1),PARAMETER :: &
      ICorrespondingImageDataTag = 2   ! specifies which 'Data' tag corresponds to image data 

    ! local variables
    INTEGER(4) :: i,j, INoOfBytes, INoOfTags, INoOfDataTags, INoOfPreDataBytes, Iy, Ix
    LOGICAL :: LFindingTags, LReadingImageData, LLookingForDataTags
    CHARACTER(36) :: STagLabel
    INTEGER(1) :: IPrevious4Bytes(4), IPreviousBytes(40), IByte
    INTEGER(4) :: I4bytePreData, IDataLengthBigEndian
    REAL(4) :: RDataBytes ! image data is assumed to be 32bit float type Little Endian
    ! intialise output variables
    RImageMatrixDM3=0 
    IErr=0
    ! intialise local variables
    LFindingTags=.TRUE.; LReadingImageData=.FALSE.; LLookingForDataTags=.TRUE.
    INoOfBytes=0; INoOfTags=0; INoOfPreDataBytes=0; INoOfDataTags=0; Iy=0; Ix=1;
    IPrevious4Bytes=0; IPreviousBytes=0; I4bytePreData=0;
    RDataBytes=0; IDataLengthBigEndian=0;

    OPEN(UNIT=1,FILE=SFilePath,STATUS='OLD',ACCESS='STREAM',ACTION='READ',IOSTAT=IErr, CONVERT='LITTLE_ENDIAN')
    IF(IErr.NE.0) THEN
      WRITE(*,'(A)') 'Error in ReadDM3TagsAndImage(). Opening .dm3 file, file path =', TRIM(SFilePath)
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
              EXIT ! non ASCII character found so no longer reading label
            ENDIF
          ENDDO
          IF(LPrintTags.AND.LPrintingAllowed) WRITE(*,'(A36,2x,I5,I15)') STagLabel,INoOfTags,INoOfBytes
          IF(LGetImageData.AND.LLookingForDataTags) THEN
            IF(STagLabel.EQ.'Data') INoOfDataTags = INoOfDataTags + 1
            IF(INoOfDataTags.EQ.ICorrespondingImageDataTag) THEN
              LLookingForDataTags=.FALSE.
              IF(LPrintingAllowed) WRITE(*,'(A)') 'Now reading image data'
              LFindingTags = .FALSE.
              LReadingImageData = .TRUE.
            ENDIF
          ENDIF
        ENDIF
  
      ELSEIF ( LReadingImageData ) THEN ! read another set of bytes of image data
        
        ! There are 16 pre-data bytes including an integer of the data array length
        ! which should match the pixel size
        IF(INoOfPreDataBytes.LT.16) THEN         
          READ(1,IOSTAT=IErr) I4bytePreData
          IF(IErr.NE.0) THEN ! error reading dm3
            WRITE(*,'(A,I0)') 'Error in ReadDM3TagsAndImage(). Reading pre data bytes, Byte number = ',INoOfBytes
            RETURN
          ENDIF
          IF(INoOfPreDataBytes.EQ.12) THEN
            ! Convert data length integer to big endian format
            CALL MVBITS( I4bytePreData, 24, 8, IDataLengthBigEndian, 0  )
            CALL MVBITS( I4bytePreData, 16, 8, IDataLengthBigEndian, 8  )
            CALL MVBITS( I4bytePreData, 8,  8, IDataLengthBigEndian, 16 )
            CALL MVBITS( I4bytePreData, 0,  8, IDataLengthBigEndian, 24 )
            IF(IDataLengthBigEndian.NE.IYPixels*IXPixels) THEN ! error
              IErr=1
              WRITE(*,'(A,I0)') 'Error in ReadDM3TagsAndImage(). Data array length does not match inputted pixel size.'
              WRITE(*,'(A,I0,A,I0)') 'Error in ReadDM3TagsAndImage(). Data array length = ',IDataLengthBigEndian,', IYPixels*IXPixels = ', IYPixels*IXPixels
              RETURN
            ENDIF 
          ENDIF
          INoOfBytes = INoOfBytes + 4
          INoOfPreDataBytes = INoOfPreDataBytes + 4
        ELSE ! read image data bytes
          READ(1,IOSTAT=IErr) RDataBytes
          IF(IErr.NE.0) THEN ! error reading dm3
            WRITE(*,'(A,I0)') 'Error in ReadDM3TagsAndImage(). while reading data, Byte number = ',INoOfBytes
            RETURN
          ENDIF
          INoOfBytes = INoOfBytes + 4
          Iy = Iy + 1
          IF(Iy .EQ. IYPixels + 1) THEN
            Iy = 1
            Ix = Ix + 1
          ENDIF
          IF(Ix .LT. IXPixels + 1) THEN
            RImageMatrixDM3(Ix,Iy) = RDataBytes
          ELSE ! reached end of image data
            LFindingTags = .TRUE.
            LReadingImageData = .FALSE.
          ENDIF                   
        ENDIF

      ELSE ! error, this should never happen
        IErr=1; WRITE(*,'(A)') 'Error in ReadDM3TagsAndImage(). Reading DM3, cannot search for tags nor read image data'; RETURN

      ENDIF

    ENDDO    
    CLOSE(1)
  
  END SUBROUTINE ReadDM3TagsAndImage


    
END MODULE
