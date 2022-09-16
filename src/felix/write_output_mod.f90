!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: 2.0
! Date: 31-08-2022
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build: cRED 
! Author:  r.beanland@warwick.ac.uk
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
!! Writes output files for each iteration


MODULE write_output_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: WriteIterationOutput

  CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Writes output for interations including simulated .bin files, 
  !! structureFactors.txt and structure.cif 
  !!
  !! Major-Authors: 'kidwhizz' (2015), Richard Beanland (2016)
  !!
  SUBROUTINE WriteIterationOutput(IThicknessIndex,IErr)

    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE IPARA, ONLY : ILN,ISizeX,ISizeY,IhklsFrame,INoOfHKLsFrame,IFrame,IByteSize
    USE RPARA, ONLY : Rhkl, RImageSimi, RInitialThickness, RDeltaThickness, RBrightField, RTempImage
    USE SPARA, ONLY : SChemicalFormula
    USE IChannels, ONLY : IChOutWIImage

    IMPLICIT NONE

    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND), INTENT(IN) :: IThicknessIndex
    INTEGER(IKIND) :: IThickness,ind,jnd,knd
    REAL(RKIND),DIMENSION(ISizeY,ISizeX) :: RImageToWrite
    CHARACTER(10) :: hString,kString,lString
    CHARACTER(4) :: fString
    CHARACTER(40) :: SPrintString
    CHARACTER(200) :: path,filename,fullpath
    
    IErr=0 

    ! folder names, NB if numbers are too big we get a output conversion error here
    IThickness = NINT(RInitialThickness +(IThicknessIndex-1)*RDeltaThickness)/10.0!in nm 
    IF(ISizeX.GT.99999)THEN
      WRITE(path,"(I3.3,A3,I6,A1,I3.3)") IThickness,"nm_",ISizeX,"x",ISizeY
    ELSEIF(ISizeX.GT.9999)THEN
      WRITE(path,"(I3.3,A3,I5,A1,I3.3)") IThickness,"nm_",ISizeX,"x",ISizeY
    ELSEIF(ISizeX.GT.999)THEN
      WRITE(path,"(I3.3,A3,I4,A1,I3.3)") IThickness,"nm_",ISizeX,"x",ISizeY
    ELSE
      WRITE(path,"(I3.3,A3,I3.3,A1,I3.3)") IThickness,"nm_",ISizeX,"x",ISizeY

    ENDIF
    path = SChemicalFormula(1:ILN) // "_" // path ! This adds chemical to folder name
    CALL system('mkdir ' // path)

    knd = 0
    ! Write Images to disk
    DO ind = 1,INoOfHKLsFrame
      RImageToWrite = RImageSimi(:,:,ind,IThicknessIndex)
      ! only output an image if it has signifcant intensity (>0.001, 0.1% of incident beam)
      IF(MAXVAL(RImageToWrite).GT.0.001D0) THEN
        knd = knd + 1
        ! Make the hkl string e.g. -2-2+10
        jnd=NINT(Rhkl(IhklsFrame(ind),1))
        IF (ABS(jnd).LT.10) THEN
          WRITE(hString,"(SP,I2.1)") jnd
        ELSEIF (ABS(jnd).LT.100) THEN
          WRITE(hString,"(SP,I3.1)") jnd
        ELSE
          WRITE(hString,"(SP,I4.1)") jnd
        ENDIF
        jnd=NINT(Rhkl(IhklsFrame(ind),2))
        IF (ABS(jnd).LT.10) THEN
          WRITE(kString,"(SP,I2.1)") jnd
        ELSEIF (ABS(jnd).LT.100) THEN
          WRITE(kString,"(SP,I3.1)") jnd
        ELSE
          WRITE(kString,"(SP,I4.1)") jnd
        ENDIF
        jnd=NINT(Rhkl(IhklsFrame(ind),3))
        IF (ABS(jnd).LT.10) THEN
          WRITE(lString,"(SP,I2.1)") jnd
        ELSEIF (ABS(jnd).LT.100) THEN
          WRITE(lString,"(SP,I3.1)") jnd
        ELSE
          WRITE(lString,"(SP,I4.1)") jnd
        ENDIF
        ! Make the path/filenames e.g. 'GaAs_-2-2+0.bin'
        WRITE(fString,"(I4.4)") IFrame
        filename = fString // "_" // SChemicalFormula(1:ILN) // "_" // &
          TRIM(ADJUSTL(hString)) // TRIM(ADJUSTL(kString)) // TRIM(ADJUSTL(lString)) // ".bin"
        fullpath = TRIM(ADJUSTL(path))//"/"//TRIM(ADJUSTL(filename))
        CALL message ( LL, dbg6, fullpath )

        ! Writes data to output image .bin files
        OPEN(UNIT=IChOutWIImage, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(fullpath)),&
          FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=ISizeX*IByteSize)
        IF(l_alert(IErr,"WriteIterationOutput","OPEN() output .bin file")) RETURN
        ! we write each X-line, remembering indexing is [row,col]=[y,x]
        DO jnd = 1,ISizeY
          WRITE(IChOutWIImage,rec=jnd) RImageToWrite(jnd,:)
        END DO
        CLOSE(IChOutWIImage,IOSTAT=IErr) 
        IF(l_alert(IErr,"WriteIterationOutput","CLOSE() output .bin file")) RETURN       
      END IF
    END DO

    !--------------------------------------------------------------------
    ! output how many useful reflections have been calculated
    IF(IFrame.GT.9999)THEN
      WRITE(SPrintString,'(A6,I3,A22,I5)') "Found ",knd," reflections in Frame ",IFrame
    ELSEIF(IFrame.GT.999)THEN
      WRITE(SPrintString,'(A6,I3,A22,I4)') "Found ",knd," reflections in Frame ",IFrame
    ELSEIF(IFrame.GT.99)THEN
      WRITE(SPrintString,'(A6,I3,A22,I3)') "Found ",knd," reflections in Frame ",IFrame
    ELSEIF(IFrame.GT.9)THEN
      WRITE(SPrintString,'(A6,I3,A22,I2)') "Found ",knd," reflections in Frame ",IFrame
    ELSE
      WRITE(SPrintString,'(A6,I3,A22,I1)') "Found ",knd," reflections in Frame ",IFrame
    END IF
    CALL message(LS,SPrintString)

    !--------------------------------------------------------------------
    ! write bright field image
    filename = fString // "_000.bin" 
    fullpath = TRIM(ADJUSTL(path))//"/"//TRIM(ADJUSTL(filename))
    OPEN(UNIT=IChOutWIImage, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(fullpath)),&
      FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=SIZE(RBrightField,DIM=2)*IByteSize)
    IF(l_alert(IErr,"WriteIterationOutput","OPEN() output 000.bin file")) RETURN
    ! we write each X-line, remembering indexing is [row,col]=[y,x]
    DO jnd = 1,ISizeY
      WRITE(IChOutWIImage,rec=jnd) RBrightField(jnd,:)
    END DO
    CLOSE(IChOutWIImage,IOSTAT=IErr)
    IF(l_alert(IErr,"WriteIterationOutput","CLOSE() output .bin file")) RETURN



    RETURN  
    
  END SUBROUTINE WriteIterationOutput


END MODULE write_output_mod

