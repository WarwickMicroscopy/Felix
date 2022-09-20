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
    USE IPARA, ONLY : ILN,ISizeX,ISizeY,IhklsFrame,INoOfHKLsFrame,IFrame,&
            IhklsAll,ILiveList,IByteSize
    USE RPARA, ONLY : RInputHKLs,Rhkl, RImageSimi, RInitialThickness, RDeltaThickness,&
            RBrightField, RTempImage,RDevPara,RDarkField_1,RDarkField_2,RDarkField_3,&
            RDarkField_4,RDarkField_5,RDarkField_6,RDarkField_7,RDarkField_8,&
            RDarkField_9,RDarkField_10,RDarkField_11,RDarkField_12,RDarkField_13,&
            RDarkField_14,RDarkField_15,RDarkField_16,RDarkField_17,RDarkField_18,&
            RDarkField_19,RDarkField_20
    USE SPARA, ONLY : SChemicalFormula
    USE IChannels, ONLY : IChOutWIImage

    IMPLICIT NONE

    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND), INTENT(IN) :: IThicknessIndex
    INTEGER(IKIND) :: IThickness,ind,jnd,knd,Iflag
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

    path = SChemicalFormula(1:ILN) // "_" // path ! This adds chemical formula to folder name
    IF (IFrame.EQ.1) CALL system('mkdir ' // path)

    knd = 0
    ! Compose images and write to disk
    DO ind = 1,INoOfHKLsFrame!the number of reflections in both felix.hkl and the beam pool
      ! this reflection is number in in RImageSimi
      RImageToWrite = RImageSimi(:,:,ind,IThicknessIndex)
      ! its index in the beam pool is IhklsFrame(ind), its g-vector is Rhkl(IhklsFrame(ind))
      ! its index in the felix.hkl list is IhklsAll(ind), its g-vector is RInputHKLs(IhklsAll(ind))
      ! we track all requested reflections using ILiveList,
      ! this reflection is ILiveList(IhklsAll(ind)
      ! ILiveList = 0 never been simulated
      ! ILiveList = n it is current, this is the nth frame it has been found
      ! ILiveList = 666666, it's been simulated once and finished
      ! ILiveList = -n it is current, second time now, this is the nth frame
      ! ILiveList = -666666, it's been simulated twice, shouldn't be appearing again!

      ! only output an image if it has a deviation parameter |Sg|<0.05 [arbitrary? to test]
      IF(ABS(RDevPara(IhklsFrame(ind))).LT.0.05D0) THEN
        knd = knd + 1

        IF (ILiveList(IhklsAll(ind)).EQ.-666666) THEN! Shouldn't happen [to test!]
          WRITE(SPrintString,'(A22,I4,A3,I3,1X,I3,1X,I3)') "oops! third time for #",IhklsAll(ind)," : ",&
                       NINT(Rhkl(IhklsFrame(ind),:))
          CALL message(LS, SPrintString)
          CYCLE!at the moment just continue, but if it's a genuine error it will be IErr=1
        END IF
        IF(ILiveList(IhklsAll(ind)).LT.0) THEN!it's still live, second time
          ILiveList(IhklsAll(ind)) = ILiveList(IhklsAll(ind))-1! increment the counter
          WRITE(SPrintString,'(I2,A16,I2,A3,I3,1X,I3,1X,I3)') ILiveList(IhklsAll(ind)),&
                  " frames(*) for #",IhklsAll(ind)," : ",NINT(Rhkl(IhklsFrame(ind),:))
          IF(my_rank.EQ.0)PRINT*,SPrintString
        ELSEIF (ILiveList(IhklsAll(ind)).EQ.666666) THEN! This is the second time round, start counting negatively
          ILiveList(IhklsAll(ind)) = -1
          WRITE(SPrintString,'(A18,I3,A3,I3,1X,I3,1X,I3)') " second time for #",IhklsAll(ind),&
                  " : ",NINT(Rhkl(IhklsFrame(ind),:))
          IF(my_rank.EQ.0)PRINT*,SPrintString
        ELSE! ILiveList>0, so increment
          ILiveList(IhklsAll(ind)) = ILiveList(IhklsAll(ind))+1! increment the counter
          WRITE(SPrintString,'(I2,A13,I2,A3,I3,1X,I3,1X,I3)') ILiveList(IhklsAll(ind)),&
                  " frames for #",IhklsAll(ind)," : ",NINT(Rhkl(IhklsFrame(ind),:))
          IF(my_rank.EQ.0)PRINT*,SPrintString
        END IF
        
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
!        OPEN(UNIT=IChOutWIImage, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(fullpath)),&
!          FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=ISizeX*IByteSize)
!        IF(l_alert(IErr,"WriteIterationOutput","OPEN() output .bin file")) RETURN
!        ! we write each X-line, remembering indexing is [row,col]=[y,x]
!        DO jnd = 1,ISizeY
!          WRITE(IChOutWIImage,rec=jnd) RImageToWrite(jnd,:)
!        END DO
!        CLOSE(IChOutWIImage,IOSTAT=IErr) 
!        IF(l_alert(IErr,"WriteIterationOutput","CLOSE() output .bin file")) RETURN       
      ELSE!not strong enough to give an output, check if we need to finish off this reflection 
        IF (ILiveList(IhklsAll(ind)).EQ.0) THEN
          CYCLE! it has not been in any frame, ignore it
        ELSEIF (ABS(ILiveList(IhklsAll(ind))).NE.666666) THEN
          WRITE(SPrintString,'(A10,I3,A1,I3,1X,I3,1X,I3,A1,I3,A7)') "finished #",IhklsAll(ind),":",&
                NINT(Rhkl(IhklsFrame(ind),:)),",",ILiveList(IhklsAll(ind))," frames"
          IF(my_rank.EQ.0)PRINT*,SPrintString
          ILiveList(IhklsAll(ind)) = SIGN(666666,ILiveList(IhklsAll(ind)))!set the flag complete
          IhklsAll(ind) = -IhklsAll(ind)!negative value is a flag to close the output
        END IF
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
    CALL message(LS,"//////////")

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

