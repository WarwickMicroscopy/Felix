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
!! Writes output files, NB runs on core 0 only


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
  !! Major-Authors: 'kidwhizz' (2015), Richard Beanland (2016-22)
  !!
  SUBROUTINE WriteIterationOutput(IErr)

    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE IPARA, ONLY : ILN,ISizeX,ISizeY,IhklsFrame,INoOfHKLsFrame,IFrame,&
            IhklsAll,ILiveList,ILACBEDList,ILACBEDFlag,IByteSize,INFrames,IThicknessCount
    USE RPARA, ONLY : RInputHKLs,Rhkl, RImageSimi, RInitialThickness, RDeltaThickness,&
            RTempImage,RDevPara,RLACBED_1,RLACBED_2,RLACBED_3,&
            RLACBED_4,RLACBED_5,RLACBED_6,RLACBED_7,RLACBED_8,&
            RLACBED_9,RLACBED_10,RLACBED_11,RLACBED_12,RLACBED_13,&
            RLACBED_14,RLACBED_15,RLACBED_16,RLACBED_17,RLACBED_18,&
            RLACBED_19,RLACBED_20
    USE SPARA, ONLY : SChemicalFormula
    USE IChannels, ONLY : IChOutIM,IChOutRC,IChOutIhkl

    IMPLICIT NONE

    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND) :: IThickness,ind,jnd,knd,lnd,Iflag,Irow,IhklFrames
    REAL(RKIND),DIMENSION(ISizeY,ISizeX,IThicknessCount) :: RImageToWrite
    REAL(RKIND) :: RStartFrame,RIntegratedIntensity
    CHARACTER(10) :: hString,kString,lString
    CHARACTER(40) :: fString,SMiller
    CHARACTER(40) :: SPrintString
    CHARACTER(200) :: path,filename,fullpath
    CHARACTER(:), ALLOCATABLE :: OutputString

    IErr=0 

    knd = 0
    DO ind = 1,INoOfHKLsFrame!the number of reflections in both felix.hkl and the beam pool
      ! this reflection is number ind in RImageSimi
      RImageToWrite = RImageSimi(:,:,ind,:)
      ! its index in the beam pool is IhklsFrame(ind), its g-vector is Rhkl(IhklsFrame(ind))
      ! its index in the felix.hkl list is IhklsAll(ind), its g-vector is RInputHKLs(IhklsAll(ind))
      ! we track all requested reflections using ILiveList,
      ! this reflection is ILiveList(IhklsAll(ind))
      ! ILiveList = 0 not active
      ! ILiveList = n it is current, this is the nth frame it has been found

      !For the final frame: force saving of any active reflections 
      IF (IFrame.EQ.INFrames .AND. ABS(ILiveList(IhklsAll(ind))).NE.0) THEN
        RDevPara(IhklsFrame(ind)) = 10.0D0! artificially set deviation parameter to a big number
      END IF

      ! only output an image if it has a deviation parameter |Sg|<0.05 [arbitrary? to test]
      IF(ABS(RDevPara(IhklsFrame(ind))).LT.0.05D0) THEN
        knd = knd + 1 ! counter for number of reflections found in this frame
        ILiveList(IhklsAll(ind)) = ILiveList(IhklsAll(ind))+1! increment counter for this reflection
        ! add the simulation into its output image
        IF(ILiveList(IhklsAll(ind)).EQ.1) THEN! new reflection, start up the output image
          DO jnd=1,20! find an empty container RLACBED_(n)
            IF (ILACBEDFlag(jnd).EQ.0) THEN
              ILACBEDFlag(jnd) = 1! this container is now taken 
              ILACBEDList(IhklsAll(ind)) = jnd! links reflection and container
              CALL SetupContainer(jnd,RImageToWrite,IErr)
              EXIT
            END IF
          END DO
        ELSE!its a continuation, find which container we're using and append the image
          CALL AppendContainer(ILACBEDList(IhklsAll(ind)),RImageToWrite,IErr)
        END IF  
        WRITE(SPrintString,'(I4,A13,I4,A3,I3,1X,I3,1X,I3,A1,I4)') ILiveList(IhklsAll(ind)),&
             " frames for #",IhklsAll(ind)," : ",NINT(Rhkl(IhklsFrame(ind),:)),&
             "|",ILACBEDList(IhklsAll(ind))
        CALL message(LL, SPrintString)
      ELSE! |Sg| is large, check if we need to finish off this reflection 
        IF (ILiveList(IhklsAll(ind)).EQ.0) CYCLE! it is not active, ignore it
        ! It is active, put the finished LACBED into RTempImage, extract rocking curve and save
        CALL CloseContainer(ILACBEDList(IhklsAll(ind)),IErr)
        WRITE(SPrintString,'(A10,I3,A1,I3,1X,I3,1X,I3,A1,I7,A7)') "finished #",IhklsAll(ind),":",&
              NINT(Rhkl(IhklsFrame(ind),:)),",",ILiveList(IhklsAll(ind))," frames"
        CALL message(LL, SPrintString)

        !--------------------------------------------------------------------
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
        SMiller = TRIM(ADJUSTL(hString)) // TRIM(ADJUSTL(kString)) // TRIM(ADJUSTL(lString))

        !--------------------------------------------------------------------
        ! loop over thicknesses for output
        DO lnd=1,IThicknessCount
          ! Get the folder name for this thickness
          IThickness = NINT(RInitialThickness +(lnd-1)*RDeltaThickness)/10.0!in nm 
          WRITE(path,"(I4.4,A2)") IThickness,"nm"
          path = SChemicalFormula(1:ILN) // "_" // path ! This adds chemical formula to folder name
          IF(my_rank.EQ.0)PRINT*,"writing ",Smiller," to ",TRIM(ADJUSTL(path))
          ! open rocking curve text file to append
          fullpath = TRIM(ADJUSTL(path)) // "/RockingCurves.txt"
          OPEN(UNIT=IChOutRC, ACTION='WRITE', POSITION='APPEND', STATUS='UNKNOWN', &
                  FILE=TRIM(ADJUSTL(fullpath)),IOSTAT=IErr)
          IF(l_alert(IErr,"Felixrefine","OPEN() RockingCurves.txt")) CALL abort
          ! open integrated intensities text file
          fullpath = TRIM(ADJUSTL(path)) // "/IntegratedIntensities.txt"
          OPEN(UNIT=IChOutIhkl, ACTION='WRITE', POSITION='APPEND', STATUS= 'UNKNOWN', &
                  FILE=TRIM(ADJUSTL(fullpath)),IOSTAT=IErr)

          !--------------------------------------------------------------------
          ! write rocking curve and integrated intensity
          ! Each frame in RTempImage adds ISizeX pixels, so its width in frames is
          IhklFrames = NINT(SIZE(RTempImage,DIM=2)/REAL(ISizeX))
          RStartFrame = IFrame - IhklFrames
          ! The central row
          Irow = NINT(HALF*REAL(ISizeY))        
          RIntegratedIntensity = 0.0D0
          ! Rocking curve
          WRITE(IChOutRC,*) TRIM(ADJUSTL(SMiller))
          DO jnd=1,SIZE(RTempImage,DIM=2)
            WRITE(fString,"(F8.2)") RStartFrame+REAL(jnd)/REAL(ISizeX)
            OutputString = TRIM(ADJUSTL(fString)) // ", "
            WRITE(fString,"(F8.5)") RTempImage(Irow,jnd,lnd)
            OutputString = OutputString // TRIM(ADJUSTL(fString))
            WRITE(IChOutRC,*) TRIM(ADJUSTL(OutputString))
            RIntegratedIntensity = RIntegratedIntensity + RTempImage(Irow,jnd,lnd)
          END DO
          WRITE(IChOutRC,*)!blank line to separate reflections
          CLOSE(IChOutRC,IOSTAT=IErr)
          ! Integrated intensity
          WRITE(fString,"(F9.5)") RIntegratedIntensity
          WRITE(IChOutIhkl,*) TRIM(ADJUSTL(SMiller)) // ":  " // TRIM(ADJUSTL(fString))
          CLOSE(IChOutIhkl,IOSTAT=IErr)

          !--------------------------------------------------------------------
          ! Write LACBED pattern to .bin file
          ! Make the path/filename e.g. 'GaAs_-2-2+0_10x100.bin'
          WRITE(fString,"(I4,A1,I2)") SIZE(RTempImage,DIM=2),"x",ISizeY
          filename = SChemicalFormula(1:ILN) // "_" // TRIM(ADJUSTL(SMiller)) // "_" &
                    // TRIM(ADJUSTL(fString)) // ".bin"
          fullpath = TRIM(ADJUSTL(path))//"/"//TRIM(ADJUSTL(filename))
          CALL message ( LL, dbg6, fullpath )
          OPEN(UNIT=IChOutIM, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(fullpath)),&
              FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=SIZE(RTempImage,DIM=2)*IByteSize)
          IF(l_alert(IErr,"WriteIterationOutput","OPEN() output .bin file")) RETURN
          ! we write each X-line, remembering indexing is [row,col]=[y,x]
          DO jnd = 1,ISizeY
            WRITE(IChOutIM,rec=jnd) RTempImage(jnd,:,lnd)
          END DO
          CLOSE(IChOutIM,IOSTAT=IErr) 
          IF(l_alert(IErr,"WriteIterationOutput","CLOSE() output .bin file")) RETURN       
        END DO

        !Tidy up and reset flags
        DEALLOCATE(RTempImage)
        ILACBEDFlag(ILACBEDList(IhklsAll(ind))) = 0! container is available again
        ILiveList(IhklsAll(ind)) = 0! the reflection is no longer active
      END IF
    END DO

    !--------------------------------------------------------------------
    ! output how many useful reflections have been calculated
    WRITE(SPrintString,'(A6,I3,A22,I5)') "Found ",knd," reflections in Frame ",IFrame
    CALL message(LS,SPrintString)
    CALL message(LS,"//////////")

    RETURN  
    
  END SUBROUTINE WriteIterationOutput

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: allocates a first-time container for a LACBED output 
  !!
  !! Major-Authors: Richard Beanland (2022)
  !!
  SUBROUTINE SetupContainer(jnd,RImageToWrite,IErr)

    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE IPARA, ONLY : ISizeX,ISizeY,IThicknessCount
    USE RPARA, ONLY : RLACBED_1,RLACBED_2,RLACBED_3,&
            RLACBED_4,RLACBED_5,RLACBED_6,RLACBED_7,RLACBED_8,&
            RLACBED_9,RLACBED_10,RLACBED_11,RLACBED_12,RLACBED_13,&
            RLACBED_14,RLACBED_15,RLACBED_16,RLACBED_17,RLACBED_18,&
            RLACBED_19,RLACBED_20
            
   IMPLICIT NONE
            
   REAL(RKIND),DIMENSION(ISizeY,ISizeX,IThicknessCount) :: RImageToWrite
   INTEGER(IKIND) :: jnd,IErr

    SELECT CASE(jnd)
    CASE(1)
      ALLOCATE(RLACBED_1(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_1 = RImageToWrite
    CASE(2)
      ALLOCATE(RLACBED_2(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_2 = RImageToWrite
    CASE(3)
      ALLOCATE(RLACBED_3(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_3 = RImageToWrite
    CASE(4)
      ALLOCATE(RLACBED_4(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_4 = RImageToWrite
    CASE(5)
      ALLOCATE(RLACBED_5(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_5 = RImageToWrite
    CASE(6)
      ALLOCATE(RLACBED_6(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_6 = RImageToWrite
    CASE(7)
      ALLOCATE(RLACBED_7(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_7 = RImageToWrite
    CASE(8)
      ALLOCATE(RLACBED_8(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_8 = RImageToWrite
    CASE(9)
      ALLOCATE(RLACBED_9(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_9 = RImageToWrite
    CASE(10)
      ALLOCATE(RLACBED_10(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_10 = RImageToWrite
    CASE(11)
      ALLOCATE(RLACBED_11(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_11 = RImageToWrite
    CASE(12)
      ALLOCATE(RLACBED_12(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_12 = RImageToWrite
    CASE(13)
      ALLOCATE(RLACBED_13(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_13 = RImageToWrite
    CASE(14)
      ALLOCATE(RLACBED_14(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_14 = RImageToWrite
    CASE(15)
      ALLOCATE(RLACBED_15(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_15 = RImageToWrite
    CASE(16)
      ALLOCATE(RLACBED_16(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_16 = RImageToWrite
    CASE(17)
      ALLOCATE(RLACBED_17(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_17 = RImageToWrite
    CASE(18)
      ALLOCATE(RLACBED_18(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_18 = RImageToWrite
    CASE(19)
      ALLOCATE(RLACBED_19(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_19 = RImageToWrite
    CASE(20)
      ALLOCATE(RLACBED_20(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_20 = RImageToWrite
    END SELECT

    RETURN   
  END SUBROUTINE SetupContainer

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: appends image onto an existing LACBED container 
  !!
  !! Major-Authors: Richard Beanland (2022)
  !!
  SUBROUTINE AppendContainer(jnd,RImageToWrite,IErr)

    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE IPARA, ONLY : ISizeX,ISizeY,IThicknessCount
    USE RPARA, ONLY : RTempImage, RLACBED_1,RLACBED_2,RLACBED_3,&
            RLACBED_4,RLACBED_5,RLACBED_6,RLACBED_7,RLACBED_8,&
            RLACBED_9,RLACBED_10,RLACBED_11,RLACBED_12,RLACBED_13,&
            RLACBED_14,RLACBED_15,RLACBED_16,RLACBED_17,RLACBED_18,&
            RLACBED_19,RLACBED_20
            
    IMPLICIT NONE
            
    REAL(RKIND),DIMENSION(ISizeY,ISizeX,IThicknessCount) :: RImageToWrite
    INTEGER(IKIND) :: jnd,IErr
   
    SELECT CASE(jnd)
    CASE(1)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_1,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_1,DIM=2),:) = RLACBED_1
      RTempImage(:,SIZE(RLACBED_1,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_1)
      ALLOCATE(RLACBED_1(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_1 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(2)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_2,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_2,DIM=2),:) = RLACBED_2
      RTempImage(:,SIZE(RLACBED_2,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_2)
      ALLOCATE(RLACBED_2(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_2 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(3)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_3,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_3,DIM=2),:) = RLACBED_3
      RTempImage(:,SIZE(RLACBED_3,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_3)
      ALLOCATE(RLACBED_3(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_3 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(4)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_4,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_4,DIM=2),:) = RLACBED_4
      RTempImage(:,SIZE(RLACBED_4,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_4)
      ALLOCATE(RLACBED_4(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_4 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(5)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_5,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_5,DIM=2),:) = RLACBED_5
      RTempImage(:,SIZE(RLACBED_5,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_5)
      ALLOCATE(RLACBED_5(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_5 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(6)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_6,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_6,DIM=2),:) = RLACBED_6
      RTempImage(:,SIZE(RLACBED_6,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_6)
      ALLOCATE(RLACBED_6(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_6 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(7)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_7,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_7,DIM=2),:) = RLACBED_7
      RTempImage(:,SIZE(RLACBED_7,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_7)
      ALLOCATE(RLACBED_7(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_7 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(8)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_8,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_8,DIM=2),:) = RLACBED_8
      RTempImage(:,SIZE(RLACBED_8,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_8)
      ALLOCATE(RLACBED_8(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_8 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(9)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_9,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_9,DIM=2),:) = RLACBED_9
      RTempImage(:,SIZE(RLACBED_9,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_9)
      ALLOCATE(RLACBED_9(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_9 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(10)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_10,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_10,DIM=2),:) = RLACBED_10
      RTempImage(:,SIZE(RLACBED_10,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_10)
      ALLOCATE(RLACBED_10(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_10 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(11)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_11,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_11,DIM=2),:) = RLACBED_11
      RTempImage(:,SIZE(RLACBED_11,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_11)
      ALLOCATE(RLACBED_11(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_11 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(12)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_12,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_12,DIM=2),:) = RLACBED_12
      RTempImage(:,SIZE(RLACBED_12,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_12)
      ALLOCATE(RLACBED_12(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_12 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(13)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_13,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_13,DIM=2),:) = RLACBED_13
      RTempImage(:,SIZE(RLACBED_13,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_13)
      ALLOCATE(RLACBED_13(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_13 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(14)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_14,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_14,DIM=2),:) = RLACBED_14
      RTempImage(:,SIZE(RLACBED_14,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_14)
      ALLOCATE(RLACBED_14(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_14 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(15)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_15,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_15,DIM=2),:) = RLACBED_15
      RTempImage(:,SIZE(RLACBED_15,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_15)
      ALLOCATE(RLACBED_15(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_15 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(16)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_16,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_16,DIM=2),:) = RLACBED_16
      RTempImage(:,SIZE(RLACBED_16,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_16)
      ALLOCATE(RLACBED_16(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_16 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(17)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_17,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_17,DIM=2),:) = RLACBED_17
      RTempImage(:,SIZE(RLACBED_17,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_17)
      ALLOCATE(RLACBED_17(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_17 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(18)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_18,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_18,DIM=2),:) = RLACBED_18
      RTempImage(:,SIZE(RLACBED_18,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_18)
      ALLOCATE(RLACBED_18(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_18 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(19)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_19,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_19,DIM=2),:) = RLACBED_19
      RTempImage(:,SIZE(RLACBED_19,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_19)
      ALLOCATE(RLACBED_19(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_19 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(20)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_20,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_20,DIM=2),:) = RLACBED_20
      RTempImage(:,SIZE(RLACBED_20,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_20)
      ALLOCATE(RLACBED_20(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_20 = RTempImage
      DEALLOCATE(RTempImage)
    END SELECT

    RETURN
  END SUBROUTINE AppendContainer

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: appends image onto an existing LACBED container 
  !!
  !! Major-Authors: Richard Beanland (2022)
  !!
  SUBROUTINE CloseContainer(jnd,IErr)

    USE MyNumbers
    USE message_mod

    ! global inputs
    USE IPARA, ONLY : ISizeY,IThicknessCount
    USE RPARA, ONLY : RTempImage,RLACBED_1,RLACBED_2,RLACBED_3,&
            RLACBED_4,RLACBED_5,RLACBED_6,RLACBED_7,RLACBED_8,&
            RLACBED_9,RLACBED_10,RLACBED_11,RLACBED_12,RLACBED_13,&
            RLACBED_14,RLACBED_15,RLACBED_16,RLACBED_17,RLACBED_18,&
            RLACBED_19,RLACBED_20

    IMPLICIT NONE

    INTEGER(IKIND) :: jnd,IErr

    SELECT CASE(jnd)

    CASE(1)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_1,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_1
      DEALLOCATE(RLACBED_1)
    CASE(2)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_2,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_2
      DEALLOCATE(RLACBED_2)
    CASE(3)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_3,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_3
      DEALLOCATE(RLACBED_3)
    CASE(4)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_4,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_4
      DEALLOCATE(RLACBED_4)
    CASE(5)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_5,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_5
      DEALLOCATE(RLACBED_5)
    CASE(6)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_6,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_6
      DEALLOCATE(RLACBED_6)
    CASE(7)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_7,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_7
      DEALLOCATE(RLACBED_7)
    CASE(8)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_8,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_8
      DEALLOCATE(RLACBED_8)
    CASE(9)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_9,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_9
      DEALLOCATE(RLACBED_9)
    CASE(10)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_10,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_10
      DEALLOCATE(RLACBED_10)
    CASE(11)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_11,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_11
      DEALLOCATE(RLACBED_11)
    CASE(12)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_12,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_12
      DEALLOCATE(RLACBED_12)
    CASE(13)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_13,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_13
      DEALLOCATE(RLACBED_13)
    CASE(14)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_14,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_14
      DEALLOCATE(RLACBED_14)
    CASE(15)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_15,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_15
      DEALLOCATE(RLACBED_15)
    CASE(16)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_16,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_16
      DEALLOCATE(RLACBED_16)
    CASE(17)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_17,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_17
      DEALLOCATE(RLACBED_17)
    CASE(18)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_18,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_18
      DEALLOCATE(RLACBED_18)
    CASE(19)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_19,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_19
      DEALLOCATE(RLACBED_19)
    CASE(20)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_20,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_20
      DEALLOCATE(RLACBED_20)
    END SELECT

    RETURN

  END SUBROUTINE CloseContainer

END MODULE write_output_mod

