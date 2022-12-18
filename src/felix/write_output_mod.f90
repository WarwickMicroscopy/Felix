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
  !! Procedure-description: Writes output for iterations including simulated .bin files
  !! RockingCurves.txt and IntegratedIntensities.txt 
  !!
  !! Major-Authors: Richard Beanland (2022)
  !!
  SUBROUTINE WriteIterationOutput(IErr)

    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE IPARA, ONLY : ILN,ISizeX,ISizeY,IhklsFrame,INoOfHKLsFrame,IFrame,&
            IhklsAll,ILiveList,ILACBEDList,ILACBEDFlag,IByteSize,INFrames,IThicknessCount
    USE RPARA, ONLY : RInputHKLs,Rhkl, RImageSimi, RInitialThickness, RDeltaThickness,&
            RTempImage,RDevPara,RarVecM,RbrVecM,RcrVecM
    USE SPARA, ONLY : SChemicalFormula
    USE IChannels, ONLY : IChOutIM,IChOutRC,IChOutIhkl

    IMPLICIT NONE

    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND) :: IThickness,ind,jnd,knd,lnd,Iflag,Irow,IhklFrames
    REAL(RKIND),DIMENSION(ISizeY,ISizeX,IThicknessCount) :: RImageToWrite
    REAL(RKIND),DIMENSION(3) :: RGvecM
    REAL(RKIND) :: RStartFrame,RIntegratedIntensity,RLorAngle
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
      ! NEEDS an error catching mechanism when all containers are full
      IF(ABS(RDevPara(IhklsFrame(ind))).LT.0.05D0) THEN
        knd = knd + 1 ! counter for number of reflections found in this frame
        ILiveList(IhklsAll(ind)) = ILiveList(IhklsAll(ind))+1! increment counter for this reflection
        ! add the simulation into its output image
        IF(ILiveList(IhklsAll(ind)).EQ.1) THEN! new reflection, start up the output image
          DO jnd=1,50! find an empty container RLACBED_(n)
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
        ! Lorentz factor
        ! x-direction in the microscope frame is [100] so g.x is just the x-component of g in the microscope frame
        RGvecM = Rhkl(IhklsFrame(ind),1)*RarVecM+Rhkl(IhklsFrame(ind),2)*RbrVecM+Rhkl(IhklsFrame(ind),3)*RcrVecM
        RLorAngle = ACOS(RGvecM(1))

        !--------------------------------------------------------------------
        ! loop over thicknesses for output
        DO lnd = 1,IThicknessCount
          ! Get the folder name for this thickness
          IThickness = NINT(RInitialThickness +(lnd-1)*RDeltaThickness)/10.0!in nm 
          WRITE(path,"(I4.4,A2)") IThickness,"nm"
          path = SChemicalFormula(1:ILN) // "_" // path ! This adds chemical formula to folder name
!DBG          IF(my_rank.EQ.0)PRINT*,"writing ",Smiller," to ",TRIM(ADJUSTL(path))
          ! open rocking curve text file to append
          fullpath = TRIM(ADJUSTL(path)) // "/RockingCurves.txt"
          OPEN(UNIT=IChOutRC, ACTION='WRITE', POSITION='APPEND', STATUS='UNKNOWN', &
                  FILE=TRIM(ADJUSTL(fullpath)),IOSTAT=IErr)
          IF(l_alert(IErr,"Felixrefine","OPEN() RockingCurves.txt")) CALL abort
          ! open integrated intensities text file to append
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
          WRITE(IChOutRC,*) TRIM(ADJUSTL(SMiller))!hkl
          WRITE(OutputString,"(F8.5)") RLorAngle!angle for Lorentz factor
          WRITE(IChOutRC,*) TRIM(ADJUSTL(OutputString))
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
          WRITE(OutputString,"(F8.5)") RLorAngle!angle for Lorentz factor
          WRITE(IChOutIhkl,*) TRIM(ADJUSTL(hString)) //","// TRIM(ADJUSTL(kString)) //","// TRIM(ADJUSTL(lString)) &
                  //","// TRIM(ADJUSTL(fString)) // TRIM(ADJUSTL(OutputString))
          CLOSE(IChOutIhkl,IOSTAT=IErr)

          !--------------------------------------------------------------------
          ! Write LACBED pattern to .bin file
          ! Make the path/filename e.g. 'GaAs_-2-2+0_10x100.bin'
          WRITE(fString,"(I,A1,I)") SIZE(RTempImage,DIM=2),"x",ISizeY
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
    IF(IFrame.GT.9999)THEN
      WRITE(SPrintString,'(I3,A22,I5.5)') knd," reflections in Frame ",IFrame
    ELSEIF(IFrame.GT.999)THEN
      WRITE(SPrintString,'(I3,A22,I4.4)') knd," reflections in Frame ",IFrame
    ELSEIF(IFrame.GT.99)THEN
      WRITE(SPrintString,'(I3,A22,I3.3)') knd," reflections in Frame ",IFrame
    ELSEIF(IFrame.GT.9)THEN
      WRITE(SPrintString,'(I3,A22,I2.2)') knd," reflections in Frame ",IFrame
    ELSE
      WRITE(SPrintString,'(I3,A22,I1.1)') knd," reflections in Frame ",IFrame
    END IF
    CALL message(LS,SPrintString)
!    CALL message(LS,"//////////")

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
            RLACBED_19,RLACBED_20,RLACBED_21,RLACBED_22,RLACBED_23,&
            RLACBED_24,RLACBED_25,RLACBED_26,RLACBED_27,RLACBED_28,&
            RLACBED_29,RLACBED_30,RLACBED_31,RLACBED_32,RLACBED_33,&
            RLACBED_34,RLACBED_35,RLACBED_36,RLACBED_37,RLACBED_38,&
            RLACBED_39,RLACBED_40,RLACBED_41,RLACBED_42,RLACBED_43,&
            RLACBED_44,RLACBED_45,RLACBED_46,RLACBED_47,RLACBED_48,&
            RLACBED_49,RLACBED_50
            
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
    CASE(21)
      ALLOCATE(RLACBED_21(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_21 = RImageToWrite
    CASE(22)
      ALLOCATE(RLACBED_22(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_22 = RImageToWrite
    CASE(23)
      ALLOCATE(RLACBED_23(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_23 = RImageToWrite
    CASE(24)
      ALLOCATE(RLACBED_24(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_24 = RImageToWrite
    CASE(25)
      ALLOCATE(RLACBED_25(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_25 = RImageToWrite
    CASE(26)
      ALLOCATE(RLACBED_26(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_26 = RImageToWrite
    CASE(27)
      ALLOCATE(RLACBED_27(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_27 = RImageToWrite
    CASE(28)
      ALLOCATE(RLACBED_28(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_28 = RImageToWrite
    CASE(29)
      ALLOCATE(RLACBED_29(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_29 = RImageToWrite
    CASE(30)
      ALLOCATE(RLACBED_20(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_30 = RImageToWrite
    CASE(31)
      ALLOCATE(RLACBED_31(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_31 = RImageToWrite
    CASE(32)
      ALLOCATE(RLACBED_32(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_32 = RImageToWrite
    CASE(33)
      ALLOCATE(RLACBED_33(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_33 = RImageToWrite
    CASE(34)
      ALLOCATE(RLACBED_34(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_34 = RImageToWrite
    CASE(35)
      ALLOCATE(RLACBED_35(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_35 = RImageToWrite
    CASE(36)
      ALLOCATE(RLACBED_36(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_36 = RImageToWrite
    CASE(37)
      ALLOCATE(RLACBED_37(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_37 = RImageToWrite
    CASE(38)
      ALLOCATE(RLACBED_38(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_38 = RImageToWrite
    CASE(39)
      ALLOCATE(RLACBED_39(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_39 = RImageToWrite
    CASE(40)
      ALLOCATE(RLACBED_20(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_40 = RImageToWrite
    CASE(41)
      ALLOCATE(RLACBED_41(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_41 = RImageToWrite
    CASE(42)
      ALLOCATE(RLACBED_42(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_42 = RImageToWrite
    CASE(43)
      ALLOCATE(RLACBED_43(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_43 = RImageToWrite
    CASE(44)
      ALLOCATE(RLACBED_44(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_44 = RImageToWrite
    CASE(45)
      ALLOCATE(RLACBED_45(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_45 = RImageToWrite
    CASE(46)
      ALLOCATE(RLACBED_46(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_46 = RImageToWrite
    CASE(47)
      ALLOCATE(RLACBED_47(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_47 = RImageToWrite
    CASE(48)
      ALLOCATE(RLACBED_48(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_48 = RImageToWrite
    CASE(49)
      ALLOCATE(RLACBED_49(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_49 = RImageToWrite
    CASE(50)
      ALLOCATE(RLACBED_20(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_50 = RImageToWrite

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
    USE RPARA, ONLY : RLACBED_1,RLACBED_2,RLACBED_3,&
            RLACBED_4,RLACBED_5,RLACBED_6,RLACBED_7,RLACBED_8,&
            RLACBED_9,RLACBED_10,RLACBED_11,RLACBED_12,RLACBED_13,&
            RLACBED_14,RLACBED_15,RLACBED_16,RLACBED_17,RLACBED_18,&
            RLACBED_19,RLACBED_20,RLACBED_21,RLACBED_22,RLACBED_23,&
            RLACBED_24,RLACBED_25,RLACBED_26,RLACBED_27,RLACBED_28,&
            RLACBED_29,RLACBED_30,RLACBED_31,RLACBED_32,RLACBED_33,&
            RLACBED_34,RLACBED_35,RLACBED_36,RLACBED_37,RLACBED_38,&
            RLACBED_39,RLACBED_40,RLACBED_41,RLACBED_42,RLACBED_43,&
            RLACBED_44,RLACBED_45,RLACBED_46,RLACBED_47,RLACBED_48,&
            RLACBED_49,RLACBED_50,RTempImage
            
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
    CASE(21)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_21,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_21,DIM=2),:) = RLACBED_21
      RTempImage(:,SIZE(RLACBED_21,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_21)
      ALLOCATE(RLACBED_21(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_21 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(22)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_22,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_22,DIM=2),:) = RLACBED_22
      RTempImage(:,SIZE(RLACBED_22,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_22)
      ALLOCATE(RLACBED_22(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_22 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(23)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_23,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_23,DIM=2),:) = RLACBED_23
      RTempImage(:,SIZE(RLACBED_23,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_23)
      ALLOCATE(RLACBED_23(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_23 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(24)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_24,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_24,DIM=2),:) = RLACBED_24
      RTempImage(:,SIZE(RLACBED_24,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_24)
      ALLOCATE(RLACBED_24(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_24 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(25)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_25,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_25,DIM=2),:) = RLACBED_25
      RTempImage(:,SIZE(RLACBED_25,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_25)
      ALLOCATE(RLACBED_25(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_25 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(26)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_26,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_26,DIM=2),:) = RLACBED_26
      RTempImage(:,SIZE(RLACBED_26,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_26)
      ALLOCATE(RLACBED_26(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_26 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(27)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_27,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_27,DIM=2),:) = RLACBED_27
      RTempImage(:,SIZE(RLACBED_27,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_27)
      ALLOCATE(RLACBED_27(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_27 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(28)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_28,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_28,DIM=2),:) = RLACBED_28
      RTempImage(:,SIZE(RLACBED_28,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_28)
      ALLOCATE(RLACBED_28(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_28 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(29)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_29,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_29,DIM=2),:) = RLACBED_29
      RTempImage(:,SIZE(RLACBED_29,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_29)
      ALLOCATE(RLACBED_29(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_29 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(30)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_30,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_30,DIM=2),:) = RLACBED_30
      RTempImage(:,SIZE(RLACBED_30,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_30)
      ALLOCATE(RLACBED_30(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_30 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(31)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_31,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_31,DIM=2),:) = RLACBED_31
      RTempImage(:,SIZE(RLACBED_31,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_31)
      ALLOCATE(RLACBED_31(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_31 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(32)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_32,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_32,DIM=2),:) = RLACBED_32
      RTempImage(:,SIZE(RLACBED_32,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_32)
      ALLOCATE(RLACBED_32(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_32 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(33)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_33,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_33,DIM=2),:) = RLACBED_33
      RTempImage(:,SIZE(RLACBED_33,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_33)
      ALLOCATE(RLACBED_33(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_33 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(34)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_34,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_34,DIM=2),:) = RLACBED_34
      RTempImage(:,SIZE(RLACBED_34,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_34)
      ALLOCATE(RLACBED_34(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_34 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(35)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_35,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_35,DIM=2),:) = RLACBED_35
      RTempImage(:,SIZE(RLACBED_35,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_35)
      ALLOCATE(RLACBED_35(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_35 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(36)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_36,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_36,DIM=2),:) = RLACBED_36
      RTempImage(:,SIZE(RLACBED_36,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_36)
      ALLOCATE(RLACBED_36(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_36 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(37)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_37,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_37,DIM=2),:) = RLACBED_37
      RTempImage(:,SIZE(RLACBED_37,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_37)
      ALLOCATE(RLACBED_37(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_37 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(38)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_38,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_38,DIM=2),:) = RLACBED_38
      RTempImage(:,SIZE(RLACBED_38,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_38)
      ALLOCATE(RLACBED_38(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_38 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(39)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_39,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_39,DIM=2),:) = RLACBED_39
      RTempImage(:,SIZE(RLACBED_39,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_39)
      ALLOCATE(RLACBED_39(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_39 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(40)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_40,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_40,DIM=2),:) = RLACBED_40
      RTempImage(:,SIZE(RLACBED_40,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_40)
      ALLOCATE(RLACBED_40(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_40 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(41)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_41,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_41,DIM=2),:) = RLACBED_41
      RTempImage(:,SIZE(RLACBED_41,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_41)
      ALLOCATE(RLACBED_41(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_41 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(42)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_42,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_42,DIM=2),:) = RLACBED_42
      RTempImage(:,SIZE(RLACBED_42,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_42)
      ALLOCATE(RLACBED_42(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_42 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(43)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_43,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_43,DIM=2),:) = RLACBED_43
      RTempImage(:,SIZE(RLACBED_43,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_43)
      ALLOCATE(RLACBED_43(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_43 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(44)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_44,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_44,DIM=2),:) = RLACBED_44
      RTempImage(:,SIZE(RLACBED_44,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_44)
      ALLOCATE(RLACBED_44(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_44 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(45)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_45,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_45,DIM=2),:) = RLACBED_45
      RTempImage(:,SIZE(RLACBED_45,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_45)
      ALLOCATE(RLACBED_45(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_45 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(46)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_46,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_46,DIM=2),:) = RLACBED_46
      RTempImage(:,SIZE(RLACBED_46,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_46)
      ALLOCATE(RLACBED_46(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_46 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(47)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_47,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_47,DIM=2),:) = RLACBED_47
      RTempImage(:,SIZE(RLACBED_47,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_47)
      ALLOCATE(RLACBED_47(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_47 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(48)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_48,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_48,DIM=2),:) = RLACBED_48
      RTempImage(:,SIZE(RLACBED_48,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_48)
      ALLOCATE(RLACBED_48(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_48 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(49)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_49,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_49,DIM=2),:) = RLACBED_49
      RTempImage(:,SIZE(RLACBED_49,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_49)
      ALLOCATE(RLACBED_49(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_49 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(50)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_50,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_50,DIM=2),:) = RLACBED_50
      RTempImage(:,SIZE(RLACBED_50,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_50)
      ALLOCATE(RLACBED_50(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_50 = RTempImage
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
    USE RPARA, ONLY : RLACBED_1,RLACBED_2,RLACBED_3,&
            RLACBED_4,RLACBED_5,RLACBED_6,RLACBED_7,RLACBED_8,&
            RLACBED_9,RLACBED_10,RLACBED_11,RLACBED_12,RLACBED_13,&
            RLACBED_14,RLACBED_15,RLACBED_16,RLACBED_17,RLACBED_18,&
            RLACBED_19,RLACBED_20,RLACBED_21,RLACBED_22,RLACBED_23,&
            RLACBED_24,RLACBED_25,RLACBED_26,RLACBED_27,RLACBED_28,&
            RLACBED_29,RLACBED_30,RLACBED_31,RLACBED_32,RLACBED_33,&
            RLACBED_34,RLACBED_35,RLACBED_36,RLACBED_37,RLACBED_38,&
            RLACBED_39,RLACBED_40,RLACBED_41,RLACBED_42,RLACBED_43,&
            RLACBED_44,RLACBED_45,RLACBED_46,RLACBED_47,RLACBED_48,&
            RLACBED_49,RLACBED_50,RTempImage

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
    CASE(21)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_21,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_21
      DEALLOCATE(RLACBED_21)
    CASE(22)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_22,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_22
      DEALLOCATE(RLACBED_22)
    CASE(23)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_23,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_23
      DEALLOCATE(RLACBED_23)
    CASE(24)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_24,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_24
      DEALLOCATE(RLACBED_24)
    CASE(25)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_25,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_25
      DEALLOCATE(RLACBED_25)
    CASE(26)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_26,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_26
      DEALLOCATE(RLACBED_26)
    CASE(27)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_27,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_27
      DEALLOCATE(RLACBED_27)
    CASE(28)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_28,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_28
      DEALLOCATE(RLACBED_28)
    CASE(29)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_29,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_29
      DEALLOCATE(RLACBED_29)
    CASE(30)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_30,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_30
      DEALLOCATE(RLACBED_30)
    CASE(31)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_31,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_31
      DEALLOCATE(RLACBED_31)
    CASE(32)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_32,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_32
      DEALLOCATE(RLACBED_32)
    CASE(33)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_33,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_33
      DEALLOCATE(RLACBED_33)
    CASE(34)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_34,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_34
      DEALLOCATE(RLACBED_34)
    CASE(35)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_35,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_35
      DEALLOCATE(RLACBED_35)
    CASE(36)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_36,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_36
      DEALLOCATE(RLACBED_36)
    CASE(37)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_37,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_37
      DEALLOCATE(RLACBED_37)
    CASE(38)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_38,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_38
      DEALLOCATE(RLACBED_38)
    CASE(39)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_39,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_39
      DEALLOCATE(RLACBED_39)
    CASE(40)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_40,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_40
      DEALLOCATE(RLACBED_40)
    CASE(41)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_41,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_41
      DEALLOCATE(RLACBED_41)
    CASE(42)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_42,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_42
      DEALLOCATE(RLACBED_42)
    CASE(43)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_43,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_43
      DEALLOCATE(RLACBED_43)
    CASE(44)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_44,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_44
      DEALLOCATE(RLACBED_44)
    CASE(45)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_45,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_45
      DEALLOCATE(RLACBED_45)
    CASE(46)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_46,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_46
      DEALLOCATE(RLACBED_46)
    CASE(47)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_47,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_47
      DEALLOCATE(RLACBED_47)
    CASE(48)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_48,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_48
      DEALLOCATE(RLACBED_48)
    CASE(49)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_49,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_49
      DEALLOCATE(RLACBED_49)
    CASE(50)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_50,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_50
      DEALLOCATE(RLACBED_50)
    END SELECT

    RETURN

  END SUBROUTINE CloseContainer

END MODULE write_output_mod
