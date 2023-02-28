!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: 2.0
! Date: 19-12-2022
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
    USE IPARA, ONLY : ILN,ISizeX,ISizeY,IhklsFrame,INoOfHKLsFrame,IFrame,IOutputFLAG,&
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
    REAL(RKIND) :: RStartFrame,RIntegratedIntensity,RLorAngle,RgMag
    CHARACTER(10) :: hString,kString,lString
    CHARACTER(40) :: fString,SMiller,Shkl,SPrintString
    CHARACTER(200) :: path,filename,fullpath,OutputString

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
          Iflag = 0!a check to ensure a container has been found
          DO jnd=1,50! find an empty container RLACBED_(n)
            IF (ILACBEDFlag(jnd).EQ.0) THEN!we've found an empty container
              ILACBEDFlag(jnd) = 1! this container is now taken 
              ILACBEDList(IhklsAll(ind)) = jnd! links reflection and container
              Iflag = 1
              CALL SetupContainer(jnd,RImageToWrite,IErr)
              EXIT
            END IF
          END DO
          IF(Iflag.EQ.0)THEN! we have more reflections to track than containers to hold them, oops
            CALL message(LS, "Error: Too many reflections to track!")
            IErr = 1
            RETURN
          END IF
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
        !WRITE(SPrintString,'(A10,I3,A1,I3,1X,I3,1X,I3,A1,I7,A7)') "finished #",IhklsAll(ind),":",&
        !      NINT(Rhkl(IhklsFrame(ind),:)),",",ILiveList(IhklsAll(ind))," frames"
        !CALL message(LL, SPrintString)

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
        Shkl = TRIM(ADJUSTL(hString)) //","// TRIM(ADJUSTL(kString)) //","// TRIM(ADJUSTL(lString))
        SMiller = TRIM(ADJUSTL(hString)) // TRIM(ADJUSTL(kString)) // TRIM(ADJUSTL(lString))

        ! Lorentz factor
        ! x-direction in the microscope frame is [100] so g.x is just the x-component of g in the microscope frame
        RGvecM = Rhkl(IhklsFrame(ind),1)*RarVecM+Rhkl(IhklsFrame(ind),2)*RbrVecM+Rhkl(IhklsFrame(ind),3)*RcrVecM
        RgMag = SQRT(DOT_PRODUCT(RGvecM,RGvecM))
        RLorAngle = ACOS(RGvecM(1)/RgMag)*180.0D0/PI

        !--------------------------------------------------------------------
        ! loop over thicknesses for output
        DO lnd = 1,IThicknessCount
          ! Get the folder name for this thickness
          IThickness = NINT(RInitialThickness +(lnd-1)*RDeltaThickness)/10.0!in nm 
          WRITE(path,"(I4.4,A2)") IThickness,"nm"
          path = SChemicalFormula(1:ILN) // "_" // path ! This adds chemical formula to folder name
!DBG          IF(my_rank.EQ.0)PRINT*,"writing ",Smiller," to ",TRIM(ADJUSTL(path))

          !--------------------------------------------------------------------
          ! Write LACBED pattern to .bin file
          ! Make the path/filename e.g. 'GaAs_-2-2+0_10x100.bin'
          ! only happens if IOutputFLAG>1
          IF(IOutputFLAG.GE.2) THEN
            WRITE(fString,"(I,A1,I2)") SIZE(RTempImage,DIM=2),"x",ISizeY
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
          END IF

          !--------------------------------------------------------------------
          ! write rocking curve and integrated intensity
          ! Integrated intensity - always done
          fullpath = TRIM(ADJUSTL(path)) // "/IntegratedIntensities.txt"
          OPEN(UNIT=IChOutIhkl, ACTION='WRITE', POSITION='APPEND', STATUS= 'UNKNOWN', &
                  FILE=TRIM(ADJUSTL(fullpath)),IOSTAT=IErr)
          RIntegratedIntensity = 0.0D0
          Irow = NINT(HALF*REAL(ISizeY))! The central row        
          DO jnd=1,SIZE(RTempImage,DIM=2)
            RIntegratedIntensity = RIntegratedIntensity + RTempImage(Irow,jnd,lnd)
          END DO
          WRITE(fString,"(E12.5,A2,F8.4)") RIntegratedIntensity,", ",RLorAngle
          WRITE(IChOutIhkl,*) TRIM(ADJUSTL(Shkl))//", "// TRIM(ADJUSTL(fString))
          CLOSE(IChOutIhkl,IOSTAT=IErr)
          !Rocking curve - only for IOutputFLAG>0
          IF(IOutputFLAG.GE.1) THEN
            ! open rocking curve text file to append
            fullpath = TRIM(ADJUSTL(path)) // "/RockingCurves.txt"
            OPEN(UNIT=IChOutRC, ACTION='WRITE', POSITION='APPEND', STATUS='UNKNOWN', &
                 FILE=TRIM(ADJUSTL(fullpath)),IOSTAT=IErr)
            IF(l_alert(IErr,"Felixrefine","OPEN() RockingCurves.txt")) CALL abort
            ! Each frame in RTempImage adds ISizeX pixels, so its width in frames is
            IhklFrames = NINT(SIZE(RTempImage,DIM=2)/REAL(ISizeX))
            RStartFrame = IFrame - IhklFrames
            WRITE(IChOutRC,*)"#"!# to separate reflections
            WRITE(IChOutRC,*) TRIM(ADJUSTL(Shkl))!hkl
            WRITE(OutputString,"(F8.2)") RLorAngle!angle for Lorentz factor
            WRITE(IChOutRC,*) TRIM(ADJUSTL(OutputString))
            DO jnd=1,SIZE(RTempImage,DIM=2)
              WRITE(OutputString,"(F8.2,A2,E9.2)") RStartFrame+REAL(jnd)/REAL(ISizeX),", ", &
                    RTempImage(Irow,jnd,lnd)
              WRITE(IChOutRC,*) TRIM(ADJUSTL(OutputString))
            RIntegratedIntensity = RIntegratedIntensity + RTempImage(Irow,jnd,lnd)
            END DO
            WRITE(IChOutRC,*)
            CLOSE(IChOutRC,IOSTAT=IErr)
          END IF
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
            RLACBED_49,RLACBED_50,RLACBED_51,RLACBED_52,RLACBED_53,&
            RLACBED_54,RLACBED_55,RLACBED_56,RLACBED_57,RLACBED_58,&
            RLACBED_59,RLACBED_60,RLACBED_61,RLACBED_62,RLACBED_63,&
            RLACBED_64,RLACBED_65,RLACBED_66,RLACBED_67,RLACBED_68,&
            RLACBED_69,RLACBED_70,RLACBED_71,RLACBED_72,RLACBED_73,&
            RLACBED_74,RLACBED_75,RLACBED_76,RLACBED_77,RLACBED_78,&
            RLACBED_79,RLACBED_80,RLACBED_81,RLACBED_82,RLACBED_83,&
            RLACBED_84,RLACBED_85,RLACBED_86,RLACBED_87,RLACBED_88,&
            RLACBED_89,RLACBED_90,RLACBED_91,RLACBED_92,RLACBED_93,&
            RLACBED_94,RLACBED_95,RLACBED_96,RLACBED_97,RLACBED_98,&
            RLACBED_99,RLACBED_100
            
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
      ALLOCATE(RLACBED_30(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
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
      ALLOCATE(RLACBED_40(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
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
      ALLOCATE(RLACBED_50(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_50 = RImageToWrite
    CASE(51)
      ALLOCATE(RLACBED_51(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_51 = RImageToWrite
    CASE(52)
      ALLOCATE(RLACBED_52(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_52 = RImageToWrite
    CASE(53)
      ALLOCATE(RLACBED_53(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_53 = RImageToWrite
    CASE(54)
      ALLOCATE(RLACBED_54(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_54 = RImageToWrite
    CASE(55)
      ALLOCATE(RLACBED_55(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_55 = RImageToWrite
    CASE(56)
      ALLOCATE(RLACBED_56(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_56 = RImageToWrite
    CASE(57)
      ALLOCATE(RLACBED_57(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_57 = RImageToWrite
    CASE(58)
      ALLOCATE(RLACBED_58(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_58 = RImageToWrite
    CASE(59)
      ALLOCATE(RLACBED_59(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_59 = RImageToWrite
    CASE(60)
      ALLOCATE(RLACBED_60(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_60 = RImageToWrite
    CASE(61)
      ALLOCATE(RLACBED_61(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_61 = RImageToWrite
    CASE(62)
      ALLOCATE(RLACBED_62(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_62 = RImageToWrite
    CASE(63)
      ALLOCATE(RLACBED_63(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_63 = RImageToWrite
    CASE(64)
      ALLOCATE(RLACBED_64(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_64 = RImageToWrite
    CASE(65)
      ALLOCATE(RLACBED_65(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_65 = RImageToWrite
    CASE(66)
      ALLOCATE(RLACBED_66(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_66 = RImageToWrite
    CASE(67)
      ALLOCATE(RLACBED_67(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_67 = RImageToWrite
    CASE(68)
      ALLOCATE(RLACBED_68(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_68 = RImageToWrite
    CASE(69)
      ALLOCATE(RLACBED_69(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_69 = RImageToWrite
    CASE(70)
      ALLOCATE(RLACBED_70(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_70 = RImageToWrite
    CASE(71)
      ALLOCATE(RLACBED_71(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_71 = RImageToWrite
    CASE(72)
      ALLOCATE(RLACBED_72(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_72 = RImageToWrite
    CASE(73)
      ALLOCATE(RLACBED_73(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_73 = RImageToWrite
    CASE(74)
      ALLOCATE(RLACBED_74(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_74 = RImageToWrite
    CASE(75)
      ALLOCATE(RLACBED_75(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_75 = RImageToWrite
    CASE(76)
      ALLOCATE(RLACBED_76(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_76 = RImageToWrite
    CASE(77)
      ALLOCATE(RLACBED_77(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_77 = RImageToWrite
    CASE(78)
      ALLOCATE(RLACBED_78(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_78 = RImageToWrite
    CASE(79)
      ALLOCATE(RLACBED_79(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_79 = RImageToWrite
    CASE(80)
      ALLOCATE(RLACBED_80(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_80 = RImageToWrite
    CASE(81)
      ALLOCATE(RLACBED_81(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_81 = RImageToWrite
    CASE(82)
      ALLOCATE(RLACBED_82(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_82 = RImageToWrite
    CASE(83)
      ALLOCATE(RLACBED_83(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_83 = RImageToWrite
    CASE(84)
      ALLOCATE(RLACBED_84(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_84 = RImageToWrite
    CASE(85)
      ALLOCATE(RLACBED_85(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_85 = RImageToWrite
    CASE(86)
      ALLOCATE(RLACBED_86(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_86 = RImageToWrite
    CASE(87)
      ALLOCATE(RLACBED_87(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_87 = RImageToWrite
    CASE(88)
      ALLOCATE(RLACBED_88(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_88 = RImageToWrite
    CASE(89)
      ALLOCATE(RLACBED_89(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_89 = RImageToWrite
    CASE(90)
      ALLOCATE(RLACBED_90(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_90 = RImageToWrite
    CASE(91)
      ALLOCATE(RLACBED_91(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_91 = RImageToWrite
    CASE(92)
      ALLOCATE(RLACBED_92(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_92 = RImageToWrite
    CASE(93)
      ALLOCATE(RLACBED_93(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_93 = RImageToWrite
    CASE(94)
      ALLOCATE(RLACBED_94(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_94 = RImageToWrite
    CASE(95)
      ALLOCATE(RLACBED_95(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_95 = RImageToWrite
    CASE(96)
      ALLOCATE(RLACBED_96(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_96 = RImageToWrite
    CASE(97)
      ALLOCATE(RLACBED_97(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_97 = RImageToWrite
    CASE(98)
      ALLOCATE(RLACBED_98(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_98 = RImageToWrite
    CASE(99)
      ALLOCATE(RLACBED_99(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_99 = RImageToWrite
    CASE(100)
      ALLOCATE(RLACBED_100(ISizeY,ISizeX,IThicknessCount), STAT=IErr)
      RLACBED_100 = RImageToWrite

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
            RLACBED_49,RLACBED_50,RLACBED_51,RLACBED_52,RLACBED_53,&
            RLACBED_54,RLACBED_55,RLACBED_56,RLACBED_57,RLACBED_58,&
            RLACBED_59,RLACBED_60,RLACBED_61,RLACBED_62,RLACBED_63,&
            RLACBED_64,RLACBED_65,RLACBED_66,RLACBED_67,RLACBED_68,&
            RLACBED_69,RLACBED_70,RLACBED_71,RLACBED_72,RLACBED_73,&
            RLACBED_74,RLACBED_75,RLACBED_76,RLACBED_77,RLACBED_78,&
            RLACBED_79,RLACBED_80,RLACBED_81,RLACBED_82,RLACBED_83,&
            RLACBED_84,RLACBED_85,RLACBED_86,RLACBED_87,RLACBED_88,&
            RLACBED_89,RLACBED_90,RLACBED_91,RLACBED_92,RLACBED_93,&
            RLACBED_94,RLACBED_95,RLACBED_96,RLACBED_97,RLACBED_98,&
            RLACBED_99,RLACBED_100,RTempImage
            
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
    CASE(51)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_51,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_51,DIM=2),:) = RLACBED_51
      RTempImage(:,SIZE(RLACBED_51,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_51)
      ALLOCATE(RLACBED_51(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_51 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(52)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_52,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_52,DIM=2),:) = RLACBED_52
      RTempImage(:,SIZE(RLACBED_52,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_52)
      ALLOCATE(RLACBED_52(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_52 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(53)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_53,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_53,DIM=2),:) = RLACBED_53
      RTempImage(:,SIZE(RLACBED_53,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_53)
      ALLOCATE(RLACBED_53(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_53 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(54)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_54,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_54,DIM=2),:) = RLACBED_54
      RTempImage(:,SIZE(RLACBED_54,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_54)
      ALLOCATE(RLACBED_54(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_54 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(55)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_55,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_55,DIM=2),:) = RLACBED_55
      RTempImage(:,SIZE(RLACBED_55,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_55)
      ALLOCATE(RLACBED_55(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_55 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(56)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_56,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_56,DIM=2),:) = RLACBED_56
      RTempImage(:,SIZE(RLACBED_56,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_56)
      ALLOCATE(RLACBED_56(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_56 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(57)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_57,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_57,DIM=2),:) = RLACBED_57
      RTempImage(:,SIZE(RLACBED_57,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_57)
      ALLOCATE(RLACBED_57(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_57 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(58)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_58,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_58,DIM=2),:) = RLACBED_58
      RTempImage(:,SIZE(RLACBED_58,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_58)
      ALLOCATE(RLACBED_58(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_58 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(59)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_59,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_59,DIM=2),:) = RLACBED_59
      RTempImage(:,SIZE(RLACBED_59,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_59)
      ALLOCATE(RLACBED_59(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_59 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(60)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_60,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_60,DIM=2),:) = RLACBED_60
      RTempImage(:,SIZE(RLACBED_60,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_60)
      ALLOCATE(RLACBED_60(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_60 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(61)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_61,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_61,DIM=2),:) = RLACBED_61
      RTempImage(:,SIZE(RLACBED_61,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_61)
      ALLOCATE(RLACBED_61(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_61 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(62)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_62,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_62,DIM=2),:) = RLACBED_62
      RTempImage(:,SIZE(RLACBED_62,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_62)
      ALLOCATE(RLACBED_62(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_62 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(63)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_63,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_63,DIM=2),:) = RLACBED_63
      RTempImage(:,SIZE(RLACBED_63,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_63)
      ALLOCATE(RLACBED_63(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_63 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(64)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_64,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_64,DIM=2),:) = RLACBED_64
      RTempImage(:,SIZE(RLACBED_64,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_64)
      ALLOCATE(RLACBED_64(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_64 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(65)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_65,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_65,DIM=2),:) = RLACBED_65
      RTempImage(:,SIZE(RLACBED_65,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_65)
      ALLOCATE(RLACBED_65(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_65 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(66)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_66,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_66,DIM=2),:) = RLACBED_66
      RTempImage(:,SIZE(RLACBED_66,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_66)
      ALLOCATE(RLACBED_66(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_66 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(67)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_67,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_67,DIM=2),:) = RLACBED_67
      RTempImage(:,SIZE(RLACBED_67,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_67)
      ALLOCATE(RLACBED_67(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_67 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(68)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_68,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_68,DIM=2),:) = RLACBED_68
      RTempImage(:,SIZE(RLACBED_68,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_68)
      ALLOCATE(RLACBED_68(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_68 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(69)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_69,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_69,DIM=2),:) = RLACBED_69
      RTempImage(:,SIZE(RLACBED_69,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_69)
      ALLOCATE(RLACBED_69(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_69 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(70)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_70,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_70,DIM=2),:) = RLACBED_70
      RTempImage(:,SIZE(RLACBED_70,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_70)
      ALLOCATE(RLACBED_70(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_70 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(71)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_71,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_71,DIM=2),:) = RLACBED_71
      RTempImage(:,SIZE(RLACBED_71,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_71)
      ALLOCATE(RLACBED_71(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_71 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(72)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_72,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_72,DIM=2),:) = RLACBED_72
      RTempImage(:,SIZE(RLACBED_72,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_72)
      ALLOCATE(RLACBED_72(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_72 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(73)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_73,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_73,DIM=2),:) = RLACBED_73
      RTempImage(:,SIZE(RLACBED_73,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_73)
      ALLOCATE(RLACBED_73(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_73 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(74)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_74,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_74,DIM=2),:) = RLACBED_74
      RTempImage(:,SIZE(RLACBED_74,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_74)
      ALLOCATE(RLACBED_74(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_74 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(75)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_75,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_75,DIM=2),:) = RLACBED_75
      RTempImage(:,SIZE(RLACBED_75,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_75)
      ALLOCATE(RLACBED_75(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_75 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(76)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_76,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_76,DIM=2),:) = RLACBED_76
      RTempImage(:,SIZE(RLACBED_76,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_76)
      ALLOCATE(RLACBED_76(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_76 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(77)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_77,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_77,DIM=2),:) = RLACBED_77
      RTempImage(:,SIZE(RLACBED_77,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_77)
      ALLOCATE(RLACBED_77(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_77 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(78)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_78,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_78,DIM=2),:) = RLACBED_78
      RTempImage(:,SIZE(RLACBED_78,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_78)
      ALLOCATE(RLACBED_78(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_78 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(79)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_79,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_79,DIM=2),:) = RLACBED_79
      RTempImage(:,SIZE(RLACBED_79,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_79)
      ALLOCATE(RLACBED_79(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_79 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(80)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_80,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_80,DIM=2),:) = RLACBED_80
      RTempImage(:,SIZE(RLACBED_80,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_80)
      ALLOCATE(RLACBED_80(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_80 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(81)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_81,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_81,DIM=2),:) = RLACBED_81
      RTempImage(:,SIZE(RLACBED_81,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_81)
      ALLOCATE(RLACBED_81(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_81 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(82)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_82,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_82,DIM=2),:) = RLACBED_82
      RTempImage(:,SIZE(RLACBED_82,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_82)
      ALLOCATE(RLACBED_82(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_82 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(83)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_83,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_83,DIM=2),:) = RLACBED_83
      RTempImage(:,SIZE(RLACBED_83,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_83)
      ALLOCATE(RLACBED_83(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_83 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(84)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_84,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_84,DIM=2),:) = RLACBED_84
      RTempImage(:,SIZE(RLACBED_84,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_84)
      ALLOCATE(RLACBED_84(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_84 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(85)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_85,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_85,DIM=2),:) = RLACBED_85
      RTempImage(:,SIZE(RLACBED_85,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_85)
      ALLOCATE(RLACBED_85(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_85 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(86)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_86,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_86,DIM=2),:) = RLACBED_86
      RTempImage(:,SIZE(RLACBED_86,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_86)
      ALLOCATE(RLACBED_86(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_86 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(87)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_87,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_87,DIM=2),:) = RLACBED_87
      RTempImage(:,SIZE(RLACBED_87,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_87)
      ALLOCATE(RLACBED_87(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_87 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(88)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_88,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_88,DIM=2),:) = RLACBED_88
      RTempImage(:,SIZE(RLACBED_88,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_88)
      ALLOCATE(RLACBED_88(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_88 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(89)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_89,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_89,DIM=2),:) = RLACBED_89
      RTempImage(:,SIZE(RLACBED_89,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_89)
      ALLOCATE(RLACBED_89(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_89 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(90)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_90,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_90,DIM=2),:) = RLACBED_90
      RTempImage(:,SIZE(RLACBED_90,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_90)
      ALLOCATE(RLACBED_90(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_90 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(91)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_91,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_91,DIM=2),:) = RLACBED_91
      RTempImage(:,SIZE(RLACBED_91,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_91)
      ALLOCATE(RLACBED_91(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_91 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(92)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_92,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_92,DIM=2),:) = RLACBED_92
      RTempImage(:,SIZE(RLACBED_92,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_92)
      ALLOCATE(RLACBED_92(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_92 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(93)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_93,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_93,DIM=2),:) = RLACBED_93
      RTempImage(:,SIZE(RLACBED_93,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_93)
      ALLOCATE(RLACBED_93(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_93 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(94)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_94,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_94,DIM=2),:) = RLACBED_94
      RTempImage(:,SIZE(RLACBED_94,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_94)
      ALLOCATE(RLACBED_94(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_94 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(95)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_95,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_95,DIM=2),:) = RLACBED_95
      RTempImage(:,SIZE(RLACBED_95,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_95)
      ALLOCATE(RLACBED_95(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_95 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(96)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_96,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_96,DIM=2),:) = RLACBED_96
      RTempImage(:,SIZE(RLACBED_96,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_96)
      ALLOCATE(RLACBED_96(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_96 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(97)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_97,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_97,DIM=2),:) = RLACBED_97
      RTempImage(:,SIZE(RLACBED_97,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_97)
      ALLOCATE(RLACBED_97(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_97 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(98)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_98,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_98,DIM=2),:) = RLACBED_98
      RTempImage(:,SIZE(RLACBED_98,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_98)
      ALLOCATE(RLACBED_98(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_98 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(99)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_99,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_99,DIM=2),:) = RLACBED_99
      RTempImage(:,SIZE(RLACBED_99,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_99)
      ALLOCATE(RLACBED_99(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_99 = RTempImage
      DEALLOCATE(RTempImage)
    CASE(100)
      ALLOCATE(RTempImage(ISizeY,ISizeX+SIZE(RLACBED_100,DIM=2),IThicknessCount), STAT=IErr)
      IF(l_alert(IErr,"write_output","allocate RTempImage")) RETURN
      RTempImage(:,1:SIZE(RLACBED_100,DIM=2),:) = RLACBED_100
      RTempImage(:,SIZE(RLACBED_100,DIM=2)+1:,:) = RImageToWrite
      DEALLOCATE(RLACBED_100)
      ALLOCATE(RLACBED_100(ISizeY,SIZE(RTempImage,DIM=2),IThicknessCount), STAT=IErr)
      RLACBED_100 = RTempImage
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
            RLACBED_49,RLACBED_50,RLACBED_51,RLACBED_52,RLACBED_53,&
            RLACBED_54,RLACBED_55,RLACBED_56,RLACBED_57,RLACBED_58,&
            RLACBED_59,RLACBED_60,RLACBED_61,RLACBED_62,RLACBED_63,&
            RLACBED_64,RLACBED_65,RLACBED_66,RLACBED_67,RLACBED_68,&
            RLACBED_69,RLACBED_70,RLACBED_71,RLACBED_72,RLACBED_73,&
            RLACBED_74,RLACBED_75,RLACBED_76,RLACBED_77,RLACBED_78,&
            RLACBED_79,RLACBED_80,RLACBED_81,RLACBED_82,RLACBED_83,&
            RLACBED_84,RLACBED_85,RLACBED_86,RLACBED_87,RLACBED_88,&
            RLACBED_89,RLACBED_90,RLACBED_91,RLACBED_92,RLACBED_93,&
            RLACBED_94,RLACBED_95,RLACBED_96,RLACBED_97,RLACBED_98,&
            RLACBED_99,RLACBED_100,RTempImage

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
    CASE(51)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_51,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_51
      DEALLOCATE(RLACBED_51)
    CASE(52)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_52,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_52
      DEALLOCATE(RLACBED_52)
    CASE(53)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_53,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_53
      DEALLOCATE(RLACBED_53)
    CASE(54)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_54,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_54
      DEALLOCATE(RLACBED_54)
    CASE(55)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_55,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_55
      DEALLOCATE(RLACBED_55)
    CASE(56)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_56,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_56
      DEALLOCATE(RLACBED_56)
    CASE(57)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_57,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_57
      DEALLOCATE(RLACBED_57)
    CASE(58)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_58,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_58
      DEALLOCATE(RLACBED_58)
    CASE(59)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_59,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_59
      DEALLOCATE(RLACBED_59)
    CASE(60)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_60,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_60
      DEALLOCATE(RLACBED_60)
    CASE(61)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_61,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_61
      DEALLOCATE(RLACBED_61)
    CASE(62)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_62,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_62
      DEALLOCATE(RLACBED_62)
    CASE(63)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_63,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_63
      DEALLOCATE(RLACBED_63)
    CASE(64)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_64,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_64
      DEALLOCATE(RLACBED_64)
    CASE(65)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_65,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_65
      DEALLOCATE(RLACBED_65)
    CASE(66)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_66,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_66
      DEALLOCATE(RLACBED_66)
    CASE(67)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_67,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_67
      DEALLOCATE(RLACBED_67)
    CASE(68)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_68,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_68
      DEALLOCATE(RLACBED_68)
    CASE(69)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_69,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_69
      DEALLOCATE(RLACBED_69)
    CASE(70)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_70,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_70
      DEALLOCATE(RLACBED_70)
    CASE(71)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_71,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_71
      DEALLOCATE(RLACBED_71)
    CASE(72)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_72,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_72
      DEALLOCATE(RLACBED_72)
    CASE(73)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_73,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_73
      DEALLOCATE(RLACBED_73)
    CASE(74)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_74,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_74
      DEALLOCATE(RLACBED_74)
    CASE(75)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_75,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_75
      DEALLOCATE(RLACBED_75)
    CASE(76)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_76,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_76
      DEALLOCATE(RLACBED_76)
    CASE(77)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_77,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_77
      DEALLOCATE(RLACBED_77)
    CASE(78)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_78,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_78
      DEALLOCATE(RLACBED_78)
    CASE(79)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_79,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_79
      DEALLOCATE(RLACBED_79)
    CASE(80)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_80,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_80
      DEALLOCATE(RLACBED_80)
    CASE(81)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_81,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_81
      DEALLOCATE(RLACBED_81)
    CASE(82)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_82,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_82
      DEALLOCATE(RLACBED_82)
    CASE(83)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_83,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_83
      DEALLOCATE(RLACBED_83)
    CASE(84)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_84,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_84
      DEALLOCATE(RLACBED_84)
    CASE(85)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_85,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_85
      DEALLOCATE(RLACBED_85)
    CASE(86)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_86,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_86
      DEALLOCATE(RLACBED_86)
    CASE(87)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_87,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_87
      DEALLOCATE(RLACBED_87)
    CASE(88)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_88,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_88
      DEALLOCATE(RLACBED_88)
    CASE(89)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_89,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_89
      DEALLOCATE(RLACBED_89)
    CASE(90)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_90,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_90
      DEALLOCATE(RLACBED_90)
    CASE(91)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_91,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_91
      DEALLOCATE(RLACBED_91)
    CASE(92)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_92,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_92
      DEALLOCATE(RLACBED_92)
    CASE(93)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_93,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_93
      DEALLOCATE(RLACBED_93)
    CASE(94)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_94,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_94
      DEALLOCATE(RLACBED_94)
    CASE(95)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_95,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_95
      DEALLOCATE(RLACBED_95)
    CASE(96)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_96,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_96
      DEALLOCATE(RLACBED_96)
    CASE(97)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_97,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_97
      DEALLOCATE(RLACBED_97)
    CASE(98)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_98,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_98
      DEALLOCATE(RLACBED_98)
    CASE(99)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_99,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_99
      DEALLOCATE(RLACBED_99)
    CASE(100)
      ALLOCATE(RTempImage(ISizeY,SIZE(RLACBED_100,DIM=2),IThicknessCount), STAT=IErr)
      RTempImage = RLACBED_100
      DEALLOCATE(RLACBED_100)
    END SELECT

    RETURN

  END SUBROUTINE CloseContainer

END MODULE write_output_mod

