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
  SUBROUTINE WriteIterationOutput(IFrame, IErr)

    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE IPARA, ONLY : ILN,ISizeX,ISizeY,IhklsFrame,INoOfHKLsFrame,IOutputFLAG,&
            IhklsAll,ILiveList,ILACBEDList,ILACBEDFlag,IByteSize,INFrames,IThicknessCount
    USE RPARA, ONLY : RInputHKL,Rhkl, RImageSimi, RInitialThickness, RDeltaThickness,&
            RTempImage,RDevC,RarVecM,RbrVecM,RcrVecM
    USE SPARA, ONLY : SChemicalFormula
    USE IChannels, ONLY : IChOutIM,IChOutRC,IChOutIhkl

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(IN) :: IFrame
    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND) :: IThickness,ind,jnd,knd,lnd,Iflag,Irow,IhklFrames
    REAL(RKIND),DIMENSION(ISizeY,ISizeX,IThicknessCount) :: RImageToWrite
    REAL(RKIND),DIMENSION(ITHREE) :: RGvecM
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
      ! its index in the felix.hkl list is IhklsAll(ind), its g-vector is RInputHKL(IhklsAll(ind))
      ! we track all requested reflections using ILiveList,
      ! this reflection is ILiveList(IhklsAll(ind))
      ! ILiveList = 0 not active
      ! ILiveList = n it is current, this is the nth frame it has been found

      !For the final frame: force saving of any active reflections 
      IF (IFrame.EQ.INFrames .AND. ABS(ILiveList(IhklsAll(ind))).NE.0) THEN
        RDevC(IhklsFrame(ind)) = 10.0D0! artificially set deviation parameter to a big number
      END IF

      ! only output an image if it has a deviation parameter |Sg|<0.05 [arbitrary? to test]
      ! NEEDS an error catching mechanism when all containers are full
      IF(ABS(RDevC(IhklsFrame(ind))).LT.0.05D0) THEN
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
              !CALL SetupContainer(jnd,RImageToWrite,IErr)
              EXIT
            END IF
          END DO
          IF(Iflag.EQ.0)THEN! we have more reflections to track than containers to hold them, oops
            CALL message(LS, "Error: Too many reflections to track!")
            IErr = 1
            RETURN
          END IF
        !ELSE!its a continuation, find which container we're using and append the image
          !CALL AppendContainer(ILACBEDList(IhklsAll(ind)),RImageToWrite,IErr)
        END IF  
        WRITE(SPrintString,'(I4,A13,I4,A3,I3,1X,I3,1X,I3,A1,I4)') ILiveList(IhklsAll(ind)),&
             " frames for #",IhklsAll(ind)," : ",NINT(Rhkl(IhklsFrame(ind),:)),&
             "|",ILACBEDList(IhklsAll(ind))
        CALL message(LL, SPrintString)
      ELSE! |Sg| is large, check if we need to finish off this reflection 
        IF (ILiveList(IhklsAll(ind)).EQ.0) CYCLE! it is not active, ignore it
        ! It is active, put the finished LACBED into RTempImage, extract rocking curve and save
        !CALL CloseContainer(ILACBEDList(IhklsAll(ind)),IErr)
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

 

END MODULE write_output_mod

