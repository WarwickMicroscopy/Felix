!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: :VERSION: RB_coord / 1.15 /
! Date:    :DATE: 16-01-2019
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build:   :BUILD: Mode F: test different lattice types" 
! Author:  :AUTHOR: r.beanland
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
  PUBLIC :: WriteIterationOutputWrapper, WriteIterationOutput, WriteOutVariables, &
        NormaliseExperimentalImagesAndWriteOut,WriteDifferenceImages,UncertBrak

  CONTAINS

  !>
  !!Procedure-description: gives an uncertainty in brackets
  !!
  !! Author: Richard Beanland (2019)
  !!
  SUBROUTINE UncertBrak(Rval,Rerr,Sout,IErr)

    USE MyNumbers
    USE MyMPI


    IMPLICIT NONE
    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND) :: Iord
    REAL(RKIND), INTENT(INOUT) :: Rval,Rerr
    REAL(RKIND) :: Rmag,RtruncErr,RtruncVal,RerrI
    CHARACTER(14), INTENT(OUT) :: Sout
    CHARACTER(10) :: Sval,Serr,Sformat

    !infinity and NaN check
    IF (ABS(Rval)-1.GT.ABS(Rval).OR.ABS(Rval).NE.ABS(Rval)) THEN
      IErr=1
!      Sout=""
      RETURN
    END IF
    IF (ABS(Rerr)-1.GT.ABS(Rerr).OR.ABS(Rerr).NE.ABS(Rerr)) THEN
      IErr=1
!      Sout=""
      RETURN
    END IF

    !Zero check
    IF (ABS(Rerr).GT.TINY) THEN
      !Make error positive if it isn't already
      IF (Rerr.LT.ZERO) Rerr=-Rerr
      !Check error is in expected limits
      IF (Rerr.GE.HUNDRED.OR.CEILING(Rerr*10000.0).LT.ONE) THEN
        IErr=1
        RETURN
      END IF

      !Find order of magnitude for error, Iord
      Rmag=TEN
      Iord=1
      DO WHILE (Rerr.LT.Rmag)
        Rmag=Rmag/TEN
        Iord=Iord-1
      END DO

      ! make an number RerrI corresponding to Rerr between 1 and 10 if it's a decimal
      IF (Iord.LT.0) THEN
        RerrI=Rerr/Rmag
      !if it is already an integer then no change
      ELSE
        RerrI=Rerr
        !round to 1 s.f. if needed
        IF (Rerr.GT.15.0) RerrI=TEN*(NINT(Rerr/TEN))
      END IF

      !Make the error string, use a 2-digit error for values below 1.5
      !Trim the value to the correct number of decimals using Sformat
      IF (Iord.LT.0) THEN!subdecimal error
        !make the error string and the format for the value string
        IF (RerrI.GE.9.5) THEN!we round up to 1
          Serr="(1)"
          Iord=Iord+1
        ELSEIF (RerrI.GE.1.05.AND.RerrI.LE.1.5) THEN!we need 2 digits
          WRITE(Serr,FMT='(A1,I2,A1)')  "(",NINT(RerrI*TEN),")"
          Iord=Iord-1
        ELSE
          WRITE(Serr,FMT='(A1,I1,A1)')  "(",NINT(RerrI),")"
        END IF
        WRITE(Sformat,FMT='(A4,I1,A1)') "(F8." , -Iord , ")"
      ELSE!1 to 99
        IF (RerrI.GE.95) THEN!we round up to 100
          Serr="(100.)"
        ELSEIF(Rerr.GE.1.05.AND.Rerr.LE.1.5) THEN!include a decimal
          WRITE(Serr,FMT='(A1,F3.1,A1)')  "(",RerrI,")"
          WRITE(Sformat,*) "(F8.1)"
        ELSE
          IF (Rerr.GE.TEN) THEN!double digit with a . to show it's >0
            WRITE(Serr,FMT='(A1,I2,A2)')  "(",NINT(RerrI),".)"
          ELSE!single digit with a . to show it's >0
            WRITE(Serr,FMT='(A1,I1,A2)')  "(",NINT(RerrI),".)"
          END IF
          WRITE(Sformat,*) "(F8.0)"
        END IF
      END IF
     !make the value string
      WRITE(Sval,FMT=Sformat) Rval
      !concatenate to produce the output
      Sout=TRIM(ADJUSTL(Sval)) // TRIM(ADJUSTL(Serr))
    ELSE!Error was zero, just output the value without an uncertainty
      WRITE(Sout,FMT='(F8.4)') Rval
    END IF

  END SUBROUTINE UncertBrak

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! A subroutine to call subroutines, not sure why it is needed 
  !!
  !! Major-Authors: Jacob Richardson (2017)
  !!
  SUBROUTINE WriteIterationOutputWrapper(Iter,IThicknessIndex,IPrintFLAG,IErr)

    USE MyNumbers
    USE message_mod
    USE MyMPI
    ! global inputs/outputs
    USE IPARA, ONLY : IPrint, IPreviousPrintedIteration
    USE SPARA, ONLY : SPrintString
    
    IMPLICIT NONE
    INTEGER(IKIND),INTENT(IN) :: Iter,IThicknessIndex,IPrintFLAG
    INTEGER(IKIND),INTENT(OUT) :: IErr

    IF(Iter.EQ.0) IErr=1
    IF(l_alert(IErr,"WriteIterationOutputWrapper","Unexpectedly recieved Iter = 0")) RETURN

    IF(my_rank.EQ.0) THEN
      ! use IPrint and IPrintFLAG to specify how often to write Iteration output
      SELECT CASE(IPrintFLAG)
      
        CASE(0)!We print a simulation only if IPrint<>0 and it's been a while
          IF(IPrint.NE.0.AND.(Iter.GE.(IPreviousPrintedIteration+IPrint))) THEN
            WRITE(SPrintString,FMT='(A16,I2,A35)') ".  Writing output; ",&
              Iter-IPreviousPrintedIteration," iterations since the previous save"
            CALL message (LS,SPrintString)
            CALL WriteIterationOutput(Iter,IThicknessIndex,IErr)
            IPreviousPrintedIteration = Iter 
          END IF
          
        CASE(1)!We print a simulation after a refinement cycle has completed
            CALL message (LS,".  Writing output; end of this refinement cycle")
            CALL WriteIterationOutput(Iter,IThicknessIndex,IErr)
          
        CASE(2)!We print the final simulation
          CALL message ( LS, ".  Writing output; final simulation" )
          CALL WriteIterationOutput(Iter,IThicknessIndex,IErr)
          
      END SELECT
      IF(l_alert(IErr,"WriteIterationOutputWrapper","WriteIterationOutput")) RETURN

    END IF

  END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Writes difference images during gradient determination
  !!
  !! Major-Authors: Richard Beanland (2019)
  !!
  SUBROUTINE WriteDifferenceImages(Iter,IThicknessIndex,IVar,Rx,Rdx,IErr)
  
    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE IPARA, ONLY : ILN,IPixelCount,ISimFLAG,IOutPutReflections,INoOfLacbedPatterns,INhkl,IByteSize
    USE SPARA, ONLY : SChemicalFormula
    USE RPARA, ONLY : Rhkl,RInitialThickness,RDeltaThickness,   RImageSimi
    USE IChannels, ONLY : IChOutWIImage, IChOut
    
    IMPLICIT NONE

    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND), INTENT(IN) :: Iter,IThicknessIndex,IVar
    INTEGER(IKIND) :: IThickness,ind,jnd
    REAL(RKIND) :: Rx,Rdx,Rdelta!The parameter being changed and the amount of change
    REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: RImageToWrite
    CHARACTER(4) :: dString
    CHARACTER(10) :: hString,kString,lString
    CHARACTER(200) :: path,filename,fullpath

    IErr=0
    IThickness = (RInitialThickness + (IThicknessIndex-1)*RDeltaThickness)/10!in nm 

    !Directory for difference image
    IF (IVar.LT.10) THEN
      WRITE(path,"(A1,I4.4,A2,I1,A1,I3.3,A3,I3.3,A1,I3.3)") &
            "I",Iter,"_D",IVar,"_",IThickness,"nm_",2*IPixelcount,"x",2*IPixelcount
    ELSE
      WRITE(path,"(A1,I4.4,A2,I2.2,A1,I3.3,A3,I3.3,A1,I3.3)") &
            "I",Iter,"_D",IVar,"_",IThickness,"nm_",2*IPixelcount,"x",2*IPixelcount
    END IF
    path = SChemicalFormula(1:ILN) // "_" // path ! This adds chemical to folder name
    CALL system('mkdir ' // path)

    ! Write Images to disk
    DO ind = 1,INoOfLacbedPatterns
      ! Make the hkl string e.g. -2-2+10
      jnd=NINT(Rhkl(IOutPutReflections(ind),1))
      IF (ABS(jnd).LT.10) THEN
        WRITE(hString,"(SP,I2.1)") jnd
      ELSEIF (ABS(jnd).LT.100) THEN
        WRITE(hString,"(SP,I3.1)") jnd
      ELSE
        WRITE(hString,"(SP,I4.1)") jnd
      ENDIF
      jnd=NINT(Rhkl(IOutPutReflections(ind),2))
      IF (ABS(jnd).LT.10) THEN
        WRITE(kString,"(SP,I2.1)") jnd
      ELSEIF (ABS(jnd).LT.100) THEN
        WRITE(kString,"(SP,I3.1)") jnd
      ELSE
        WRITE(kString,"(SP,I4.1)") jnd
      ENDIF
      jnd=NINT(Rhkl(IOutPutReflections(ind),3))
      IF (ABS(jnd).LT.10) THEN
        WRITE(lString,"(SP,I2.1)") jnd
      ELSEIF (ABS(jnd).LT.100) THEN
        WRITE(lString,"(SP,I3.1)") jnd
      ELSE
        WRITE(lString,"(SP,I4.1)") jnd
      ENDIF
      ! Make the path/filenames e.g. 'GaAs_-2-2+0.bin'
      filename = SChemicalFormula(1:ILN) // "_" // &
        TRIM(ADJUSTL(hString)) // TRIM(ADJUSTL(kString)) // TRIM(ADJUSTL(lString)) // ".bin"
      fullpath = TRIM(ADJUSTL(path))//"/"//TRIM(ADJUSTL(filename))
      CALL message ( LL, dbg6, fullpath )
      RImageToWrite = RImageSimi(:,:,ind,IThicknessIndex)
      ! Writes data to output image .bin files
      OPEN(UNIT=IChOutWIImage, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(fullpath)),&
          FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=2*IPixelCount*IByteSize)
      IF(l_alert(IErr,"WriteIterationOutput","OPEN() output .bin file")) RETURN      
      DO jnd = 1,2*IPixelCount
        WRITE(IChOutWIImage,rec=jnd) RImageToWrite(jnd,:)
      END DO
      CLOSE(IChOutWIImage,IOSTAT=IErr) 
      IF(l_alert(IErr,"WriteIterationOutput","CLOSE() output .bin file")) RETURN       
    END DO

    !Output of delta values - not necessary since it is already captured in
    !iteration_log.txt
!    Rdelta=Rdx/(Rx-Rdx)!Fractional change in the parameter
!    OPEN(UNIT=IChOut,FILE='differences.txt',FORM='formatted',STATUS='unknown',&
!          POSITION='append')
!    WRITE(UNIT=IChOut,FMT='(I4,1X,F12.9)') Iter,Rdelta
!    CLOSE(IChOut)
    
  END SUBROUTINE WriteDifferenceImages

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Writes output for interations including simulated .bin files, 
  !! structureFactors.txt and structure.cif 
  !!
  !! Major-Authors: 'kidwhizz' (2015), Richard Beanland (2016)
  !!
  SUBROUTINE WriteIterationOutput(Iter,IThicknessIndex,IErr)

    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE IPARA, ONLY : ILN,IPixelCount,ISimFLAG,IOutPutReflections,INoOfLacbedPatterns,INhkl,IByteSize
    USE CPARA, ONLY : CUgMat
    USE RPARA, ONLY : Rhkl,RgPool, RImageSimi, RInitialThickness, RDeltaThickness
    USE SPARA, ONLY : SChemicalFormula
    USE IChannels, ONLY : IChOutWIImage, IChOut

    IMPLICIT NONE

    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND), INTENT(IN) :: Iter,IThicknessIndex
    INTEGER(IKIND) :: IThickness,ind,jnd
    REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: RImageToWrite
    CHARACTER(10) :: hString,kString,lString
    CHARACTER(200) :: path,filename,fullpath
    
    IErr=0 

    IF (ISimFLAG.EQ.0) THEN !felixrefine output
      IThickness = NINT((RInitialThickness +(IThicknessIndex-1)*RDeltaThickness)/TEN)!in nm 
      WRITE(path,"(A1,I4.4,A1,I3.3,A3,I3.3,A1,I3.3)") &
            "I",Iter,"_",IThickness,"nm_",2*IPixelcount,"x",2*IPixelcount
    ELSE ! Sim Output
    IThickness = NINT(RInitialThickness +(IThicknessIndex-1)*RDeltaThickness)!in A 
      WRITE(path,"(A4,I4.4,A2,I3.3,A1,I3.3)") &
            "Sim_",IThickness,"A_",2*IPixelcount,"x",2*IPixelcount
    END IF
    path = SChemicalFormula(1:ILN) // "_" // path ! This adds chemical to folder name
    CALL system('mkdir ' // path)

    ! Write Images to disk
    DO ind = 1,INoOfLacbedPatterns
      ! Make the hkl string e.g. -2-2+10
      jnd=NINT(Rhkl(IOutPutReflections(ind),1))
      IF (ABS(jnd).LT.10) THEN
        WRITE(hString,"(SP,I2.1)") jnd
      ELSEIF (ABS(jnd).LT.100) THEN
        WRITE(hString,"(SP,I3.1)") jnd
      ELSE
        WRITE(hString,"(SP,I4.1)") jnd
      ENDIF
      jnd=NINT(Rhkl(IOutPutReflections(ind),2))
      IF (ABS(jnd).LT.10) THEN
        WRITE(kString,"(SP,I2.1)") jnd
      ELSEIF (ABS(jnd).LT.100) THEN
        WRITE(kString,"(SP,I3.1)") jnd
      ELSE
        WRITE(kString,"(SP,I4.1)") jnd
      ENDIF
      jnd=NINT(Rhkl(IOutPutReflections(ind),3))
      IF (ABS(jnd).LT.10) THEN
        WRITE(lString,"(SP,I2.1)") jnd
      ELSEIF (ABS(jnd).LT.100) THEN
        WRITE(lString,"(SP,I3.1)") jnd
      ELSE
        WRITE(lString,"(SP,I4.1)") jnd
      ENDIF
      ! Make the path/filenames e.g. 'GaAs_-2-2+0.bin'
      filename = SChemicalFormula(1:ILN) // "_" // TRIM(ADJUSTL(hString)) // TRIM(ADJUSTL(kString)) // TRIM(ADJUSTL(lString)) // ".bin"
      fullpath = TRIM(ADJUSTL(path))//"/"//TRIM(ADJUSTL(filename))
      CALL message ( LL, dbg6, fullpath )
      RImageToWrite = RImageSimi(:,:,ind,IThicknessIndex)
!DBG    IF(my_rank.EQ.0)PRINT*,TRIM(ADJUSTL(fullpath))!DEBUG
      ! Writes data to output image .bin files
      OPEN(UNIT=IChOutWIImage, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(fullpath)),&
          FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=2*IPixelCount*IByteSize)
      IF(l_alert(IErr,"WriteIterationOutput","OPEN() output .bin file")) RETURN      
      DO jnd = 1,2*IPixelCount
        WRITE(IChOutWIImage,rec=jnd) RImageToWrite(jnd,:)
      END DO
      CLOSE(IChOutWIImage,IOSTAT=IErr) 
      IF(l_alert(IErr,"WriteIterationOutput","CLOSE() output .bin file")) RETURN       
    END DO

    ! writes out structure.cif
    CALL WriteIterationCIF(Iter,path,IErr)
    IF(l_alert(IErr,"WriteIterationOutput","WriteIterationCIF")) RETURN   
    IErr=0
    ! write out StructureFactors.txt
    WRITE(filename,*) "StructureFactors.txt"
    WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))
    OPEN(UNIT=IChOut,STATUS='UNKNOWN',FILE=TRIM(ADJUSTL(fullpath)))

    DO ind = 1,INhkl
      WRITE(IChOut,FMT='(3I5.1,2(1X,F13.9),2(1X,E14.6))') NINT(Rhkl(ind,:)),&
              RgPool(ind,1),RgPool(ind,2),CUgMat(ind,1)
    END DO

    CLOSE(IChOut)    
    
    RETURN  
    
  END SUBROUTINE WriteIterationOutput

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Write out structure.cif containing non symmetrically relate
  !! atomic positions.
  !!
  !! Major-Authors: 'kidwhizz' (2015), Richard Beanland (2016)
  !!
  SUBROUTINE WriteIterationCIF(Iter,path,IErr)

    USE MyNumbers
    USE message_mod
    USE setup_space_group_mod!required for the subroutine ConvertSpaceGroupToNumber
 
    ! global inputs
    USE IPARA, ONLY : ISpaceGrp,ILN,IIndependentVariableAtom,IIndependentVariableType
    USE RPARA, ONLY : RLengthX, RLengthY, RLengthZ, RAlpha, RBeta, RGamma, &
          RBasisAtomPosition, RBasisAtomDelta, RIndependentDelta, RBasisIsoDW, RBasisOccupancy, RVolume
    USE SPARA, ONLY : SSpaceGrp, SBasisAtomLabel, SBasisAtomName,SChemicalFormula,&
          SSymString
    USE IChannels, ONLY : IChOutSimplex

    IMPLICIT NONE

    CHARACTER(200), INTENT(IN) :: path
    INTEGER(IKIND), INTENT(IN) :: Iter
    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND) :: ind,jnd
    CHARACTER(200) :: filename, fullpath
    CHARACTER(100) :: String
    CHARACTER(14) :: Sout

    IErr=0
    ! Write out unique atomic positions
    WRITE(filename,"(A1,I4.4,A4)") "_",Iter,".cif"
    filename=SChemicalFormula(1:ILN) // filename!gives e.g. SrTiO3_0001.cif 
    WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))

    OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',FILE=TRIM(ADJUSTL(fullpath)))
 
    !Introductory lines
    WRITE(IChOutSimplex,FMT='(A31)') "#(C) 2019 University of Warwick"
    WRITE(IChOutSimplex,FMT='(A16)') "data_felixrefine"
    WRITE(IChOutSimplex,FMT='(A31)') "_audit_creation_date 2019-08-06"!need to find out how to get a date!
    WRITE(IChOutSimplex,FMT='(A22,A)') "_chemical_formula_sum ",TRIM(ADJUSTL(SChemicalFormula))

! Citation data would go here

    !unit cell
    WRITE(IChOutSimplex,FMT='(A5)') "loop_"
    WRITE(IChOutSimplex,FMT='(A14,1X,F7.4)') "_cell_length_a",RLengthX
    WRITE(IChOutSimplex,FMT='(A14,1X,F7.4)') "_cell_length_b",RLengthY
    WRITE(IChOutSimplex,FMT='(A14,1X,F7.4)') "_cell_length_c",RLengthZ
    WRITE(IChOutSimplex,FMT='(A17,1X,F7.2)') "_cell_angle_alpha",RAlpha*180/PI
    WRITE(IChOutSimplex,FMT='(A17,1X,F7.2)') "_cell_angle_beta ",RBeta*180/PI
    WRITE(IChOutSimplex,FMT='(A17,1X,F7.2)') "_cell_angle_gamma",RGamma*180/PI
    WRITE(IChOutSimplex,FMT='(A12,1X,F7.2)') "_cell_volume",RVolume
    !symmetry
    WRITE(IChOutSimplex,FMT='(A,A,A)') "_symmetry_space_group_name_H-M  '",TRIM(ADJUSTL(SSpaceGrp)),"'"
    WRITE(IChOutSimplex,FMT='(A27,1X,I3)') "_symmetry_Int_Tables_number",ISpaceGrp
    WRITE(IChOutSimplex,FMT='(A5)') "loop_"
    WRITE(IChOutSimplex,FMT='(A27)') "_symmetry_equiv_pos_site_id"
    WRITE(IChOutSimplex,FMT='(A26)') "_symmetry_equiv_pos_as_xyz"
    DO jnd = 1,SIZE(SSymString,DIM=1)!
      WRITE(IChOutSimplex,FMT='(I3,1X,A30)') jnd,SSymString(jnd)
    END DO
    !atom coordinates etc
    WRITE(IChOutSimplex,FMT='(A5)') "loop_"
    WRITE(IChOutSimplex,FMT='(A16)') "_atom_site_label"
    WRITE(IChOutSimplex,FMT='(A22)') "_atom_site_type_symbol"
!    WRITE(IChOutSimplex,FMT='(A25)') "_atom_site_symmetry_multiplicity"
!    WRITE(IChOutSimplex,FMT='(A25)') "_atom_site_Wyckoff_symbol"
    WRITE(IChOutSimplex,FMT='(A18)') "_atom_site_fract_x"
    WRITE(IChOutSimplex,FMT='(A18)') "_atom_site_fract_y"
    WRITE(IChOutSimplex,FMT='(A18)') "_atom_site_fract_z"
    WRITE(IChOutSimplex,FMT='(A25)') "_atom_site_B_iso_or_equiv"
    WRITE(IChOutSimplex,FMT='(A20)') "_atom_site_occupancy"

    !Make output string by appending each part
    DO jnd = 1,SIZE(RBasisAtomPosition,DIM=1)
      !Label and name
      WRITE(String,FMT='(A3,1X,A3,1X)')SBasisAtomLabel(jnd), SBasisAtomName(jnd)
      !Atom coords using the uncertainties in RBasisAtomDelta
      DO ind = 1,3
        CALL UncertBrak(RBasisAtomPosition(jnd,ind),RBasisAtomDelta(jnd,ind),Sout,IErr)
        String = TRIM(ADJUSTL(String)) // "  " // TRIM(ADJUSTL(Sout))!append onto String
      END DO
      !Isotropic DWF
      WRITE(Sout,FMT='(F7.4)') RBasisIsoDW(jnd)
      !replace Sout if it is being refined
      DO ind = 1,SIZE(IIndependentVariableAtom)
        IF(IIndependentVariableType(ind).EQ.4.AND.jnd.EQ.IIndependentVariableAtom(ind))THEN
          CALL UncertBrak(RBasisIsoDW(jnd),RIndependentDelta(ind),Sout,IErr)
        END IF
      END DO
      String = TRIM(ADJUSTL(String)) // "  " // TRIM(ADJUSTL(Sout))!append onto String
      !Occupancy
      WRITE(Sout,FMT='(F7.4)') RBasisOccupancy(jnd)
      DO ind = 1,SIZE(IIndependentVariableAtom)
        IF(IIndependentVariableType(ind).EQ.3.AND.jnd.EQ.IIndependentVariableAtom(ind))THEN
          CALL UncertBrak(RBasisOccupancy(jnd),RIndependentDelta(ind),Sout,IErr)
        END IF
      END DO
       String = TRIM(ADJUSTL(String)) // "  " // TRIM(ADJUSTL(Sout))
      WRITE(IChOutSimplex,FMT='(A)') String
    END DO
    WRITE(IChOutSimplex,FMT='(A)') "#End of felixrefine cif"
    
    CLOSE(IChOutSimplex)

!    Write out full atomic positions

!    WRITE(filename,*) "StructureFull.txt"
!    WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))
!    CALL message( LL, "RAtomPosition,SAtomName" )
!    OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',FILE=TRIM(ADJUSTL(fullpath)))
!    DO jnd = 1,SIZE(RAtomPosition,DIM=1)
!      WRITE(IChOutSimplex,FMT='(A2,1X,3(F9.6,1X))') SAtomName(jnd),RAtomPosition(jnd,1:3)
!    END DO
!    CLOSE(IChOutSimplex)

  END SUBROUTINE WriteIterationCIF

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Adds the current fit and simulation parameters to IterationLog.txt
  !!
  !! Major-Authors: 'kidwhizz' (2015), Richard Beanland (2016)
  !!
  SUBROUTINE WriteOutVariables(Iter,IErr)

    USE MyNumbers
    USE message_mod

    ! global inputs
    USE IPARA, ONLY : IAbsorbFLAG, IRefineMode, INoofUgs, IUgOffset, &
                      IRefinementVariableTypes
    USE RPARA, ONLY : RBasisAtomPosition, RBasisOccupancy, RBasisIsoDW, &
                      RAnisotropicDebyeWallerFactorTensor, RFigureofMerit, &
                      RAbsorptionPercentage, RLengthX, RLengthY, RLengthZ, RAlpha, RBeta, &
                      RGamma, RConvergenceAngle, RAcceleratingVoltage                    
    USE CPARA, ONLY : CUniqueUg
    USE IChannels, ONLY : IChOutSimplex     

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(IN) :: Iter
    INTEGER(IKIND) :: IErr,ind,IStart,IEnd,jnd,ITotalOutputVariables
    CHARACTER(200) :: SFormat, STotalOutputVariables
    INTEGER(IKIND),DIMENSION(IRefinementVariableTypes) :: IOutputVariables
    REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RDataOut

    ! Need to Determine total no. of variables to be written out
    ! this is different from the no. of refinement variables
    
    IF (IAbsorbFLAG.EQ.2) THEN ! Structure Factors are complex 
      ! so require two output variables each
      IOutputVariables(1) = IRefineMode(1)*2*INoofUgs+1 
      ! plus one for proportional absorption     
    ELSE
      IOutputVariables(1) = IRefineMode(1)*2*INoofUgs     
    END IF
    ! Atom Coordinates
    IOutputVariables(2) = IRefineMode(2)*SIZE(RBasisAtomPosition,DIM=1)* &
          SIZE(RBasisAtomPosition,DIM=2)
    ! Occupancies
    IOutputVariables(3) = IRefineMode(3)*SIZE(RBasisOccupancy,DIM=1) 
    ! Isotropic Debye Waller Factors
    IOutputVariables(4) = IRefineMode(4)*SIZE(RBasisIsoDW,DIM=1) 
    ! Anisotropic Debye Waller Factors
    IOutputVariables(5) = IRefineMode(5)*SIZE(RAnisotropicDebyeWallerFactorTensor)
    IOutputVariables(6) = IRefineMode(6) * 3 ! Lattice Parameters (a,b,c) 
    IOutputVariables(7) = IRefineMode(7) * 3 ! Lattice Angles (alpha,beta,gamma)
    IOutputVariables(8) = IRefineMode(8) ! Convergence angle
    IOutputVariables(9) = IRefineMode(9) ! Absorption
    IOutputVariables(10) = IRefineMode(10) ! Accelerating Voltage
    ITotalOutputVariables = SUM(IOutputVariables) ! Total Output

    ALLOCATE(RDataOut(ITotalOutputVariables),STAT=IErr)
    DO jnd = 1,IRefinementVariableTypes
      IF(IRefineMode(jnd).EQ.0) THEN
        CYCLE ! The refinement variable type is not being refined, skip
      END IF
      IF(jnd.EQ.1) THEN ! It's an atom coordinate refinement
        IStart = 1
      ELSE
        IStart = SUM(IOutputVariables(1:(jnd-1)))+1 
        !?? RB there is probably a better way of doing this
      END IF
      IEND = SUM(IOutputVariables(1:jnd))

      SELECT CASE(jnd)
      CASE(1)!A, Ug's
        DO ind = 1,INoofUgs
           IStart = (ind*2)-1
           IEnd = ind*2
           RDataOut(IStart:IEnd) = [REAL(CUniqueUg(ind+IUgOffset)), &
                  REAL(AIMAG(CUniqueUg(ind+IUgOffset)),RKIND)]
        END DO
        RDataOut(IEnd+1) = RAbsorptionPercentage!RB last variable is absorption
      CASE(2)!B, Atom coords
        RDataOut(IStart:IEnd) = &
              RESHAPE(TRANSPOSE(RBasisAtomPosition),SHAPE(RDataOut(IStart:IEnd)))
      CASE(3)!C, Occupancy
        RDataOut(IStart:IEnd) = RBasisOccupancy
      CASE(4)!D, Isotropic DWFs
        RDataOut(IStart:IEnd) = RBasisIsoDW
      CASE(5)!E, Anisotropic DWFs
        RDataOut(IStart:IEnd) = &
              RESHAPE(RAnisotropicDebyeWallerFactorTensor,SHAPE(RDataOut(IStart:IEnd)))
      CASE(6)!F, Lattice parameter
        RDataOut(IStart:IEnd) = [RLengthX, RLengthY, RLengthZ]
      CASE(7)!G, Unit cell angles
        RDataOut(IStart:IEnd) = [RAlpha, RBeta, RGamma]
      CASE(8)!H, convergence angle
        RDataOut(IStart:IEnd) = RConvergenceAngle
      CASE(9)!I, kV
        RDataOut(IStart:IEnd) = RAcceleratingVoltage
      END SELECT
    END DO

    WRITE(STotalOutputVariables,*) ITotalOutputVariables
    WRITE(SFormat,*) "(I5.1,1X,F13.9,1X,"//TRIM(ADJUSTL(STotalOutputVariables))//"(F13.9,1X))"

    OPEN(UNIT=IChOutSimplex,FILE='iteration_log.txt',FORM='formatted',STATUS='unknown',&
          POSITION='append')
    WRITE(UNIT=IChOutSimplex,FMT=SFormat) Iter,RFigureofMerit,RDataOut
    CLOSE(IChOutSimplex)

  END SUBROUTINE WriteOutVariables

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: This can be used as a quick way to visually compare .img input 
  !! images and simulated .bin images. This scales each experimental image such that the max value
  !! matches the corresponding final simulated image and writes it out in an equaivalent
  !! .bin format in the current directory.
  !!
  !! Major-Authors: Jacob Richardson (2017)
  !!
  SUBROUTINE NormaliseExperimentalImagesAndWriteOut(IThicknessIndex,IErr)

    USE MyNumbers
    USE message_mod

    !global inputs
    USE RPara, ONLY : RImageExpi, RImageSimi, Rhkl, RInitialThickness, RDeltaThickness
    USE IPara, ONLY : IPixelcount, IOutPutReflections, INoOfLacbedPatterns,ILN
    USE Spara, ONLY : SChemicalFormula
    USE IChannels, ONLY : IChOutWIImage, IChOut

    IMPLICIT NONE   
    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND), INTENT(IN) :: IThicknessIndex
    CHARACTER(200) :: path, filename, fullpath
    INTEGER(IKIND) :: IThickness, ind, jnd
    REAL(RKIND), ALLOCATABLE :: RImageToWrite(:,:)

    CALL message(LS,"Writing out normalised experimental images")

    IThickness = (RInitialThickness + (IThicknessIndex-1)*RDeltaThickness)/10!in nm 

    WRITE(path,"(A,A6,I3.3,A1,I3.3)") &
           SChemicalFormula(1:ILN),"_expt_",2*IPixelcount,"x",2*IPixelcount

    CALL system('mkdir ' // path)

    ! Write Images to disk
    DO ind = 1,INoOfLacbedPatterns
      ! Make the path/filenames  
      WRITE(filename,"(A,A,I3.3,A1,I3.3,A1,SP,3(I2.1),A)")&
             SChemicalFormula(1:ILN),"_experi_",2*IPixelcount,"x",2*IPixelcount,&
            "_",NINT(Rhkl(IOutPutReflections(ind),1:3)),'.bin'
      fullpath=TRIM(path)//'/'//filename

      RImageToWrite = RImageExpi(:,:,ind)
      RImageToWrite = RImageToWrite/ MAXVAL(RImageToWrite)*MAXVAL(RImageSimi(:,:,ind,IThicknessIndex))

      ! Writes data to output image .bin files
      OPEN(UNIT=IChOutWIImage, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(fullpath)),&
            FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=2*IPixelCount*8)
      IF(l_alert(IErr,"WriteIterationOutput","OPEN() output .bin file")) RETURN       
      DO jnd = 1,2*IPixelCount
        WRITE(IChOutWIImage,rec=jnd) RImageToWrite(jnd,:)
      END DO
      CLOSE(IChOutWIImage,IOSTAT=IErr) 
      IF(l_alert(IErr,"WriteIterationOutput","CLOSE() output .bin file")) RETURN       
    END DO   

  END SUBROUTINE NormaliseExperimentalImagesAndWriteOut


END MODULE write_output_mod
