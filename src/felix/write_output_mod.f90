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
!! Module-description: 
!!
MODULE write_output_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: WriteIterationOutput, WriteIterationCIF, WriteOutVariables

  CONTAINS

  !>
  !! Procedure-description: 
  !!
  !! Major-Authors: 'kidwhizz' (2015), Richard Beanland (2016)
  !!
  SUBROUTINE WriteIterationOutput(Iter,IThicknessIndex,IExitFlag,IErr)

    !?? called in felixrefine & SimulateAndFit

    USE MyNumbers
    USE message_mod
    
    ! global outputs
    USE IPARA, ONLY : IPreviousPrintedIteration
    
    ! global inputs
    USE IPARA, ONLY : IPixelCount, ISimFLAG, IOutPutReflections, INoOfLacbedPatterns, &
                      nReflections
    USE CPARA, ONLY : CUgMat
    USE RPARA, ONLY : Rhkl, RImageSimi, RInitialThickness, RDeltaThickness
    USE SPARA, ONLY : SChemicalFormula
    USE IChannels, ONLY : IChOutWIImage, IChOut

    IMPLICIT NONE

    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND), INTENT(IN) :: Iter, IThicknessIndex, IExitFLAG
    INTEGER(IKIND) :: IThickness,ind,jnd
    REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: RImageToWrite
    CHARACTER*200 :: path, filename, fullpath
    CHARACTER*20 :: SIntString
    
    IThickness = (RInitialThickness + (IThicknessIndex-1)*RDeltaThickness)/10!in nm 

    IF (ISimFLAG.EQ.0) THEN !felixrefine output
      WRITE(path,"(A1,I4.4,A1,I3.3,A3,I3.3,A1,I3.3)") &
            "I",Iter,"_",IThickness,"nm_",2*IPixelcount,"x",2*IPixelcount
      path = TRIM(SChemicalFormula) // "_" // path ! This adds chemical to folder name
    ELSE ! Sim Output
      WRITE(path,"(A4,I3.3,A3,I3.3,A1,I3.3)") &
            "Sim_",IThickness,"nm_",2*IPixelcount,"x",2*IPixelcount
    END IF
    CALL system('mkdir ' // path)

    IF (ISimFLAG.EQ.0.AND.IExitFLAG.EQ.0) THEN ! felixrefine output
      IF (IPreviousPrintedIteration.EQ.0) THEN
        CALL message ( LS, "Writing output; baseline simulation" )
      ELSE
        CALL message ( LS, "Writing output; iterations since the previous save = ", &
              Iter-IPreviousPrintedIteration)
      END IF
    END IF
    
    ! Write Images to disk
    DO ind = 1,INoOfLacbedPatterns
      ! Make the path/filenames  
      ! Iterates over 3 vector components to make filename e.g. 'GaAs-2-2+0.bin.
      WRITE(filename,"(A1,I3.3,A3,I3.3,A1,I3.3,A1)")"_",IThickness,"nm_",&
                  2*IPixelcount,"x",2*IPixelcount,"_"
      filename = TRIM(ADJUSTL(path))//"/"//TRIM(ADJUSTL(SChemicalFormula))&
                  //TRIM(ADJUSTL(filename))    
      DO jnd = 1,3
        WRITE(SIntString,*) NINT(Rhkl(IOutPutReflections(ind),jnd))
        IF (NINT(Rhkl(IOutPutReflections(ind),jnd)) >= 0) THEN    
          filename = TRIM(filename) // '+'
        END IF
        filename = TRIM(filename) // TRIM(ADJUSTL(SIntString))
      END DO
      filename = TRIM(filename) // '.bin'

      CALL message ( LL, dbg6, filename )
      RImageToWrite = RImageSimi(:,:,ind,IThicknessIndex)	

      ! Writes data to output image .bin files
      OPEN(UNIT=IChOutWIImage, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),&
	          FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=2*IPixelCount*8)
      IF(l_alert(IErr,"WriteIterationOutput","OPEN() output .bin file")) RETURN       
      DO jnd = 1,2*IPixelCount
        WRITE(IChOutWIImage,rec=jnd) RImageToWrite(jnd,:)
      END DO
      CLOSE(IChOutWIImage,IOSTAT=IErr) 
      IF(l_alert(IErr,"WriteIterationOutput","CLOSE() output .bin file")) RETURN       
    END DO

    CALL WriteIterationCIF(path,IErr) 
    IF(l_alert(IErr,"WriteIterationOutput","WriteIterationCIF")) RETURN   

    WRITE(filename,*) "StructureFactors.txt"
    WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))
    OPEN(UNIT=IChOut,STATUS='UNKNOWN',FILE=TRIM(ADJUSTL(fullpath)))

    DO ind = 1,nReflections
       WRITE(IChOut,FMT='(3I5.1,2F13.9)') NINT(Rhkl(ind,:)),CUgMat(ind,1)
    END DO

    CLOSE(IChOut)    

    IPreviousPrintedIteration = Iter 
    
    RETURN  
    
  END SUBROUTINE WriteIterationOutput

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: This scales each experimental image such that the max value
  !! matches the corresponding final simulated image and writes it out in an equaivalent
  !! .bin format.
  !!
  !! Major-Authors: Jacob Richardson (2017)
  !!
  SUBROUTINE NormaliseExperimentalImagesAndWriteOut(IThicknessIndex,IErr)

    USE MyNumbers
    USE message_mod

    !global inputs
    USE RPara, ONLY : RImageExpi, RImageSimi, Rhkl, RInitialThickness, RDeltaThickness
    USE IPara, ONLY : IPixelcount, IOutPutReflections, INoOfLacbedPatterns
    USE Spara, ONLY : SChemicalFormula
    USE IChannels, ONLY : IChOutWIImage, IChOut

    IMPLICIT NONE   
    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND), INTENT(IN) :: IThicknessIndex
    CHARACTER(200) :: path, filename
    CHARACTER(20) :: SIntString
    INTEGER(IKIND) :: IThickness, ind, jnd
    REAL(RKIND), ALLOCATABLE :: RImageToWrite(:,:)

    IThickness = (RInitialThickness + (IThicknessIndex-1)*RDeltaThickness)/10!in nm 

    CALL message('printing experimental images')
    WRITE(path,"(A6,I3.3,A3,I3.3,A1,I3.3)") &
          "approx",IThickness,"nm_",2*IPixelcount,"x",2*IPixelcount
    path = "experi_images_"// TRIM(SChemicalFormula) // "_" // path

    CALL system('mkdir ' // path)
    
    ! Write Images to disk
    DO ind = 1,INoOfLacbedPatterns
      ! Make the path/filenames  
      ! Iterates over 3 vector components to make filename e.g. 'GaAs-2-2+0.bin.
      WRITE(filename,"(A1,I3.3,A3,I3.3,A1,I3.3,A1)")"_approx",IThickness,"nm_",&
                  2*IPixelcount,"x",2*IPixelcount,"_"
      filename = TRIM(ADJUSTL(path))//"/"//TRIM(ADJUSTL(SChemicalFormula))&
                  //TRIM(ADJUSTL(filename))
      DO jnd = 1,3
        WRITE(SIntString,*) NINT(Rhkl(IOutPutReflections(ind),jnd))
        IF (NINT(Rhkl(IOutPutReflections(ind),jnd)) >= 0) THEN    
          filename = TRIM(filename) // '+'
        END IF
        filename = TRIM(filename) // TRIM(ADJUSTL(SIntString))
      END DO
      filename = TRIM(filename) // '.bin'

      RImageToWrite = RImageExpi(:,:,ind)
      RImageToWrite = RImageToWrite / MAXVAL(RImageToWrite) * MAXVAL(RImageSimi(:,:,ind,IThicknessIndex))

      ! Writes data to output image .bin files
      OPEN(UNIT=IChOutWIImage, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),&
            FORM='UNFORMATTED',ACCESS='DIRECT',IOSTAT=IErr,RECL=2*IPixelCount*8)
      IF(l_alert(IErr,"WriteIterationOutput","OPEN() output .bin file")) RETURN       
      DO jnd = 1,2*IPixelCount
        WRITE(IChOutWIImage,rec=jnd) RImageToWrite(jnd,:)
      END DO
      CLOSE(IChOutWIImage,IOSTAT=IErr) 
      IF(l_alert(IErr,"WriteIterationOutput","CLOSE() output .bin file")) RETURN       
    END DO   

  END SUBROUTINE NormaliseExperimentalImagesAndWriteOut

  !>
  !! Procedure-description: Write out non symmetrically related atomic positions
  !!
  !! Major-Authors: 'kidwhizz' (2015), Richard Beanland (2016)
  !!
  SUBROUTINE WriteIterationCIF(path,IErr)

    !?? called in felixrefine & SimulateAndFit, via WriteIterationOutput() 

    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE RPARA, ONLY : RLengthX, RLengthY, RLengthZ, RAlpha, RBeta, RGamma, &
                      RBasisAtomPosition, RBasisIsoDW, RBasisOccupancy
    USE SPARA, ONLY : SSpaceGrp, SBasisAtomLabel, SBasisAtomName
    USE IChannels, ONLY : IChOutSimplex

    IMPLICIT NONE

    CHARACTER*200, INTENT(IN) :: path
    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND) :: jnd
    CHARACTER*200 :: filename, fullpath

    ! Write out non symmetrically related atomic positions

    WRITE(filename,*) "structure.cif"
    WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))

    OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',FILE=TRIM(ADJUSTL(fullpath)))
    ! RB
    WRITE(IChOutSimplex,FMT='(A16)') "data_felixrefine"
    WRITE(IChOutSimplex,FMT='(A5)') "loop_"
    WRITE(IChOutSimplex,FMT='(A14,1X,F7.4)') "_cell_length_a",RLengthX
    WRITE(IChOutSimplex,FMT='(A14,1X,F7.4)') "_cell_length_b",RLengthY
    WRITE(IChOutSimplex,FMT='(A14,1X,F7.4)') "_cell_length_c",RLengthZ
    WRITE(IChOutSimplex,FMT='(A17,1X,F7.2)') "_cell_angle_alpha",RAlpha*180/PI
    WRITE(IChOutSimplex,FMT='(A16,1X,F7.2)') "_cell_angle_beta",RBeta*180/PI
    WRITE(IChOutSimplex,FMT='(A17,1X,F7.2)') "_cell_angle_gamma",RGamma*180/PI
    WRITE(IChOutSimplex,FMT='(A32,A10,A1)') "_symmetry_space_group_name_H-M '",SSpaceGrp,"'"
    WRITE(IChOutSimplex,FMT='(A5)') " "
    WRITE(IChOutSimplex,FMT='(A5)') "loop_"
    WRITE(IChOutSimplex,FMT='(A22)') "_atom_site_label"
    WRITE(IChOutSimplex,FMT='(A22)') "_atom_site_type_symbol"
!    WRITE(IChOutSimplex,FMT='(A25)') "_atom_site_symmetry_multiplicity"
!    WRITE(IChOutSimplex,FMT='(A25)') "_atom_site_Wyckoff_symbol"
    WRITE(IChOutSimplex,FMT='(A18)') "_atom_site_fract_x"
    WRITE(IChOutSimplex,FMT='(A18)') "_atom_site_fract_y"
    WRITE(IChOutSimplex,FMT='(A18)') "_atom_site_fract_z"
    WRITE(IChOutSimplex,FMT='(A25)') "_atom_site_B_iso_or_equiv"
    WRITE(IChOutSimplex,FMT='(A20)') "_atom_site_occupancy"

    DO jnd = 1,SIZE(RBasisAtomPosition,DIM=1)!RB only gives refined atoms, needs work
      WRITE(IChOutSimplex,FMT='(2(A3,1X),3(F7.4,1X),2(F5.2,1X))') &
	          SBasisAtomLabel(jnd), SBasisAtomName(jnd), RBasisAtomPosition(jnd,:), &
            RBasisIsoDW(jnd),RBasisOccupancy(jnd)
    END DO
    WRITE(IChOutSimplex,FMT='(A22)') "#End of refinement cif"
    
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
                      RGamma, RConvergenceAngle, RAcceleratingVoltage, RRSoSScalingFactor                    
    USE CPARA, ONLY : CUniqueUg
    USE IChannels, ONLY : IChOutSimplex     

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(IN) :: Iter
    INTEGER(IKIND) :: IErr,ind,IStart,IEnd,jnd,ITotalOutputVariables
    CHARACTER*200 :: SFormat, STotalOutputVariables
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
      CASE(1)
        DO ind = 1,INoofUgs
           IStart = (ind*2)-1
           IEnd = ind*2
           RDataOut(IStart:IEnd) = [REAL(CUniqueUg(ind+IUgOffset)), &
                  REAL(AIMAG(CUniqueUg(ind+IUgOffset)),RKIND)]
        END DO
        RDataOut(IEnd+1) = RAbsorptionPercentage!RB last variable is absorption
      CASE(2)
        RDataOut(IStart:IEnd) = &
              RESHAPE(TRANSPOSE(RBasisAtomPosition),SHAPE(RDataOut(IStart:IEnd)))
      CASE(3)
        RDataOut(IStart:IEnd) = RBasisOccupancy
      CASE(4)
        RDataOut(IStart:IEnd) = RBasisIsoDW
      CASE(5)
        RDataOut(IStart:IEnd) = &
              RESHAPE(RAnisotropicDebyeWallerFactorTensor,SHAPE(RDataOut(IStart:IEnd)))
      CASE(6)
        RDataOut(IStart:IEnd) = [RLengthX, RLengthY, RLengthZ]
      CASE(7)
        RDataOut(IStart:IEnd) = [RAlpha, RBeta, RGamma]
      CASE(8)
        RDataOut(IStart:IEnd) = RConvergenceAngle
      CASE(9)
        RDataOut(IStart:IEnd) = RAbsorptionPercentage
      CASE(10)
        RDataOut(IStart:IEnd) = RAcceleratingVoltage
      CASE(11)
        RDataOut(IStart:IEnd) = RRSoSScalingFactor
      CASE(12)
          DO ind = 1,INoofUgs
             IStart = (ind*2)-1
             IEnd = ind*2
             RDataOut(IStart:IEnd) = &
                  [REAL(CUniqueUg(ind+IUgOffset)), REAL(AIMAG(CUniqueUg(ind+IUgOffset)),RKIND)]
          END DO
          IF (IAbsorbFLAG.EQ.1) THEN
		        RDataOut(IEnd+1) = RAbsorptionPercentage 
            ! RB last variable is proportional absorption
          END IF
      END SELECT
    END DO

    WRITE(STotalOutputVariables,*) ITotalOutputVariables
    WRITE(SFormat,*) "(I5.1,1X,F13.9,1X,"//TRIM(ADJUSTL(STotalOutputVariables))//"(F13.9,1X))"

    OPEN(UNIT=IChOutSimplex,file='iteration_log.txt',form='formatted',status='unknown',&
          position='append')
    WRITE(UNIT=IChOutSimplex,FMT=SFormat) Iter-1,RFigureofMerit,RDataOut
    CLOSE(IChOutSimplex)

  END SUBROUTINE WriteOutVariables

END MODULE write_output_mod
