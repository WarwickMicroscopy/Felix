!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixrefine
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all right reserved
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
!  This file is part of felixrefine.
!
!  felixrefine is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  felixrefine is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with felixrefine.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: Felixrefine.f90,v 1.89 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE WriteIterationOutput(Iter,IThicknessIndex,IExitFlag,IErr)
!This code needs to be taken up a subroutine level and combined with CalculateFigureofMeritandDetermineThickness to avoid repeated calculation of the image to output and BlurG
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,Iter,IThickness,ind,jnd,gnd,hnd
  INTEGER(IKIND),INTENT(IN) :: IThicknessIndex,IExitFLAG
  REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: RImageToWrite
  REAL(RKIND) :: Rradius
  CHARACTER*200 :: path,SPrintString
!PRINT*,"Iter",Iter,IPreviousPrintedIteration,IPrint
  IF(IExitFLAG.EQ.1.OR.(Iter.GE.(IPreviousPrintedIteration+IPrint))) THEN
     IThickness = (RInitialThickness + (IThicknessIndex-1)*RDeltaThickness)/10!RB in nm 
     WRITE(path,"(A10,I4.4,A1,I3.3,A3,I3.3,A1,I3.3)") &
          "Iteration",Iter,&
          "_",IThickness,&
          "nm_",2*IPixelcount,&
          "x",2*IPixelcount
     
     call system('mkdir ' // path)

     IF (IExitFLAG.EQ.0) THEN
	   IF (IPreviousPrintedIteration.LT.0) THEN
         WRITE(SPrintString,FMT='(A35)') "Writing output; baseline simulation"
	   ELSE
         WRITE(SPrintString,FMT='(A16,I4,A35)') "Writing output; ",&
		 Iter-IPreviousPrintedIteration," iterations since the previous save"
	   END IF
     ELSE
        WRITE(SPrintString,FMT='(A28,I4,A35)') "Exiting and writing output; ",&
		Iter-IPreviousPrintedIteration," iterations since the previous save" 
	 END IF
        PRINT*,TRIM(ADJUSTL(SPrintString))
     
     IPreviousPrintedIteration = Iter

!    Was WriteIterationImages
  DO ind = 1,INoOfLacbedPatterns
     CALL OpenReflectionImage(IChOutWIImage,path,IErr,ind,2*IPixelCount,2_IKIND)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationImages(",my_rank,") error in OpenReflectionImage()"
        RETURN
     ENDIF
	 !convert vector RSimulatedPatterns into 2D RImageToWrite (again, was done once in CalculateFigureofMeritandDetermineThickness
     RImageToWrite = ZERO
     DO jnd = 1,IPixelTotal
        gnd = IPixelLocations(jnd,1)
        hnd = IPixelLocations(jnd,2)
        RImageToWrite(gnd,hnd) = RSimulatedPatterns(ind,IThicknessIndex,jnd)
     END DO
	 
	 !Apply blur again, temp fix until all subroutines combined into one
	 IF (IImageProcessingFLAG.EQ.4) THEN
       Rradius=0.95_RKIND!!!*+*+ will need to be added as a line in felix.inp +*+*!!!
       CALL BlurG(RImageToWrite,Rradius,IErr)
	 END IF

     CALL WriteReflectionImage(IChOutWIImage,&
          RImageToWrite,IErr,2*IPixelCount,2*IPixelCount,2_IKIND)       
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationImages(", my_rank, ") error in WriteReflectionImage()"
        RETURN
     ENDIF
     
     CLOSE(IChOutWIImage,IOSTAT=IErr)       
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationImages(", my_rank, ") error Closing Reflection Image()"
        RETURN
     ENDIF
     
  END DO

     CALL WriteIterationStructure(path,IErr) 
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationOutput(",my_rank,")error in WriteIterationStructure"
        RETURN
     ENDIF

     CALL WriteStructureFactors(path,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationOutput(",my_rank,")error in WriteStructureFactors()"
        RETURN
     ENDIF

     CALL WriteOutVariables(Iter,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationOutput(",my_rank,")error in WriteOutVariables()"
        RETURN
     ENDIF
  END IF
  
END SUBROUTINE WriteIterationOutput

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE WriteIterationImages(path,IThicknessIndex,IErr)
!now redundant
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,jnd,hnd,gnd,IThicknessIndex
  REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: RImageToWrite
  CHARACTER*200,INTENT(IN) :: path

  DO ind = 1,INoOfLacbedPatterns
     CALL OpenReflectionImage(IChOutWIImage,path,IErr,ind,2*IPixelCount,2_IKIND)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationImages(",my_rank,") error in OpenReflectionImage()"
        RETURN
     ENDIF

     RImageToWrite = ZERO
     DO jnd = 1,IPixelTotal
        gnd = IPixelLocations(jnd,1)
        hnd = IPixelLocations(jnd,2)
        RImageToWrite(gnd,hnd) = RSimulatedPatterns(ind,IThicknessIndex,jnd)
     END DO

     CALL WriteReflectionImage(IChOutWIImage,&
          RImageToWrite,IErr,2*IPixelCount,2*IPixelCount,2_IKIND)       
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationImages(", my_rank, ") error in WriteReflectionImage()"
        RETURN
     ENDIF
     
     CLOSE(IChOutWIImage,IOSTAT=IErr)       
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationImages(", my_rank, ") error Closing Reflection Image()"
        RETURN
     ENDIF
     
  END DO

END SUBROUTINE WriteIterationImages

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE WriteIterationStructure(path,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,jnd
  CHARACTER*200,INTENT(IN) :: path
  CHARACTER*200 :: SPrintString,filename,fullpath

!!$  Write out non symmetrically related atomic positions

  WRITE(filename,*) "Structure.cif"
  WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))

  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',FILE=TRIM(ADJUSTL(fullpath)))
 !RB
  WRITE(IChOutSimplex,FMT='(A16)') "data_felixrefine"
  WRITE(IChOutSimplex,FMT='(A5)') "loop_"
  WRITE(IChOutSimplex,FMT='(A14,1X,F9.6)') "_cell_length_a",RLengthX
  WRITE(IChOutSimplex,FMT='(A14,1X,F9.6)') "_cell_length_b",RLengthY
  WRITE(IChOutSimplex,FMT='(A14,1X,F9.6)') "_cell_length_c",RLengthZ
  WRITE(IChOutSimplex,FMT='(A17,1X,F9.6)') "_cell_angle_alpha",RAlpha*180/PI
  WRITE(IChOutSimplex,FMT='(A16,1X,F9.6)') "_cell_angle_beta",RBeta*180/PI
  WRITE(IChOutSimplex,FMT='(A17,1X,F9.6)') "_cell_angle_gamma",RGamma*180/PI
  WRITE(IChOutSimplex,FMT='(A30,1X,A10)') "_symmetry_space_group_name_H-M '",SSpaceGrp,"'"
  WRITE(IChOutSimplex,FMT='(A5)') " "
  WRITE(IChOutSimplex,FMT='(A5)') "loop_"
  WRITE(IChOutSimplex,FMT='(A22)') "_atom_site_type_symbol"
!  WRITE(IChOutSimplex,FMT='(A25)') "_atom_site_Wyckoff_symbol"
  WRITE(IChOutSimplex,FMT='(A18)') "_atom_site_fract_x"
  WRITE(IChOutSimplex,FMT='(A18)') "_atom_site_fract_y"
  WRITE(IChOutSimplex,FMT='(A18)') "_atom_site_fract_z"
!  WRITE(IChOutSimplex,FMT='(A25)') "_atom_site_B_iso_or_equiv"
!  WRITE(IChOutSimplex,FMT='(A20)') "_atom_site_occupancy"

  DO jnd = 1,SIZE(RBasisAtomPosition,DIM=1)!RB only gives refined atoms, needs work
     WRITE(IChOutSimplex,FMT='(A2,1X,3(F9.6,1X))') &
	 SBasisAtomName(jnd),RBasisAtomPosition(jnd,:)
!     WRITE(IChOutSimplex,FMT='(A2,1X,A1,1X,3(F9.6,1X),F5.3,1X,F5.3)') &
!	 SBasisAtomName(jnd),SWyckoffSymbols(jnd),RBasisAtomPosition(jnd,:), &
!	 RBasisIsoDW(jnd),RBasisOccupancy(jnd)
  END DO
  WRITE(IChOutSimplex,FMT='(A22)') "#End of refinement cif"
  
  CLOSE(IChOutSimplex)

!!$  Write out full atomic positions

!XX  WRITE(filename,*) "StructureFull.txt"
!XX  WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))
!XXPRINT*,"RAtomPosition,SAtomName"  
!XX  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
!XX        FILE=TRIM(ADJUSTL(fullpath)))
!XX    DO jnd = 1,SIZE(RAtomPosition,DIM=1)
!XX     WRITE(IChOutSimplex,FMT='(A2,1X,3(F9.6,1X))') SAtomName(jnd),RAtomPosition(jnd,1:3)
!XX    END DO
!XX  CLOSE(IChOutSimplex)

END SUBROUTINE WriteIterationStructure

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE WriteOutVariables(Iter,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,IStart,IEnd,jnd,ITotalOutputVariables
  INTEGER(IKIND),INTENT(IN) :: Iter
  CHARACTER*200 :: SFormat,STotalOutputVariables
  INTEGER(IKIND),DIMENSION(IRefinementVariableTypes) :: IOutputVariables
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RDataOut

  ! Need to Determine total no. of variables to be written out, this is different from the no. of refinement variables
  
  IOutputVariables(1) =  IRefineMode(1) * &
       2*INoofUgs+1 ! Structure Factors are Complex so require two output variables each     
  IOutputVariables(2) = IRefineMode(2) * & !Structural Coordinates
       (SIZE(RBasisAtomPosition,DIM=1) * SIZE(RBasisAtomPosition,DIM=2))
  IOutputVariables(3) = &
       IRefineMode(3) * & !Atomic Site Occupancies
       SIZE(RBasisOccupancy,DIM=1)
  IOutputVariables(4) = &
       IRefineMode(4) * & !Isotropic Debye Waller Factors
       SIZE(RBasisIsoDW,DIM=1)
  IOutputVariables(5) = &
       IRefineMode(5) * & !Anisotropic Debye Waller Factors
       SIZE(RAnisotropicDebyeWallerFactorTensor)
  IOutputVariables(6) = &    
       IRefineMode(6) * 3 !Lattice Parameters (a,b,c) 
  IOutputVariables(7) = &
       IRefineMode(7) * 3 !Lattice Angles (alpha,beta,gamma)
  IOutputVariables(8) = & 
       IRefineMode(8) !Convergence angle
  IOutputVariables(9) = &
       IRefineMode(9) !Absorption
  IOutputVariables(10) = &
       IRefineMode(10) !Accelerating Voltage
  IOutputVariables(11) = &
       IRefineMode(11) !Residual Sum of Squares Scaling Factor
  IOutputVariables(12) =  IRefineMode(12) * &
       2*INoofUgs+1 ! Structure Factors are Complex so require two output variables each     
  
  ITotalOutputVariables = SUM(IOutputVariables) ! Total Output
  
  ALLOCATE(RDataOut(ITotalOutputVariables),STAT=IErr)

  DO jnd = 1,IRefinementVariableTypes
     IF(IRefineMode(jnd).EQ.0) THEN
        CYCLE !The refinement variable type is not being refined, skip
     END IF
     IF(jnd.EQ.1) THEN!It's an atom coordinate refinement
        IStart = 1
     ELSE
        IStart = SUM(IOutputVariables(1:(jnd-1)))+1
     END IF
     IEND = SUM(IOutputVariables(1:jnd))

     SELECT CASE(jnd)
     CASE(1)
        DO ind = 1,INoofUgs
           IStart = (ind*2)-1
           IEnd = ind*2
           RDataOut(IStart:IEnd) = [REAL(CUniqueUg(ind+IUgOffset)), REAL(AIMAG(CUniqueUg(ind+IUgOffset)),RKIND)]
        END DO
		RDataOut(IEnd+1) = RAbsorptionPercentage!RB last variable is absorption
     CASE(2)
        RDataOut(IStart:IEnd) = RESHAPE(TRANSPOSE(RBasisAtomPosition),SHAPE(RDataOut(IStart:IEnd)))
     CASE(3)
        RDataOut(IStart:IEnd) = RBasisOccupancy
     CASE(4)
        RDataOut(IStart:IEnd) = RBasisIsoDW
     CASE(5)
        RDataOut(IStart:IEnd) = RESHAPE(RAnisotropicDebyeWallerFactorTensor,SHAPE(RDataOut(IStart:IEnd)))
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
           RDataOut(IStart:IEnd) = [REAL(CUniqueUg(ind+IUgOffset)), REAL(AIMAG(CUniqueUg(ind+IUgOffset)),RKIND)]
        END DO
		RDataOut(IEnd+1) = RAbsorptionPercentage!RB last variable is absorption
    END SELECT
  END DO

  WRITE(STotalOutputVariables,*) ITotalOutputVariables
  WRITE(SFormat,*) "(I5.1,1X,F13.9,1X,"//TRIM(ADJUSTL(STotalOutputVariables))//"(F13.9,1X))"

  OPEN(UNIT=IChOutSimplex,file='IterationLog.txt',form='formatted',status='unknown',position='append')

  WRITE(UNIT=IChOutSimplex,FMT=SFormat) Iter,RCrossCorrelation,RDataOut

  CLOSE(IChOutSimplex)

END SUBROUTINE WriteOutVariables

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE WriteStructureFactors(path,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind
  CHARACTER*200,INTENT(IN) :: path
  CHARACTER*200 :: filename,fullpath

  WRITE(filename,*) "StructureFactors.txt"
  WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))
  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
       FILE=TRIM(ADJUSTL(fullpath)))

  DO ind = 1,nReflections
     WRITE(IChOutSimplex,FMT='(3I5.1,2F13.9)') NINT(Rhkl(ind,:)),CUgMat(ind,1)
  END DO

  CLOSE(IChOutSimplex)

END SUBROUTINE WriteStructureFactors

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
