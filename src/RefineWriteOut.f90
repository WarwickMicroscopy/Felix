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


SUBROUTINE WriteOutVariables(IIterationCount,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,IStart,IEnd,jnd,ITotalOutputVariables
  INTEGER(IKIND),INTENT(IN) :: IIterationCount
  CHARACTER*200 :: SFormat,STotalOutputVariables
  INTEGER(IKIND),DIMENSION(IRefinementVariableTypes) :: IOutputVariables
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RDataOut

  ! Need to Determine total no. of variables to be written out, this is different from the no. of refinement variables
  
  IOutputVariables(1) =  IRefineModeSelectionArray(1) * &
       2*nReflections ! Structure Factors are Complex so require two output variables each     
  IOutputVariables(2) = IRefineModeSelectionArray(2) * & !Structural Coordinates
       (SIZE(RAtomSiteFracCoordVec,DIM=1) * SIZE(RAtomSiteFracCoordVec,DIM=2))
  IOutputVariables(3) = &
       IRefineModeSelectionArray(3) * & !Atomic Site Occupancies
       SIZE(RAtomicSitePartialOccupancy,DIM=1)
  IOutputVariables(4) = &
       IRefineModeSelectionArray(4) * & !Isotropic Debye Waller Factors
       SIZE(RIsotropicDebyeWallerFactors,DIM=1)
  IOutputVariables(5) = &
       IRefineModeSelectionArray(5) * & !Anisotropic Debye Waller Factors
       SIZE(RAnisotropicDebyeWallerFactorTensor)
  IOutputVariables(6) = &    
       IRefineModeSelectionArray(6) * 3 !Lattice Parameters (a,b,c) 
  IOutputVariables(7) = &
       IRefineModeSelectionArray(7) * 3 !Lattice Angles (alpha,beta,gamma)
  IOutputVariables(8) = & 
       IRefineModeSelectionArray(8) !Convergence angle
  IOutputVariables(9) = &
       IRefineModeSelectionArray(9) !Absorption
  IOutputVariables(10) = &
       IRefineModeSelectionArray(10) !Accelerating Voltage
  IOutputVariables(11) = &
       IRefineModeSelectionArray(11) !Residual Sum of Squares Scaling Factor
  
  ITotalOutputVariables = SUM(IOutputVariables) ! Total Output
  
  ALLOCATE(RDataOut(ITotalOutputVariables), &
       STAT=IErr)

  DO jnd = 1,IRefinementVariableTypes

     IF(IRefineModeSelectionArray(jnd).EQ.0) THEN
        CYCLE !The refinement variable type is not being refined, skip
     END IF
     
     IF(jnd.EQ.1) THEN
        IStart = 1
     ELSE
        IStart = SUM(IOutputVariables(1:(jnd-1)))+1
     END IF

     IEND = SUM(IOutputVariables(1:jnd))

     SELECT CASE(jnd)
     CASE(1)
        DO ind = 1,15!RB nReflections
           IStart = (ind*2)-1
           IEnd = ind*2
           RDataOut(IStart:IEnd) = [REAL(REAL(CUgMat(ind,1)),RKIND), REAL(AIMAG(CUgMat(ind,1)),RKIND)]
        END DO
		RDataOut(IEnd+1) = RAbsorptionPercentage!RIndependentVariableValues(2*INoofUgs+1)!RB last variable is absorption
     CASE(2)
        RDataOut(IStart:IEnd) = RESHAPE(TRANSPOSE(RAtomSiteFracCoordVec),SHAPE(RDataOut(IStart:IEnd)))
     CASE(3)
        RDataOut(IStart:IEnd) = RAtomicSitePartialOccupancy
     CASE(4)
        RDataOut(IStart:IEnd) = RIsotropicDebyeWallerFactors
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
     END SELECT
  END DO

  WRITE(STotalOutputVariables,*) ITotalOutputVariables
  WRITE(SFormat,*) "(I5.1,1X,F13.9,1X,"//TRIM(ADJUSTL(STotalOutputVariables))//"(F13.9,1X))"

  OPEN(UNIT=IChOutSimplex,file='IterationLog.txt',form='formatted',status='unknown',position='append')

  WRITE(UNIT=IChOutSimplex,FMT=SFormat) IIterationCount, RCrossCorrelation,RDataOut

  CLOSE(IChOutSimplex)

END SUBROUTINE WriteOutVariables

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE WriteIterationOutput(IIterationCount,IThicknessIndex,IExitFlag,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,IIterationCount,IThickness
  INTEGER(IKIND),INTENT(IN) :: IThicknessIndex,IExitFLAG
  CHARACTER*200 :: path
  
  IF(IExitFLAG.EQ.1.OR.(IIterationCount.GE.(IPreviousPrintedIteration+IPrint))) THEN
     

     IThickness = (RInitialThickness + (IThicknessIndex-1)*RDeltaThickness)/10!RB in nm 
     
!RB     WRITE(path,"(A10,I5.5,A3,I1.1,I1.1,I1.1,I1.1,A2,I5.5,A2,I5.5,A2,I5.5)") &    
     WRITE(path,"(A10,I4.4,A1,I3.3,A3,I3.3,A1,I3.3)") &
          "Iteration",IIterationCount,&
!RB          "-f-",&
!RB          IScatterFactorMethodFLAG, &
!RB          IZolzFLAG, &
!RB          IAbsorbFLAG, &
!RB          IAnisoDebyeWallerFactorFlag,&
          "_",IThickness,&
          "nm_",2*IPixelcount,&
          "x",2*IPixelcount
     
     call system('mkdir ' // path)
     
     PRINT*,"IExitFLAG = ",IExitFLAG,"; there have been",&
          IIterationCount-IPreviousPrintedIteration,"iterations from my last print"
     
     IPreviousPrintedIteration = IIterationCount
      
     CALL WriteIterationImages(path,IThicknessIndex,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationOutput(", my_rank, ") error ", IErr, &
             " in WriteIterationImages"
        RETURN
     END IF

     CALL WriteIterationStructure(path,IErr) 
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationOutput(", my_rank, ") error ", IErr, &
             " in WriteIterationStructure"
        RETURN
     ENDIF
           
     CALL WriteStructureFactors(path,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ",IErr,"in WriteStructureFactors()"
        RETURN
     ENDIF

     CALL WriteOutVariables(IIterationCount,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Felixfunction(", my_rank, ") error ",IErr,"in WriteOutVariables()"
        RETURN
     ENDIF
  END IF
  
END SUBROUTINE WriteIterationOutput

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
        FILE=TRIM(ADJUSTL(fullpath)))
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

  DO jnd = 1,SIZE(RAtomSiteFracCoordVec,DIM=1)!RB only gives refined atoms, needs work
     WRITE(IChOutSimplex,FMT='(A2,1X,3(F9.6,1X))') &
	 SAtomName(jnd),RAtomSiteFracCoordVec(jnd,:)
!     WRITE(IChOutSimplex,FMT='(A2,1X,A1,1X,3(F9.6,1X),F5.3,1X,F5.3)') &
!	 SAtomName(jnd),SWyckoffSymbols(jnd),RAtomSiteFracCoordVec(jnd,:), &
!	 RIsotropicDebyeWallerFactors(jnd),RAtomicSitePartialOccupancy(jnd)
  END DO
  WRITE(IChOutSimplex,FMT='(A22)') "#End of refinement cif"
  
  CLOSE(IChOutSimplex)

!!$  Write out full atomic positions

  !CALL ExperimentalSetup(IErr)!RB Really? wtf?!?
  !IF( IErr.NE.0 ) THEN
  !   PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
  !        " in ExperimentalSetup"
  !   RETURN
  !ENDIF

!XX  WRITE(filename,*) "StructureFull.txt"
!XX  WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))
!XXPRINT*,"MNP,SMNP"  
!XX  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
!XX        FILE=TRIM(ADJUSTL(fullpath)))
!XX    DO jnd = 1,SIZE(MNP,DIM=1)
!XX     WRITE(IChOutSimplex,FMT='(A2,1X,3(F9.6,1X))') SMNP(jnd),MNP(jnd,1:3)
!XX    END DO
!XX  CLOSE(IChOutSimplex)

END SUBROUTINE WriteIterationStructure

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE WriteIterationImages(path,IThicknessIndex,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,jnd,hnd,gnd,IThicknessIndex
  REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: RImage
  CHARACTER*200,INTENT(IN) :: path

  DO ind = 1,INoOfLacbedPatterns
     CALL OpenReflectionImage(IChOutWIImage,path,IErr,ind,2*IPixelCount,2_IKIND)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationImages(", my_rank, ") error in OpenReflectionImage()"
        RETURN
     ENDIF
     
     RImage = ZERO
     DO jnd = 1,IPixelTotal
        gnd = IPixelLocations(jnd,1)
        hnd = IPixelLocations(jnd,2)
        RImage(gnd,hnd) = RIndividualReflections(ind,IThicknessIndex,jnd)
     END DO
     
     CALL WriteReflectionImage(IChOutWIImage,&
          RImage,IErr,2*IPixelCount,2*IPixelCount,2_IKIND)       
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
