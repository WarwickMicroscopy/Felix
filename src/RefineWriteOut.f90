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

  INTEGER(IKIND) :: &
       IErr,ind,IStart,IEnd,jnd, &
       ITotalOutputVariables
  INTEGER(IKIND),INTENT(IN) :: &
       IIterationCount
  CHARACTER*200 :: &
       SFormat,STotalOutputVariables
  INTEGER(IKIND),DIMENSION(IRefinementVariableTypes) :: &
       IOutputVariables
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       RDataOut

  ! Need to Determine total no. of variables to be written out, this is different from the no. of refinement variables
  
  IOutputVariables(1) =  IRefineModeSelectionArray(1) * &
       2*SIZE(CUgMat,DIM=1) ! Structure Factors are Complex so require two output variables each     
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
  
  ALLOCATE(&
       RDataOut(&
       ITotalOutputVariables), &
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
        DO ind = 1,SIZE(CUgMat,DIM=1)
           IStart = (ind*2)-1
           IEnd = ind*2
           RDataOut(IStart:IEnd) = [REAL(REAL(CUgMat(ind,1)),RKIND), REAL(AIMAG(CUgMat(ind,1)),RKIND)]
        END DO
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

SUBROUTINE WriteIterationOutput(IIterationCount,IThicknessIndex,IExitFlag,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,IIterationCount,IThickness
  INTEGER(IKIND),INTENT(IN) :: &
       IThicknessIndex,IExitFLAG
  CHARACTER*200 :: &
       path
  
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
     ENDIF
     
     DEALLOCATE( &
          Rhkl,&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteIterationOutput(", my_rank, ") error ", IErr, &
             " in Deallocation Rhkl"
        RETURN
     ENDIF

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


SUBROUTINE WriteStructureFactors(path,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind
  CHARACTER*200,INTENT(IN) :: &
       path
  CHARACTER*200 :: &
       filename,fullpath

  WRITE(filename,*) "StructureFactors.txt"
  WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))

  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
       FILE=TRIM(ADJUSTL(fullpath)))
  DO ind = 1,SIZE(CUgMat,DIM=1)
     WRITE(IChOutSimplex,FMT='(3I5.1,2F13.9)') NINT(RHKL(ind,:)),CUgMat(ind,1)
  END DO

  CLOSE(IChOutSimplex)

END SUBROUTINE WriteStructureFactors

SUBROUTINE WriteIterationStructure(path,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,jnd
  CHARACTER*200,INTENT(IN) :: &
       path
  CHARACTER*200 :: &
       SPrintString,filename,fullpath

!!$  Write out non symmetrically related atomic positions

  WRITE(filename,*) "StructureCif.txt"
  WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))

  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
        FILE=TRIM(ADJUSTL(fullpath)))
  
  DO jnd = 1,SIZE(RAtomSiteFracCoordVec,DIM=1)
     WRITE(IChOutSimplex,FMT='(A2,1X,3(F9.6,1X))') SAtomName(jnd),RAtomSiteFracCoordVec(jnd,:)
  END DO
  
  CLOSE(IChOutSimplex)

!!$  Write out full atomic positions

  CALL ExperimentalSetup(IErr)!RB Really? wtf?!?
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in ExperimentalSetup"
     RETURN
  ENDIF

  WRITE(filename,*) "StructureFull.txt"
  WRITE(fullpath,*) TRIM(ADJUSTL(path)),'/',TRIM(ADJUSTL(filename))

  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
        FILE=TRIM(ADJUSTL(fullpath)))
  
  DO jnd = 1,SIZE(MNP,DIM=1)
     WRITE(IChOutSimplex,FMT='(A2,1X,3(F9.6,1X))') SMNP(jnd),MNP(jnd,1:3)
  END DO
  
  CLOSE(IChOutSimplex)
  
  DEALLOCATE( &
       MNP,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation MNP"
     RETURN
  ENDIF
      
  DEALLOCATE( &
        SMNP, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation SMNP"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RFullAtomicFracCoordVec, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RFullAtomicFracCoordVec"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       SFullAtomicNameVec,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation SFullAtomicNameVec"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IFullAnisotropicDWFTensor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation IFullAnisotropicDWFTensor"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IFullAtomNumber,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation IFullAtomNumber"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RFullIsotropicDebyeWallerFactor,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RFullIsotropicDebyeWallerFactor"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RFullPartialOccupancy,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RFullPartialOccupancy"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RDWF,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RDWF"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       ROcc,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation ROcc"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IAtoms,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation IAtoms"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       IAnisoDWFT,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation IAnisoDWFT"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RgVecMatT,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMatT"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RgVecMag,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMag"
     RETURN
  ENDIF
       
  DEALLOCATE( &
       RgVecVec,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, &
          " in Deallocation RgVecMag"
     RETURN
  ENDIF

  DEALLOCATE( &
       RrVecMat, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteIterationStructure(", my_rank, ") error ", IErr, " in DEALLOCATE of RrVecMat"
     RETURN
  ENDIF

END SUBROUTINE WriteIterationStructure

SUBROUTINE WriteIterationImages(path,IThicknessIndex,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind,jnd,hnd,gnd,IThicknessIndex
  REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: &
       RImage
  CHARACTER*200,INTENT(IN) :: &
       path

  DO ind = 1,IReflectOut
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
