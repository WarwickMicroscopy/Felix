!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! FelixSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
!
! (C) 2013/14, all right reserved
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This file is part of FelixSim.
!
!  FelixSim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  FelixSim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with FelixSim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: inout.f90,v 1.59 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Log: inout.f90,v $
! Revision 1.59  2014/04/28 12:26:19  phslaz
! Fixed minor bugs with the new reflection pool request
!
! Revision 1.58  2014/04/24 17:24:33  phslaz
! ReflectionPool, StrongBeams and WeakBeams are now chosen pixel by pixel, also fixed an odd bug with CountTotalAtoms where it would double count duplicated positions, Sim and Draw make and run
!
! Revision 1.57  2014/04/09 17:06:30  phslaz
! IImageFLAG = 2 now saves amplitude and phase images
!
! Revision 1.55  2014/03/27 21:39:14  phslaz
! Added two new flags IImageFlag and IOutputFlag
! IImageFLAG = 0 (Montage) 1 (Montage and Reflections)
! IOutputFLAG = 0 (nothing) 1 (EigenSpectra) 2 (UgMat) 3 (Wavefunctions)
! Have also put many Print statments into IWriteflag control
! code compiles and runs
!
! Revision 1.53  2014/03/27 16:01:02  phsht
! BER->Felix
!
! Revision 1.52  2014/03/27 14:55:48  phslaz
! Writing now works again, still not MPI
!
! Revision 1.49  2014/03/25 19:17:52  phslaz
! Individual images and montages now right out in correctly generated folders
!
! Revision 1.47  2014/03/25 15:35:34  phsht
! included consistent start of each source file and GPL statements
!
! Revision 1.46  2014/03/24 12:50:31  phslaz
! fixed write flag 2 and rewrote lacbed to read file once instead of every thickness
!
! Revision 1.43  2014/03/07 11:01:15  phsht
! added numbers to correct calculation of IMAXCBuffer in WriteEigenSystem_MPI()
!
! Revision 1.42  2014/03/07 10:49:44  phslaz
! Corrected issues with inpcif, should now work with badly structured cifs
!
! Revision 1.41  2014/02/21 15:26:42  phslaz
! collapsed a few parts of main.f90 into subroutines in util.f90
!
! Revision 1.39  2014/02/20 13:49:56  phslaz
! removed old write statement
!
! Revision 1.38  2014/02/20 13:34:40  phsht
! Ug matrix now written via MPI
! changed format of EigenSystem file to first have EVal and then EVec
!
! Revision 1.37  2014/02/20 13:17:31  phsht
! removed WF/WI/EX/EV/etc outputs from BER-BLOCH
! combined EX+EV output into ES output and made MPI compatible
! restructured ES file
!
! Revision 1.36  2014/02/17 14:08:55  phslaz
! Lacbed now works
!
! Revision 1.35  2014/02/07 14:33:05  phslaz
! LACBED code now reads eigen spectra output
!
! Revision 1.34  2014/02/04 15:17:41  phsht
! added code to input the EigenSystem
!
! Revision 1.33  2014/01/24 15:13:41  phsht
! new .sca file with scattering factors for up to element 103
!
! Revision 1.32  2014/01/17 16:57:27  phslaz
! InpCif now reads in isotropic debye waller factors but there are not currently used
!
! Revision 1.31  2014/01/16 16:12:42  phsht
! work on scattering factors
!
! Revision 1.30  2014/01/10 16:24:13  phslaz
! Renamed Crystal, Orthogonal and Microscope reference frame variables C,O and M respectively
!
! Revision 1.29  2014/01/09 18:04:43  phslaz
! ZOLZ and Absorption installed
!
! Revision 1.28  2013/11/29 16:46:20  phsht
! removed thickness loop again, was unwieldly
!
! Revision 1.27  2013/11/27 12:30:21  phsht
! more consistency in MPI output routines and their error handling
!
! Revision 1.26  2013/11/27 10:21:48  phsht
! now all output files work with MPI
!
! Revision 1.25  2013/11/26 22:25:56  phsht
! MPIIO for the rest of the files, but something not right yet
!
! Revision 1.24  2013/11/26 18:27:18  phsht
! MPI file IO now working for WI (intensities)
!
! Revision 1.23  2013/11/25 18:26:33  phsht
! progress on the MPI version
!
! Revision 1.22  2013/11/25 10:04:01  phsht
! Scattering factors now in BER.sca
!
! Revision 1.21  2013/11/15 16:03:33  phsht
! files for CIFTBX to read in .CIF files
!
! Revision 1.20  2013/11/13 17:48:02  phsht
! 1st attempt at LACBED code
!
! Revision 1.19  2013/10/21 15:56:31  phslaz
! Changed Variable names
!
! Revision 1.18  2013/10/03 15:48:57  phsht
! added routine HKLMake() in util.f90 to make BWM compatible;
! checked for 64 pixels with input file of this version
!
! Revision 1.17  2013/10/03 13:41:41  phsht
! typo in inout.f90 let to wrong thickness values being used;
! now works again
!
! Revision 1.16  2013/10/03 11:15:16  phsht
! new input file structure
!
! Revision 1.15  2013/10/02 20:43:28  phsht
! replaced old input names with new names;
! structure of input file still needs changing
!
! Revision 1.14  2013/09/23 16:52:29  phsht
! OPEN, write and CLOSE of data files now implemented
!
! Revision 1.13  2013/09/09 13:24:32  phsht
! InputScatteringFactors() now works
!
! Revision 1.12  2013/09/09 10:58:59  phsht
! subroutines in util.f90 now seem to work
!
! Revision 1.11  2013/09/04 07:26:30  phsht
! minor reformatting of code
!
! Revision 1.10  2013/09/03 15:24:39  phslaz
! Added The Beam Selection Criteria
!
! Revision 1.9  2013/09/02 15:42:42  phsht
! code checked with matlab version up to and including BraggCentral
!
! Revision 1.8  2013/09/02 15:13:14  phsht
! added two more flags to the input file and inout.f90
!
! Revision 1.7  2013/06/11 14:53:08  phsht
! more work in the Diffraction defs part
!
! Revision 1.6  2013/06/11 09:56:10  phsht
! InputHKL() works
!
! Revision 1.5  2013/06/10 15:08:01  phsht
! work this morning
!
! Revision 1.4  2013/06/10 08:20:28  phsht
! more work before realizing that MATLAB file is not quite the latest version
!
! Revision 1.3  2013/06/07 07:57:32  phsht
! some constants and parameters defined
!
! Revision 1.2  2013/06/07 07:11:28  phsht
! more ground work
!
! Revision 1.1  2013/04/03 19:35:45  phsht
! first installation of basic Fortran routines/structure
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! -----------------------------------------------------------------------
!Input: Read the input file
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE Input( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER IErr, ILine
  REAL(KIND=RKIND) ROfIter
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: Input()"
  END IF

  OPEN(UNIT= IChInp, ERR= 120, FILE= "Felix.inp",&
       STATUS= 'OLD')
  ILine= 1

  ! ----------------------------------------------------------------------
  ! introductory comment lines

  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ! ----------------------------------------------------------------------
  ! BLOCH method input

  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)') 
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ! ----------------------------------------------------------------------
  ! control flags

  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IWriteFLAG
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IWriteFLAG = ", IWriteFLAG
  END IF
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IImageFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IImageFLAG = ", IImageFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IOutputFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IOutputFLAG = ", IOutputFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IBinorTextFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IBinorTextFLAG = ", IBinorTextFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IScatterFactorMethodFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IScatterFactorMethodFLAG = ", IScatterFactorMethodFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) ICentralBeamFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ICentralBeamFLAG = ", ICentralBeamFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMaskFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IMaskFLAG = ", IMaskFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IZolzFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IZolzFLAG = ", IZolzFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IAbsorbFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IAbsorbFLAG = ", IAbsorbFLAG
  ENDIF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IAnisoDebyeWallerFactorFlag
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IAnisoDebyeWallerFactorFlag = ", IAnisoDebyeWallerFactorFlag
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IBeamConvergenceFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IBeamConvergenceFLAG = ", IBeamConvergenceFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IPseudoCubicFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IPseudoCubicFLAG = ", IPseudoCubicFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IXDirectionFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IXDirectionFLAG = ", IXDirectionFLAG
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IDevFLAG
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IDevFLAG = ", IDevFLAG
  END IF

  ! ----------------------------------------------------------------------
  ! beam details
  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IPixelCount
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IPixelCount = ", ILine,IPixelCount
  END IF

  ! ----------------------------------------------------------------------
  ! beam selection criteria
  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMinReflectionPool
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IMinReflectionPool = ", IMinReflectionPool
  END IF
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMinStrongBeams
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IMinStrongBeams = ", IMinStrongBeams
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IMinWeakBeams
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IMinWeakBeams = ", IMinWeakBeams
  END IF

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RBSBMax
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RBSBMax = ", RBSBMax
  END IF

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RBSPMax
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RBSPMax = ", RBSPMax
  END IF

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RConvergenceTolerance
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RConvergenceTolerance = ", RConvergenceTolerance
  END IF

  ! ----------------------------------------------------------------------
  ! crystal settings
  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RDebyeWallerConstant
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RDebyeWallerConstant = ", RDebyeWallerConstant
  END IF

  RMeanSquaredDisplacement= RDebyeWallerConstant/(8*PI**2)
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RAbsorptionPercentage
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RAbsorptionPercentage = ", RAbsorptionPercentage
  END IF

  ! ----------------------------------------------------------------------
  ! microscopy settings

  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RConvergenceAngle
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RConvergenceAngle = ", RConvergenceAngle
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IIncidentBeamDirectionX
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IIncidentBeamDirectionX = ", IIncidentBeamDirectionX
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IIncidentBeamDirectionY
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IIncidentBeamDirectionY = ", IIncidentBeamDirectionY
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IIncidentBeamDirectionZ
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IIncidentBeamDirectionZ = ", IIncidentBeamDirectionZ
  END IF

  RZDirC(1)= REAL(IIncidentBeamDirectionX,RKIND)
  RZDirC(2)= REAL(IIncidentBeamDirectionY,RKIND)
  RZDirC(3)= REAL(IIncidentBeamDirectionZ,RKIND)
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IXDirectionX
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IXDirectionX = ", IXDirectionX
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IXDirectionY
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IXDirectionY = ", IXDirectionY
  END IF
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IXDirectionZ
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IXDirectionZ = ", IXDirectionZ
  END IF

  RXDirC(1)= IXDirectionX
  RXDirC(2)= IXDirectionY
  RXDirC(3)= IXDirectionZ  

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) INormalDirectionX
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"INormalDirectionX = ", INormalDirectionX
  END IF

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) INormalDirectionY
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"INormalDirectionY = ", INormalDirectionY
  END IF
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) INormalDirectionZ
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"INormalDirectionZ = ", INormalDirectionZ
  END IF

  RNormDirC(1)= INormalDirectionX
  RNormDirC(2)= INormalDirectionY
  RNormDirC(3)= INormalDirectionZ
  
  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RAcceleratingVoltage

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RAcceleratingVoltage = ", RAcceleratingVoltage
  END IF

  !-----------------------------------------------------------------------
  ! Iterative Ug input

  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) INoofUgs
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"INoofUgs = ", INoofUgs
  END IF

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RPercentageUgChange
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RPercentageUgChange = ", RPercentageUgChange
  END IF


  ! ----------------------------------------------------------------------
  ! LACBED method input
  
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')
  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ! ----------------------------------------------------------------------
  ! sample thickness loop

  ILine= ILine+1; READ(IChInp,ERR=20,END=30,FMT='(A)')

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RInitialThickness
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RInitialThickness = ", RInitialThickness
  END IF

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RFinalThickness
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RFinalThickness = ", RFinalThickness
  END IF

  ILine= ILine+1
  READ(IChInp,15,ERR=20,END=30) RDeltaThickness
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RDeltaThickness = ", RDeltaThickness
  END IF
  
  ILine= ILine+1
  READ(IChInp,10,ERR=20,END=30) IReflectOut
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"IReflectOut = ", IReflectOut
  END IF

10 FORMAT(27X,I15.1)
  ! 10	FORMAT("IMAXIteration= ",I15.1)
15 FORMAT(27X,F18.9)
  ! 15     FORMAT("IMAXIteration= ",F18.9)
  
  CLOSE(IChInp, ERR=130)
  
  ! check the parameters for validity
   
  RETURN
  
  !	error in OPEN detected
120 PRINT*,"Input(): ERR in OPEN"
  GOTO 1000
  
  !	error in CLOSE detected
130 PRINT*,"Input(): ERR in CLOSE"
  GOTO 1000
  
  !	error in READ detected
20 PRINT*,"Input(): ERR in READ at line", ILine
  GOTO 1000
  
  !	EOF in READ occured prematurely
30 PRINT*,"Input(): EOF in READ at line", ILine
  
  ! dump the input help
  
1000 &
  PRINT*,"Input parameters:          ; explanation:"
  PRINT*,"--------------------------------------------------------------------"
  PRINT*,"IWriteFLAG          = 1          ; (12) 0/1/2/3/4 = no/log/category/wave fcn/RGamma/RHO output"
  
  IErr= 1
  RETURN
  
END SUBROUTINE Input


! -----------------------------------------------------------------------
!
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE InputScatteringFactors( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  CHARACTER*200 dummy
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       rdummy

  INTEGER IErr, ILine, ILength, ind
  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN

     PRINT*,"DBG: InputScatteringFactors()"
     
  END IF
  ILine= 0
  
!!$  OPEN(UNIT= IChInp, ERR= 120, FILE= "FelixDoyle.sca",&
!!$       STATUS= 'OLD')
  OPEN(UNIT= IChInp, ERR= 120, FILE= "Felix.sca",&
       STATUS= 'OLD')
!!$  
  DO
     READ(UNIT= IChInp, END=100, ERR=20, FMT='(a)') dummy
     ILine=ILine+1
  ENDDO
100 ILength=ILine
  
  IF(IWriteFLAG.GE.10) THEN
     
     PRINT*,"DBG: InputScatteringFactors(): ILength=", ILength
     
  END IF
  
  REWIND(UNIT=IChInp)
  IF(IWriteFLAG.GE.10) THEN
     
     PRINT*,"DBG: actual reading of data"
  END IF

  SELECT CASE (IScatterFactorMethodFLAG)

  CASE(0)

     ALLOCATE( &
          rdummy(1),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"InputScatteringFactors(): error in memory ALLOCATE()"
        RETURN
     ENDIF

     ALLOCATE(RScattFactors(ILength,12), STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"InputScatteringFactors(): error in memory ALLOCATE()"
        RETURN
     ENDIF

     DO ILine=1,ILength,1
        READ(UNIT= IChInp, FMT='(13(E15.11,1X))') &
             rdummy, RScattFactors(ILine,:)     
        !PRINT*,rdummy, RScattFactors(ILine,:)
     ENDDO

  CASE(1)

     ALLOCATE( &
          rdummy(13),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"InputScatteringFactors(): error in memory ALLOCATE()"
        RETURN
     ENDIF

     ALLOCATE(RScattFactors(ILength,8), STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"InputScatteringFactors(): error in memory ALLOCATE()"
        RETURN
     ENDIF
     
     DO ILine=1,ILength,1
        READ(UNIT= IChInp, FMT='(21(E15.11,1X))') &
             rdummy(:), RScattFactors(ILine,:)     
        !PRINT*,rdummy, RScattFactors(ILine,:)
     ENDDO

  CASE(2)

     ALLOCATE( &
          rdummy(21),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"InputScatteringFactors(): error in memory ALLOCATE()"
        RETURN
     ENDIF

     ALLOCATE(RScattFactors(ILength,8), STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"InputScatteringFactors(): error in memory ALLOCATE()"
        RETURN
     ENDIF
     
     DO ILine=1,ILength,1
        READ(UNIT= IChInp, FMT='(29(E15.11,1X))') &
             rdummy(:), RScattFactors(ILine,:)     
        !PRINT*,rdummy, RScattFactors(ILine,:)
     ENDDO

  END SELECT

  DEALLOCATE( &
       rdummy,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"InputScatteringFactors(): error in memory DEALLOCATE()"
     RETURN
  ENDIF  
  
  CLOSE(IChInp)
  
  RETURN
  
  !	error in OPEN detected
120 PRINT*,"InputScatteringFactors(): ERR in OPEN"
  GOTO 1000
  
  !	error in CLOSE detected
130 PRINT*,"InputScatteringFactors(): ERR in CLOSE"
  GOTO 1000
  
  !	error in READ detected
20 PRINT*,"InputScatteringFactors(): ERR in READ at line", ILine
  GOTO 1000
  
  !	EOF in READ occured prematurely
30 PRINT*,"InputScatteringFactors(): EOF in READ at line", ILine
  
1000 &
  PRINT*,"InputScatteringFactors expects a file such as"
  PRINT*,"--------------------------------------------------------------------"
  PRINT*,"22      0.737291729     9.96300042      0.355417565     9.96300042 ", &
       "0.496982599     0.072659586     0.016335034     0.360199179     1.42171303 ", &
       "     0.073173566     0.957232308      15.8512114"
  PRINT*,"23      0       0       0       0       0       0       0       0  ", &
       "     0       0       0       0"

  IErr= 1
  RETURN
END SUBROUTINE InputScatteringFactors
  
! -----------------------------------------------------------------------
!
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE ReadEigenSystemChunk( IAllocationChunk,IErr )
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara 
  USE IChannels

  USE MPI
  USE MyMPI
    
  IMPLICIT NONE

  CHARACTER*25 EIGENFORMAT
  INTEGER(IKIND) :: IRank, &
       Iindex, Ijndex, IWriteLine, IStrongBeamIndex,IInputBeams, &
       IAllocationChunk,InChunks,IindPrevious, IjndPrevious, &
       IErr, ILength, ind,IReadLine,INewLineFLAG,index
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: &
       REigenVectorRealTemp, REigenVectorImagTemp
  REAL(RKIND) :: &
       REigenValueRealTemp, REigenValueImagTemp

  PRINT*,"InputEigenSystem()"

  !-----------------------------------------------------------
  ! eigenvalues
  !-----------------------------------------------------------

  InBeams = 0
  ind = 0

  IF(IWriteFLAG.GE.10) THEN
     
     PRINT*,"DBG: actual reading of data"
  
  END IF
  IReadLine = 1
  INewLineFLAG = 1
  DO 
     READ(UNIT= IChInp, END=100, ERR=20, FMT='(7(I3.1,1X))',ADVANCE='NO') &
          IRank,Iindex, Ijndex, nReflections,&
          IInputBeams, IWriteLine, IStrongBeamIndex

     IF (IInputBeams.EQ.0) THEN
        IReadLine = IReadLine-1
        EXIT
     END IF
        
     InBeams(IReadLine) = IInputBeams
     ILACBEDStrongBeamList(IReadLine,IWriteLine) = IStrongBeamIndex
     IPixelLocation(IReadLine,1) = Iindex
     IPixelLocation(IReadLine,2) = Ijndex
      
     IF(INewLineFlag.EQ.1) THEN
        IPixelCountTotal = IPixelCountTotal + 1
        ALLOCATE( &
             REigenVectorRealTemp(IInputBeams),&
             REigenVectorImagTemp(IInputBeams),&
             STAT=IErr)
        IF(IErr.NE.0) THEN
           PRINT*,"ERROR IN ALLOCATE OF LACBED Temps"
           STOP
        ENDIF
        INewLineFlag=0
     ENDIF
     
     WRITE(EIGENFORMAT,*) IInputBeams

     READ(ICHInp,END=100,ERR=20,FMT="((1F13.10,1X),(1F13.10,1X),&
          "//TRIM(ADJUSTL(TRIM(EIGENFORMAT)))//"(1F13.10,1X), &
          "//TRIM(ADJUSTL(TRIM(EIGENFORMAT)))//"(1F13.10,1X))", &
          ADVANCE='YES') REigenValueRealTemp, REigenValueImagTemp,&     
          REigenVectorRealTemp,REigenVectorImagTemp
     CEigenVectorsChunk(IReadLine,IWriteLine,1:IInputBeams) = &
          REigenVectorRealTemp+REigenVectorImagTemp*CIMAGONE
     CEigenValuesChunk(IReadLine,IWriteLine) = &
          REigenValueRealTemp+ REigenValueImagTemp*CIMAGONE
     
     IF(IReadLine.EQ.IAllocationChunk.AND.IWriteLine.EQ.IInputBeams) EXIT
     IF(IWriteLine.EQ.IInputBeams) THEN
        IReadLine =  IReadLine + 1
        INewLineFLAG = 1
        DEALLOCATE(&
             REigenVectorRealTemp,&
             REigenVectorImagTemp)
     ENDIF
     
  ENDDO

  100 PRINT*,"Done with the reading already"

  RETURN
  
  !	error in OPEN detected
120 PRINT*,"InputEigenSystem(): ERR in OPEN"
  !GOTO 1000
  
  !	error in CLOSE detected
130 PRINT*,"InputEigenSystem(): ERR in CLOSE"
  !GOTO 1000
  
  !	error in READ detected
20 PRINT*,"InputEigenSystem(): ERR in READ at line", IReadLine
  !GOTO 1000
  
  !	EOF in READ occured prematurely
30 PRINT*,"InputEigenSystem(): EOF in READ at line", IReadLine
  
  IErr= 1
  RETURN
END SUBROUTINE ReadEigenSystemChunk
  
! --------------------------------------------------------------------
! OpenData
! --------------------------------------------------------------------


SUBROUTINE OpenData(IChOutWrite, prefix, surname, IErr)

  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  CHARACTER*27 surname, surnamelength
  CHARACTER*2 prefix,postfix
  INTEGER(KIND=IKIND) IChOutWrite, IErr

  CHARACTER*34 filename
  INTEGER index
  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN

     PRINT*,"OpenData()"

  END IF
  
  WRITE(surnamelength,*) LEN_TRIM(surname)

  WRITE(filename,"(A2,A2,A1,A"//TRIM(surnamelength)//",A4)") "F-",prefix,"-",surname,".txt"

  IF(IWriteFLAG.GE.10) THEN
     
     PRINT*,filename

  END IF
  IF (IWriteFLAG.GE.10) THEN
     SELECT CASE(IChOutWrite)
     CASE(IChOutWF)
        PRINT*, "OpenData: opening channel", IChOutWF, &
             "for WAVE FUNCTIONS (WF*.txt)"
     CASE(IChOutWI)
        PRINT*, "OpenData: opening channel", IChOutWI, &
             "for WAVE INTENSITIES (WI*.txt)"
     CASE(IChOutEV)
        PRINT*, "OpenData: opening channel", IChOutEV, &
             "for EIGENVALUES of UgMat (EV*.txt)"
     CASE(IChOutEX)
        PRINT*, "OpenData: opening channel", IChOutEX, &
             "for EIGENVECTORS of UgMat (EX*.txt)"
     CASE(IChOutUM)
        PRINT*, "OpenData: opening channel", IChOutUM, &
             "for UgMat (UM*.txt)"
     CASE DEFAULT
        PRINT*, "OpenData: opening UNKNOWN", IChOutWrite, &
             "channel "
     END SELECT
  END IF

  OPEN(UNIT= IChOutWrite, ERR= 10, STATUS= 'UNKNOWN', FILE=TRIM(filename))

  RETURN

  ! error in OPEN detected
10 PRINT*,"WriteDataC(): ERR in OPEN()"
  !PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenData
  
! --------------------------------------------------------------------
! OpenData
! --------------------------------------------------------------------


SUBROUTINE OpenImageForReadIn(IErr,filename)

  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(KIND=IKIND) IChOutWrite, IErr

  CHARACTER*34 filename

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN

     PRINT*,"OpenImageForReadIn()"

  END IF

  !filename = "Felix.img"

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,filename

  END IF

  OPEN(UNIT= IChInImage, ERR= 10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),FORM='UNFORMATTED',&
       ACCESS='DIRECT',IOSTAT=Ierr,RECL=IImageSizeXY(1)*8)

  RETURN

  ! error in OPEN detected
10 PRINT*,"WriteDataC(): ERR in OPEN()"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenImageForReadIn

SUBROUTINE ReadImageForRefinement(IErr)

  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN

     PRINT*,"ReadImageForRefinement()"

  END IF

  DO ind=1,IImageSizeXY(2)
     READ(IChInImage,rec=ind) RImageIn(ind,:)
  END DO
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadImageForRefinement (", my_rank, ") error in READ()",IErr
     RETURN
  ENDIF
  
END SUBROUTINE ReadImageForRefinement


! --------------------------------------------------------------------
! OpenData
! --------------------------------------------------------------------

SUBROUTINE OpenDataForAppend(IChOutWrite, prefix, surname, IErr)

  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  CHARACTER*27 surname, surnamelength
  CHARACTER*2 prefix,postfix
  INTEGER(KIND=IKIND) IChOutWrite, IErr

  CHARACTER*34 filename
  INTEGER index

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"OpenData()"

  END IF
  
  WRITE(surnamelength,*) LEN_TRIM(surname)
  
  WRITE(filename,"(A2,A2,A1,A"//TRIM(surnamelength)//",A4)") "F-",prefix,"-",surname,".txt"

  IF (IWriteFLAG.GE.10) THEN
     SELECT CASE(IChOutWrite)
     CASE(IChOutWF)
        PRINT*, "OpenData: opening channel", IChOutWF, &
             "for WAVE FUNCTIONS (WF*.txt)"
     CASE(IChOutWI)
        PRINT*, "OpenData: opening channel", IChOutWI, &
             "for WAVE INTENSITIES (WI*.txt)"
     CASE(IChOutEV)
        PRINT*, "OpenData: opening channel", IChOutEV, &
             "for EIGENVALUES of UgMat (EV*.txt)"
     CASE(IChOutEX)
        PRINT*, "OpenData: opening channel", IChOutEX, &
             "for EIGENVECTORS of UgMat (EX*.txt)"
     CASE(IChOutUM)
        PRINT*, "OpenData: opening channel", IChOutUM, &
             "for UgMat (UM*.txt)"
     CASE DEFAULT
        PRINT*, "OpenData: opening UNKNOWN", IChOutWrite, &
             "channel "
     END SELECT
  END IF
  
  OPEN(UNIT= IChOutWrite, ERR= 10, STATUS= 'UNKNOWN',&
       FILE=TRIM(filename),ACCESS='APPEND')

  RETURN

  ! error in OPEN detected
10 PRINT*,"WriteDataC(): ERR in OPEN()"
  !PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenDataForAppend

! --------------------------------------------------------------------
! Open Reflection Image
! --------------------------------------------------------------------

SUBROUTINE OpenReflectionImage(IChOutWrite, surname, IErr,IReflectWriting,IImageSizeX)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  USE IChannels

  CHARACTER*27 surname
  CHARACTER*20 prefix,postfix,h,k,l
  INTEGER(KIND=IKIND) IChOutWrite, IErr,IReflectWriting,IImageSizeX

  CHARACTER*50 filename
  CHARACTER*40 fileext
  INTEGER index

  IF ((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"OpenReflectionImage()"
  END IF
  

  SELECT CASE(IChOutWrite)
  CASE(MontageOut)
  CASE DEFAULT
     WRITE(h,*)  NINT(RHKL(IReflectWriting,1))
     WRITE(k,*)  NINT(RHKL(IReflectWriting,2))
     WRITE(l,*)  NINT(RHKL(IReflectWriting,3))
  END SELECT

  IF (IWriteFLAG.GE.10) THEN
     PRINT*,filename
  END IF
  
  SELECT CASE (IBinorTextFLAG)
  CASE(0)
     WRITE(fileext,*) TRIM(ADJUSTL(".bin")) 
  CASE(1)
     WRITE(fileext,*) TRIM(ADJUSTL(".txt"))
  END SELECT
  
  SELECT CASE(IChOutWrite)
  CASE(IChOutWFImageReal)        
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"/F-WF-A_",&
          TRIM(ADJUSTL(h)),&
          TRIM(ADJUSTL(k)),&
          TRIM(ADJUSTL(l)),&
          TRIM(ADJUSTL(fileext))
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE FUNCTION REAL PART (WR*.txt)"
     END IF
  CASE(IChOutWFImagePhase)        
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"/F-WF-P_",&
          TRIM(ADJUSTL(h)),&
          TRIM(ADJUSTL(k)),&
          TRIM(ADJUSTL(l)),&
          TRIM(ADJUSTL(fileext))
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE FUNCTION PHASE PART (WP*.txt)"
     END IF
  CASE(IChOutWIImage)        
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"/F-WI_",&
          TRIM(ADJUSTL(h)),&
          TRIM(ADJUSTL(k)),&
          TRIM(ADJUSTL(l)),&
          TRIM(ADJUSTL(fileext))
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE INTENSITIES"
     END IF
  CASE(MontageOut)        
     WRITE(filename,*) "F-WI-",TRIM(ADJUSTL(surname)),&
          TRIM(ADJUSTL(fileext))
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE INTENSITIES"
     END IF
  CASE DEFAULT
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening UNKNOWN channel ", IChOutWrite
     END IF
  END SELECT
  
  
  SELECT CASE (IBinorTextFLAG)
     CASE(0)
        OPEN(UNIT=IChOutWrite, ERR=10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),FORM='UNFORMATTED',&
             ACCESS='DIRECT',IOSTAT=Ierr,RECL=IImageSizeX*8)
        
     CASE(1)
        OPEN(UNIT=IChOutWrite, ERR=10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)))
     END SELECT
     RETURN
     
  ! error in OPEN detected
10 PRINT*,"OpenReflectionImage(): ERR in OPEN()",IErr
  IErr= 1
  RETURN
  
END SUBROUTINE OpenReflectionImage

!-----------------------------------------------------------------
! Write Reflection Images
!-----------------------------------------------------------------

SUBROUTINE WriteReflectionImage( IChOutWrite, data, IErr,IImageSizeX,IImageSizeY)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  INTEGER(KIND=IKIND) IErr,IImageSizeY,IImageSizeX
  REAL(KIND=RKIND),DIMENSION(IImageSizeX,IImageSizeY) :: data
  CHARACTER*100 CSizeofData
  INTEGER ind, IChOutWrite
  CHARACTER*100 SFormatString

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"WriteReflectionImage"
  END IF

  SELECT CASE (IBinorTextFLAG)
     
  CASE(0)
     
     DO ind = 1,(IImageSizeY)
        WRITE(IChOutWrite,rec=ind) data(ind,:)
     END DO

  CASE(1)
     
     DO ind = 1,(2*IPixelCount)
        WRITE(CSizeofData,*) 2*IPixelCount
        WRITE(SFormatString,*) "("//TRIM(ADJUSTL(CSizeofData))//"(1F6.3,1X),A1)"
        WRITE(IChOutWrite,FMT=SFormatString,ERR=20) data(ind,:)
     END DO
     
  END SELECT

  RETURN
  ! error in WRITE detected
20 PRINT*,"WriteReflectionImage(): ERR in WRITE()",Ierr
  IErr= 1
  RETURN

  ! buffer too short
30 PRINT*,"WriteImageR_MPI(): EOF in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteReflectionImage

! --------------------------------------------------------------------
! WriteDataC

SUBROUTINE WriteDataC( IChOutWrite, ipos,jpos, data, size, step, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MyMPI

  INTEGER(KIND=IKIND) ipos,jpos, size, step, IErr
  COMPLEX(KIND=CKIND) data(size)

  INTEGER ind, IChOutWrite
  DO ind=1,size,step
     IF (ABS(data(ind)).GE. TINY) THEN
        WRITE(IChOutWrite,100) my_rank, ipos,jpos,ind, DBLE(data(ind)),AIMAG(data(ind)), &
             ARG(DBLE(data(ind)),AIMAG(data(ind)))
     ELSE
        WRITE(IChOutWrite,100) my_rank, ipos,jpos,ind, DBLE(data(ind)),AIMAG(data(ind)), &
             0.0D0
     ENDIF
     
100  FORMAT(4I4,3(G25.15))
  ENDDO
  
  RETURN
  
  ! error in WRITE detected
20 PRINT*,"WriteDataC(): ERR in WRITE()"
  IErr= 1
  RETURN
  
END SUBROUTINE WriteDataC

!---------------------------------------------------------------------
!Output WriteEigenSystem_MPI
!---------------------------------------------------------------------

SUBROUTINE WriteEigenSystem_MPI( IChOutWrite, &
     ipos,jpos,nReflect,nbeamout, CdataEVal,CdataEVec,ISbeamlist, rows, cols, step, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  INTEGER(KIND=IKIND) ipos,jpos, step,rows,cols ,IErr,nReflect,nbeamout
  COMPLEX(KIND=CKIND), DIMENSION(nbeamout,nbeamout) :: CdataEVec
  COMPLEX(CKIND), DIMENSION(nbeamout) :: CdataEVal
  INTEGER(IKIND), DIMENSION(nbeamout) :: ISbeamlist

  INTEGER my_status(MPI_STATUS_SIZE)

  CHARACTER*25 FORMATstring
  CHARACTER(IMAXCBuffer) DATAstring

  INTEGER ind, IChOutWrite
  
  IF ( 2*13*SIZE(CdataEVec)+2*13*SIZE(CdataEVal)+3*SIZE(ISbeamlist)+3*6*ADD_OUT_INFO .GT. IMAXCBuffer ) THEN
     IErr=1
     PRINT*, "WriteEigenSystem_MPI(", my_rank, ") error ", IErr, &
          ", output buffer size IMAXCBuffer=", IMAXCBuffer, &
          " smaller than SIZEOF(EVec+EVal+Sbeam)=", &
          SIZEOF(CdataEVec)+SIZEOF(CdataEVal)+SIZEOF(ISbeamlist), &
          SIZE(CdataEVec) + SIZE(CdataEVal) + SIZE(ISbeamlist), &
          2*13*SIZE(CdataEVec)+2*13*SIZE(CdataEVal)+3*SIZE(ISbeamlist)+3*6*ADD_OUT_INFO
     GOTO 20
  ENDIF
  
  DO ind=1,rows,step
     
     WRITE(FORMATstring,*) nbeamout
     WRITE(DATAstring,&
          "(7(I3.1,1X),"//TRIM(ADJUSTL(TRIM(FORMATstring)))//"(1F13.10,1X), &
          "//TRIM(ADJUSTL(TRIM(FORMATstring)))//"(1F13.10,1X),(1F13.10,1X),(1F13.10,1X),A1)") &
          my_rank, ipos,jpos,nReflect,nbeamout,ind, ISbeamlist(ind), &
          REAL(CdataEVal(ind)), AIMAG(CdataEVal(ind)),&
          REAL(CdataEVec(ind,:)), AIMAG(CdataEVec(ind,:)), &
          CHAR(10)
     
     CALL MPI_FILE_WRITE_SHARED(IChOutWrite, TRIM(DATAstring), LEN_TRIM(DATAstring), &
          MPI_CHARACTER, my_status, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteEigenSystem_MPI(", my_rank, ") error ", IErr, &
             " in MPI_FILE_WRITE_SHARED() for file handle ",IChOutWrite
        RETURN
     ENDIF
     
  ENDDO
  
  RETURN

  ! error in WRITE detected
20 PRINT*,"WriteEigenSystem_MPI(): ERR in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteEigenSystem_MPI
!---------------------------------------------------------------------
!Output WriteEigenSystem_MPI
!---------------------------------------------------------------------

SUBROUTINE WriteEigenSystemBinary_MPI( IChOutWrite, &
     ipos,jpos,nReflect,nbeamout, CdataEVal,CdataEVec,ISbeamlist, rows, cols, step, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  INTEGER(KIND=IKIND) ipos,jpos, step,rows,cols ,IErr,nReflect,nbeamout
  COMPLEX(KIND=CKIND), DIMENSION(nbeamout,nbeamout) :: CdataEVec
  COMPLEX(CKIND), DIMENSION(nbeamout) :: CdataEVal
  INTEGER(IKIND), DIMENSION(nbeamout) :: ISbeamlist

  INTEGER my_status(MPI_STATUS_SIZE)

  CHARACTER*25 FORMATstring
  CHARACTER(IMAXCBuffer) DATAstring
  
  INTEGER ind, IChOutWrite
  
  IF ( 2*13*SIZE(CdataEVec)+2*13*SIZE(CdataEVal)+3*SIZE(ISbeamlist)+3*6*ADD_OUT_INFO .GT. IMAXCBuffer ) THEN
     IErr=1
     PRINT*, "WriteEigenSystem_MPI(", my_rank, ") error ", IErr, &
          ", output buffer size IMAXCBuffer=", IMAXCBuffer, &
          " smaller than SIZEOF(EVec+EVal+Sbeam)=", &
          SIZEOF(CdataEVec)+SIZEOF(CdataEVal)+SIZEOF(ISbeamlist), &
          SIZE(CdataEVec) + SIZE(CdataEVal) + SIZE(ISbeamlist), &
          2*13*SIZE(CdataEVec)+2*13*SIZE(CdataEVal)+3*SIZE(ISbeamlist)+3*6*ADD_OUT_INFO
     GOTO 20
  ENDIF

  DO ind=1,rows,step

     CALL MPI_FILE_WRITE_SHARED(IChOutWrite, TRIM(DATAstring), LEN_TRIM(DATAstring), &
          MPI_CHARACTER, my_status, IErr)

     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteEigenSystem_MPI(", my_rank, ") error ", IErr, &
             " in MPI_FILE_WRITE_SHARED() for file handle ",IChOutWrite
        RETURN
     ENDIF

  ENDDO
  
  RETURN

  ! error in WRITE detected
20 PRINT*,"WriteEigenSystem_MPI(): ERR in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteEigenSystemBinary_MPI

! --------------------------------------------------------------------
! WriteDataR
!--------------------------------------------------------------------

SUBROUTINE WriteDataR( IChOutWrite, ipos,jpos, data, size, step, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MyMPI

  INTEGER(KIND=IKIND) ipos,jpos, size, step, IErr
  REAL(KIND=RKIND) data(size)

  INTEGER ind, IChOutWrite

  DO ind=1,size,step
     WRITE(IChOutWrite,100) my_rank, ipos,jpos,ind,data(ind)
100  FORMAT(4I4,G25.15)
  ENDDO
     
  RETURN

  ! error in WRITE detected
20 PRINT*,"WriteDataR(): ERR in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteDataR

! --------------------------------------------------------------------
! OpenData_MPI

SUBROUTINE OpenDataForAppend_MPI(IChOutWrite, prefix, surname, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE IChannels

  USE MPI
  USE MyMPI

  CHARACTER*27 surname
  CHARACTER*2 prefix,postfix
  INTEGER IChOutWrite, IErr

  CHARACTER*35 filename
  INTEGER index

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"OpenDataForAppend_MPI()"

  END IF

  WRITE(filename,'(A2,A2,A1,A12,A4)') "F-",prefix,"-",surname,".txt"
  
  CALL MPI_FILE_OPEN( MPI_COMM_WORLD, TRIM(filename), &
       MPI_MODE_CREATE+MPI_MODE_APPEND, MPI_INFO_NULL, IChOutWrite, IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"OpenDataForAppend_MPI(", my_rank, ") error ", IErr, &
          " in MPI_FILE_OPEN()"
     STOP
  ENDIF
  
  IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenData: opened channel", IChOutWrite, &
             "for ", filename
  END IF

  RETURN

  ! error in OPEN detected
10 PRINT*,"OpenDataForAppend_MPI(): ERR in OPEN()"
  !PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenDataForAppend_MPI

! --------------------------------------------------------------------
! OpenData_MPI

SUBROUTINE OpenData_MPI(IChOutWrite, prefix, surname, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE IChannels

  USE MPI
  USE MyMPI

  CHARACTER*27 surname
  CHARACTER*2 prefix,postfix
  INTEGER IChOutWrite, IErr

  CHARACTER*35 filename
  INTEGER index

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"OpenData_MPI()"
  END IF

  WRITE(filename,'(A2,A2,A1,A12,A4)') "F-",prefix,"-",surname,".txt"
  
  CALL MPI_FILE_OPEN( MPI_COMM_WORLD, TRIM(filename), &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, IChOutWrite, IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"OpenDataMPI(", my_rank, ") error ", IErr, &
          " in MPI_FILE_OPEN()"
     STOP
  ENDIF
  
  IF (IWriteFLAG.GE.10) THEN
     PRINT*, "OpenData: opened channel", IChOutWrite, &
          "for ", filename
  ENDIF
  RETURN
  
  ! error in OPEN detected
10 PRINT*,"OpenDatMPI(): ERR in OPEN()"
  !PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenData_MPI

! --------------------------------------------------------------------
! OpenImage_MPI

SUBROUTINE OpenImage_MPI(IChOutWrite, surname, IErr,IReflectWriting)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE IChannels

  USE MPI
  USE MyMPI

  CHARACTER*27 surname
  CHARACTER*2 prefix,postfix
  INTEGER IChOutWrite, IErr

  CHARACTER*100 filename,h,k,l
  INTEGER index,IReflectWriting

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"OpenImage_MPI()"
     
  END IF
  
  WRITE(h,*)  NINT(RHKL(IReflectWriting,1))
  WRITE(k,*)  NINT(RHKL(IReflectWriting,2))
  WRITE(l,*)  NINT(RHKL(IReflectWriting,3))

  WRITE(filename,*) TRIM(ADJUSTL(surname)),"/F-WI_",&
       TRIM(ADJUSTL(h)),&
       TRIM(ADJUSTL(k)),&
       TRIM(ADJUSTL(l)),&
       ".txt"
  
  CALL MPI_FILE_OPEN( MPI_COMM_WORLD, TRIM(filename), &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, IChOutWrite, IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"OpenImage_MPI(", my_rank, ") error ", IErr, &
          " in MPI_FILE_OPEN()"
     STOP
  ENDIF
  
  IF (IWriteFLAG.GE.10) THEN
     PRINT*, "OpenData: opened channel", IChOutWrite, &
          "for ", filename
  END IF
  
  RETURN
  
  ! error in OPEN detected
10 PRINT*,"OpenDatMPI(): ERR in OPEN()"
  !PRINT*, "file ", filename, " does not exist --- REOPEN not possible!"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenImage_MPI

! --------------------------------------------------------------------
! WriteDataRMPI

SUBROUTINE WriteDataR_MPI( IChOutWrite, ipos,jpos, data, size, step, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  INTEGER(KIND=IKIND) ipos,jpos, size, step, IErr
  INTEGER my_status(MPI_STATUS_SIZE)
  REAL(KIND=RKIND) data(size)

  CHARACTER(IMAXRBuffer) outstring
  INTEGER ind, IChOutWrite

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"WriteDataR_MPI()"
  END IF
  ! line by line

  IF ( SIZEOF(data)+ADD_OUT_INFO .GT. IMAXRBuffer ) THEN
     IErr=1
     PRINT*, "WriteDataR_MPI(", my_rank, ") error ", IErr, &
          ", output buffer size IMAXRBuffer=", IMAXRBuffer, &
          " smaller than SIZEOF(data)=", SIZEOF(data)
     GOTO 20
  ENDIF

  WRITE(outstring,*,ERR=20) my_rank,ipos,jpos,size,(data(ind),ind=1,size), CHAR(10)

  CALL MPI_FILE_WRITE_SHARED(IChOutWrite, TRIM(outstring), LEN_TRIM(outstring), &
       MPI_CHARACTER, my_status, IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteDataRMPI(", my_rank, ") error ", IErr, &
          " in MPI_FILE_WRITE_SHARED()"
     STOP
  ENDIF
  
  GOTO 9

  ! element by element
  DO ind=1,size

     CALL MPI_FILE_WRITE_SHARED(IChOutWrite, TRIM(outstring), LEN_TRIM(outstring), &
          MPI_CHARACTER, my_status, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteDataR_MPI(", my_rank, ") error ", IErr, &
             " in MPI_FILE_WRITE_SHARED() for file handle ",IChOutWrite
        RETURN
     ENDIF
  END DO

9 RETURN

  ! error in WRITE detected
20 PRINT*,"WriteDataR_MPI(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! buffer too short
30 PRINT*,"WriteDataR_MPI(): EOF in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteDataR_MPI! --------------------------------------------------------------------
! WriteDataRMPI

SUBROUTINE WriteImageR_MPI( IChOutWrite, data, IErr,ILocalPixelCountMin,ILocalPixelCountMax)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  INTEGER(KIND=IKIND) IErr,ILocalPixelCountMin,ILocalPixelCountMax
  INTEGER my_status(MPI_STATUS_SIZE)
  REAL(KIND=RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: data

  CHARACTER(IMAXRBuffer) outstring,CSizeofData
  INTEGER ind, IChOutWrite
  CHARACTER*100 SFormatString
  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"WriteImageR_MPI"
  END IF

  IF ( 6*SIZE(data)+ADD_OUT_INFO .GT. IMAXRBuffer ) THEN
     IErr=1
     PRINT*, "WriteImageR_MPI(", my_rank, ") error ", IErr, &
          ", output buffer size IMAXRBuffer=", IMAXRBuffer, &
          " smaller than SIZEOF(data)=", SIZE(data)
  ENDIF

  DO ind = 1,(2*IPixelCount)
     WRITE(CSizeofData,*) 2*IPixelCount
     WRITE(SFormatString,*) "("//TRIM(ADJUSTL(CSizeofData))//"(1F6.3,1X),A1)"
     WRITE(outstring,FMT=SFormatString,ERR=20) data(ind,:), CHAR(10)

     CALL MPI_FILE_WRITE_SHARED(IChOutWrite,TRIM(outstring), &
          LEN_TRIM(outstring),MPI_CHARACTER, my_status, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"WriteImageRMPI(", my_rank, ") error ", IErr, &
             " in MPI_FILE_WRITE_SHARED()"
        STOP
     ENDIF
  END DO
  
  RETURN
  ! error in WRITE detected
20 PRINT*,"WriteImageR_MPI(): ERR in WRITE()"
  IErr= 1
  RETURN

  ! buffer too short
30 PRINT*,"WriteImageR_MPI(): EOF in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteImageR_MPI

! --------------------------------------------------------------------
! WriteDataCMPI

SUBROUTINE WriteDataC_MPI( IChOutWrite, ipos,jpos, Cdata, Isize, step, IErr)

  USE MyNumbers

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara
  USE CPara
  USE SPara

  USE MPI
  USE MyMPI

  INTEGER(KIND=IKIND) ipos,jpos, Isize, step, IErr
  INTEGER my_status(MPI_STATUS_SIZE)
  COMPLEX(CKIND), DIMENSION(ISize):: &
       Cdata
  CHARACTER(IMAXCBuffer) outstring
  INTEGER ind, IChOutWrite,ITotalSize
  CHARACTER*150 SFormatstring,Sisize
  
  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"WriteDataC_MPI()"
  END IF
  
  WRITE(sisize,*) isize
  
  ! line by line

  IF ( (2*13*SIZE(Cdata)+3*6*ADD_OUT_INFO ).GT. IMAXCBuffer ) THEN
     IErr=1
     PRINT*, "WriteDataC_MPI(", my_rank, ") error ", IErr, &
          ", output buffer size IMAXCBuffer=", IMAXCBuffer, &
          " smaller than SIZEOF(data)=", SIZE(Cdata)
     GOTO 20
  ENDIF

  WRITE(SFormatString,*) &
       "(4(I7.1,1X),"//TRIM(ADJUSTL(TRIM(Sisize)))//"(1F13.10,1X), &
       "//TRIM(ADJUSTL(TRIM(sisize)))//"(1F13.10,1X),A1)"
  WRITE(outstring,FMT=SFormatString, ERR=20) my_rank,ipos,jpos,Isize,&
       REAL(Cdata,RKIND),AIMAG(Cdata), CHAR(10)
  CALL MPI_FILE_WRITE_SHARED(IChOutWrite, TRIM(outstring), LEN_TRIM(outstring), &
       MPI_CHARACTER, my_status, IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"WriteDataC_MPI(", my_rank, ") error ", IErr, &
          " in MPI_FILE_WRITE_SHARED() for file handle ",IChOutWrite
     RETURN
  ENDIF
  RETURN

  ! error in WRITE detected
20 PRINT*,"WriteDataC_MPI(): ERR in WRITE()"
  IErr= 1
  RETURN

END SUBROUTINE WriteDataC_MPI

