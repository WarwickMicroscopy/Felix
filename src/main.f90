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
! $Id: main.f90,v 1.89 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Log: main.f90,v $
! Revision 1.89  2014/04/28 12:26:19  phslaz
! Fixed minor bugs with the new reflection pool request
!
! Revision 1.87  2014/04/14 16:51:12  phslaz
! Seemingly fixed the rhombahedral problem, turns out theres was a mistake in inpcif where the 3rd angle was being read in incorrectly, have also written a new hklmake which is more understandable and applies selection rules directly rather than mathematically
!
! Revision 1.86  2014/04/09 17:06:30  phslaz
! IImageFLAG = 2 now saves amplitude and phase images
!
! Revision 1.84  2014/03/27 21:39:14  phslaz
! Added two new flags IImageFlag and IOutputFlag
! IImageFLAG = 0 (Montage) 1 (Montage and Reflections)
! IOutputFLAG = 0 (nothing) 1 (EigenSpectra) 2 (UgMat) 3 (Wavefunctions)
! Have also put many Print statments into IWriteflag control
! code compiles and runs
!
! Revision 1.82  2014/03/27 18:23:17  phsht
! solved MPI_REDUCE problem by adding the IErr as last argument
!
! Revision 1.81  2014/03/27 18:14:17  phslaz
! Conflicts in header resolved
!
! Revision 1.80  2014/03/27 16:26:28  phsht
! removed makefile.GF-LAPACK
!
! Revision 1.79  2014/03/27 14:55:48  phslaz
! Writing now works again, still not MPI
!
! Revision 1.75  2014/03/25 15:37:30  phsht
! more on GPL
!
! Revision 1.74  2014/03/25 15:35:34  phsht
! included consistent start of each source file and GPL statements
!
! Revision 1.73  2014/03/13 18:10:32  phslaz
! Seg Fault due to 32bit suspected constraint
!
! Revision 1.72  2014/03/07 10:49:45  phslaz
! Corrected issues with inpcif, should now work with badly structured cifs
!
! Revision 1.69  2014/02/20 13:34:40  phsht
! Ug matrix now written via MPI
! changed format of EigenSystem file to first have EVal and then EVec
!
! Revision 1.68  2014/02/20 13:17:31  phsht
! removed WF/WI/EX/EV/etc outputs from BER-BLOCH
! combined EX+EV output into ES output and made MPI compatible
! restructured ES file
!
! Revision 1.67  2014/02/17 14:08:56  phslaz
! Lacbed now works
!
! Revision 1.65  2014/01/23 18:50:55  phslaz
! Absorption and 2Pi convention installed and checked with matlab
!
! Revision 1.64  2014/01/20 20:20:44  phslaz
! Added a check for the existence of the debye waller factor in the cif file and use the constant stated int he input file if its not there
!
! Revision 1.63  2014/01/20 18:33:59  phslaz
! Isotropic Debye Waller Factor and Atomic Site Partial Occupancy done
!
! Revision 1.62  2014/01/16 16:12:42  phsht
! work on scattering factors
!
! Revision 1.61  2014/01/13 13:46:20  phslaz
! *** empty log message ***
!
! Revision 1.58  2014/01/07 17:11:52  phslaz
! Bug Fix : Hklmake was being called before the basevector assignment and as such, being called with lots of zeros...
!
! Revision 1.57  2014/01/07 11:49:17  phsht
! transformation to microscope reference included, needs testing
!
! Revision 1.56  2013/12/19 16:30:27  phsht
! new version of HKLMake(), using .cif information
!
! Revision 1.55  2013/12/19 14:58:57  phsht
! symmetry operations now correctly interpreted from .cif;
! structure correctly read in
!
! Revision 1.54  2013/11/29 16:46:20  phsht
! removed thickness loop again, was unwieldly
!
! Revision 1.53  2013/11/27 12:30:21  phsht
! more consistency in MPI output routines and their error handling
!
! Revision 1.52  2013/11/27 10:21:48  phsht
! now all output files work with MPI
!
! Revision 1.51  2013/11/26 22:25:56  phsht
! MPIIO for the rest of the files, but something not right yet
!
! Revision 1.50  2013/11/26 18:27:18  phsht
! MPI file IO now working for WI (intensities)
!
! Revision 1.49  2013/11/25 20:10:36  phsht
! added MPI split in main pixel loop
!
! Revision 1.48  2013/11/25 18:26:33  phsht
! progress on the MPI version
!
! Revision 1.47  2013/11/22 23:53:28  phsht
! added MPI commands
!
! Revision 1.46  2013/11/15 16:03:33  phsht
! files for CIFTBX to read in .CIF files
!
! Revision 1.45  2013/11/14 14:55:52  phsht
! included weak beams perturbatively
!
! Revision 1.44  2013/11/13 17:48:02  phsht
! 1st attempt at LACBED code
!
! Revision 1.43  2013/11/13 15:44:04  phsht
! fixed an error in the diagonal elements of CUgMatEffective which
! happens after strong beam selection; code now agrees again
! with matlab
!
! Revision 1.42  2013/11/08 17:51:41  phsht
! included strong beam selection; seems to work fast, but needs checking
!
! Revision 1.41  2013/10/30 13:00:44  phslaz
! Rearranged Variable declarations, allocations and deallocations
!
! Revision 1.40  2013/10/21 15:56:31  phslaz
! Changed Variable names
!
! Revision 1.39  2013/10/03 15:48:57  phsht
! added routine HKLMake() in util.f90 to make BWM compatible;
! checked for 64 pixels with input file of this version
!
! Revision 1.38  2013/10/03 14:06:30  phsht
! outcommented a DBG line
!
! Revision 1.37  2013/10/03 13:56:38  phsht
! replaced nReflections with nBeams in wave function calculation;
! led to errors, will have to be worked on later
!
! Revision 1.36  2013/10/03 13:41:41  phsht
! typo in inout.f90 let to wrong thickness values being used;
! now works again
!
! Revision 1.35  2013/10/03 12:52:25  phsht
! two errors corrected for beam selection parameter
!
! Revision 1.34  2013/10/03 11:15:16  phsht
! new input file structure
!
! Revision 1.33  2013/10/02 20:43:28  phsht
! replaced old input names with new names;
! structure of input file still needs changing
!
! Revision 1.32  2013/10/02 10:53:35  phsht
! minor formatting change
!
! Revision 1.31  2013/10/02 08:43:03  phsht
! removed (ind,jnd) declarations in a number of variables;
! reproduces previous results and with higher accuracy
!
! Revision 1.30  2013/10/01 20:38:01  phsht
! trying to write as a single loop over all pixels
!
! Revision 1.29  2013/09/24 16:29:10  phslaz
! Ugmateffective Fixed Now, images now match matlab :)
!
! Revision 1.28  2013/09/23 16:52:29  phsht
! OPEN, write and CLOSE of data files now implemented
!
! Revision 1.27  2013/09/19 11:19:33  phsht
! eigenvector calculation now working
!
! Revision 1.26  2013/09/18 16:08:42  phsht
! work with Keiths to include changes done with Richard yesterday and
! also to continue until intensities; reached these, but don't quite agree
! yet. Nevertheless, up to the diagonalization, things do work! nd eigenvalues
! are also correct, just eigenvectors not yet.
!
! Revision 1.25  2013/09/17 17:00:21  phsht
! diagonalization by LAPACK routines now included
!
! Revision 1.24  2013/09/11 16:35:12  phslaz
! Added RGamma Values and Wavefunction Calculations DOESNT compile due to wrongly indexed 4D matrix and missing invert function (its INV in matlab)
!
! Revision 1.23  2013/09/11 09:39:27  phslaz
! Added Calculation of 2KSh into UgMatEffective
!
! Revision 1.22  2013/09/10 16:51:54  phsht
! before implementation of UgMatEffective
!
! Revision 1.21  2013/09/10 14:19:13  phsht
! works up to and including BigK
!
! Revision 1.20  2013/09/10 09:00:58  phslaz
! Began calculating Potentials, (it'll not compile yet)
!
! Revision 1.19  2013/09/09 13:43:16  phsht
! gMatMat/gMatMag defined and working
!
! Revision 1.18  2013/09/09 12:15:25  phsht
! allocated Mask correctly
!
! Revision 1.17  2013/09/09 10:58:59  phsht
! subroutines in util.f90 now seem to work
!
! Revision 1.16  2013/09/06 09:38:05  phslaz
! i think i was defining ind and jnd as the wrong type of integer?
!
! Revision 1.15  2013/09/06 09:33:01  phslaz
! I broke everything when i declared the new variables :)
!
! Revision 1.14  2013/09/06 09:24:51  phslaz
! Bug Fix : Actually defined all the variables i need for the mask definition
!
! Revision 1.13  2013/09/06 09:21:01  phslaz
! Added the creation of the circular mask (needed later) made allocatable in smodules, allocated in main
!
! Revision 1.12  2013/09/05 15:10:33  phslaz
! the allocation og hklpositions(nReflections,2) is apparently bad syntax so thats now commented
!
! Revision 1.11  2013/09/05 15:08:23  phslaz
! hklpositions is causing issues and isnt massive important right now so ive commented it out for the time being
!
! Revision 1.10  2013/09/05 14:15:02  phslaz
! Bug Fix : think ive got it now, hadnt allocated any memory
!
! Revision 1.9  2013/09/04 15:26:00  phslaz
! Bug fix : Didnt allocate memory for gvecmag
!
! Revision 1.8  2013/09/02 15:13:14  phsht
! added two more flags to the input file and inout.f90
!
! Revision 1.7  2013/08/29 18:54:31  phsht
! minor stuff
!
! Revision 1.6  2013/06/11 14:53:08  phsht
! more work in the Diffraction defs part
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
!! Revision 1.1  2013/04/03 19:35:45  phsht
! first installation of basic Fortran routines/structure
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PROGRAM FelixSim
 
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  !--------------------------------------------------------------------
  ! local variable definitions
  !--------------------------------------------------------------------
  
  IMPLICIT NONE

  REAL(RKIND) time, norm
  COMPLEX(CKIND) sumC, sumD
  
  !--------------------------------------------------------------------
  ! image related variables	
  REAL(RKIND) Rx0,Ry0, RImageRadius,Rradius, Rthickness
  
  INTEGER(IKIND) ILocalPixelCountMin, ILocalPixelCountMax
  COMPLEX(RKIND) CVgij
  
  INTEGER(IKIND) ind,jnd,hnd,knd,pnd, &
       IHours,IMinutes,ISeconds
  INTEGER, DIMENSION(:), ALLOCATABLE :: &
       IWeakBeamVec,IRankArraySize,IRankArraySizeRoot
  REAL(RKIND),DIMENSION(:,:,:,:),ALLOCATABLE :: &
       RIndividualReflectionsRoot
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RFinalMontageImageRoot
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RImage
  COMPLEX(CKIND),DIMENSION(:,:,:,:), ALLOCATABLE :: &
       CAmplitudeandPhaseRoot
  INTEGER IRootArraySize, IPixelPerRank
  CHARACTER*40 surname, path
  CHARACTER*25 CThickness 
  CHARACTER*25 CThicknessLength
 
  INTEGER(IKIND),DIMENSION(2,2) :: ITest
  
  INTEGER(IKIND):: IErr, IThickness, IThicknessIndex, ILowerLimit, &
       IUpperLimit
  REAL(RKIND) StartTime, CurrentTime, Duration, TotalDurationEstimate

  !-------------------------------------------------------------------
  ! constants
  !-------------------------------------------------------------------

  CALL Init_Numbers
  
  !-------------------------------------------------------------------
  ! set the error value to zero, will change upon error
  !-------------------------------------------------------------------

  IErr=0

  !--------------------------------------------------------------------
  ! MPI initialization
  !--------------------------------------------------------------------

  ! Initialise MPI  
  CALL MPI_Init(IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in MPI_Init()"
     GOTO 9999
  ENDIF

  ! Get the rank of the current process
  CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in MPI_Comm_rank()"
     GOTO 9999
  ENDIF

  ! Get the size of the current communicator
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in MPI_Comm_size()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! protocal feature startup
  !--------------------------------------------------------------------
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"FelixSim: ", RStr,DStr,AStr, ", process ", my_rank, " of ", p
     PRINT*,"--------------------------------------------------------------"
  END IF

  !--------------------------------------------------------------------
  ! timing startup
  !--------------------------------------------------------------------

  CALL cpu_time(StartTime)

  !--------------------------------------------------------------------
  ! INPUT section
  !--------------------------------------------------------------------

  CALL Input( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in Input()"
     GOTO 9999
  ENDIF

  CALL InputScatteringFactors( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in InputScatteringFactors()"
     GOTO 9999
  ENDIF

  CALL InpCIF(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in InpCIF()"
     GOTO 9999
  ENDIF

  IF (ITotalAtoms.EQ.0) THEN
     CALL CountTotalAtoms(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error in CountTotalAtoms()"
        GOTO 9999
     ENDIF
  END IF
     
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ITotalAtoms = ",ITotalAtoms
  END IF
  
  !--------------------------------------------------------------------
  ! open outfiles
  !--------------------------------------------------------------------

  WRITE(surname,'(A1,I1.1,A1,I1.1,A1,I1.1,A2,I4.4)') &
       "S", IScatterFactorMethodFLAG, &
       "B", ICentralBeamFLAG, &
       "M", IMaskFLAG, &
       "_P", IPixelCount
  
  ! eigensystem
  IF(IOutputFLAG.GE.1) THEN
     CALL OpenData_MPI(IChOutES_MPI, "ES", surname, IErr)
  ENDIF
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in OpenData_MPI(EigenSystem)"
     GOTO 9999
  ENDIF
  
  ! UgMatEffective
  IF(IOutputFLAG.GE.2) THEN
     CALL OpenData_MPI(IChOutUM_MPI, "UM", surname, IErr)
  ENDIF
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in OpenDataMPI()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! Allocate Crystallography Variables
  !--------------------------------------------------------------------
       
  ALLOCATE( &
       RrVecMat(ITotalAtoms,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, " in ALLOCATE()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! microscopy settings
  !--------------------------------------------------------------------

  CALL MicroscopySettings( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in MicroscopySettings()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! crystallography settings
  !-------------------------------------------------------------------
  CALL Crystallography( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in Crystallography()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! diffraction initialization
  !--------------------------------------------------------------------

  CALL DiffractionPatternDefinitions( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in DiffractionPatternDefinitions()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! allocate memory for DYNAMIC variables according to nReflections
  !--------------------------------------------------------------------

  ! Image initialisation 
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: nReflections=", nReflections
  END IF
  
  ALLOCATE( &
       Rhklpositions(nReflections,2), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variable Rhklpositions"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! image initialization
  !--------------------------------------------------------------------

  CALL ImageInitialization( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in ImageInitializtion()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! define image masks
  !--------------------------------------------------------------------
      
  !Allocate Memory for Masking Image

  ALLOCATE( &
       RMask(2*IPixelCount,2*IPixelCount),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variable RMask"
     GOTO 9999
  ENDIF

  CALL ImageMask(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ImageMask"
     GOTO 9999
  END IF

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: main(", my_rank, ") IPixelTotal=", IPixelTotal
  END IF
 
  !--------------------------------------------------------------------
  ! MAIN section
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! Calculate Reflection Matrix
  !--------------------------------------------------------------------

  ALLOCATE( &  
       RgMatMat(nReflections,nReflections,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     GOTO 9999
  ENDIF
       
  ALLOCATE( &  
       RgMatMag(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     GOTO 9999
  ENDIF

  CALL GMatrixInitialisation (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in GMatrixInitialisation"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! calculating Ug matrix
  !--------------------------------------------------------------------

  !Allocate memory for Ug Matrix

  ALLOCATE( & 
       CUgMat(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     GOTO 9999
  ENDIF  

  ALLOCATE( & 
       CUgMatPrime(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Reflection Matrix"
     GOTO 9999
  ENDIF       

  CALL UgCalculation (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in UgCalculation"
     GOTO 9999
  ENDIF
  CALL UgAddAbsorption(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in UgCalculation"
     GOTO 9999
  ENDIF

  IF(IAbsorbFLAG.EQ.1) THEN
     CUgMat =  CUgMat+CUgMatPrime
  end IF

  Deallocate( &
       RgMatMat,RgMatMag,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in Deallocation of RgMat"
     GOTO 9999
  ENDIF
  
  !--------------------------------------------------------------------
  ! high-energy approximation (not HOLZ compatible)
  !--------------------------------------------------------------------
  
  RBigK= SQRT(RElectronWaveVectorMagnitude**2 + RMeanInnerCrystalPotential)

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: main(", my_rank, ") BigK=", RBigK
  END IF
  
         
  !--------------------------------------------------------------------
  ! reserve memory for effective eigenvalue problem
  !--------------------------------------------------------------------

  !Kprime Vectors and Deviation Parameter
  
  ALLOCATE( &
       RDevPara(nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RDevPara"
     GOTO 9999
  ENDIF

  ALLOCATE( &
       IStrongBeamList(nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables IStrongBeamList"
     GOTO 9999
  ENDIF

  ALLOCATE( &
       IWeakBeamList(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables IWeakBeamList"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! MAIN LOOP: solve for each (ind,jnd) pixel
  !--------------------------------------------------------------------

  ILocalPixelCountMin= (IPixelTotal*(my_rank)/p)+1
  ILocalPixelCountMax= (IPixelTotal*(my_rank+1)/p) 

  IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"main(", my_rank, "): starting the eigenvalue problem"
     PRINT*,"main(", my_rank, "): for lines ", ILocalPixelCountMin, &
          " to ", ILocalPixelCountMax
  ENDIF
  
  IThicknessCount= (RFinalThickness- RInitialThickness)/RDeltaThickness + 1
  
  IF(IImageFLAG.LE.1) THEN
     ALLOCATE( &
          RIndividualReflections(2*IPixelCount,&
          2*IPixelCount,IReflectOut,IThicknessCount),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Individual Images"
        GOTO 9999
     ENDIF
     
     RIndividualReflections = ZERO
  ELSE
     ALLOCATE( &
          CAmplitudeandPhase(2*IPixelCount,&
          2*IPixelCount,IReflectOut,IThicknessCount),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Amplitude and Phase"
        GOTO 9999
     ENDIF
     CAmplitudeandPhase = CZERO
  END IF
  
  ALLOCATE( &
       CFullWaveFunctions(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CFullWaveFunctions"
     GOTO 9999
  ENDIF
  
  ALLOCATE( &
       RFullWaveIntensity(nReflections), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RFullWaveIntensity"
     GOTO 9999
  ENDIF  

  IMAXCBuffer = 200000
  IPixelComputed= 0
  
  DEALLOCATE( &
       RScattFactors, &
       RrVecMat, Rsg, &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in DEALLOCATE() "
     GOTO 9999
  ENDIF
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.6) THEN
     PRINT*,"main(",my_rank,") Entering BlochLoop()"
  END IF

  DO knd = ILocalPixelCountMin,ILocalPixelCountMax,1
     ind = IPixelLocations(knd,1)
     jnd = IPixelLocations(knd,2)
     CALL BlochCoefficientCalculation(ind,jnd,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " in BlochCofficientCalculation"
        GOTO 9999
     ENDIF
  END DO
  
  IF((IWriteFLAG.GE.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"MAIN : ",my_rank," is exiting calculation loop"
  END IF

  !--------------------------------------------------------------------
  ! close outfiles
  !--------------------------------------------------------------------
  
  ! eigensystem
  IF(IOutputFLAG.GE.1) THEN
     CALL MPI_FILE_CLOSE(IChOutES_MPI, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " Closing IChOutES"
        GOTO 9999
     ENDIF     
  ENDIF
    
  ! UgMatEffective
  IF(IOutputFLAG.GE.2) THEN
     CALL MPI_FILE_CLOSE(IChOutUM_MPI, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " Closing IChOutUM"
        GOTO 9999
     ENDIF     
  ENDIF

  ALLOCATE( &
       RIndividualReflectionsRoot(2*IPixelCount,&
       2*IPixelCount,IReflectOut,IThicknessCount),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Root Reflections"
     GOTO 9999
  ENDIF
  
  IF(IImageFLAG.GE.2) THEN
     ALLOCATE(&
          CAmplitudeandPhaseRoot(2*IPixelCount,&
          2*IPixelCount,IReflectOut,IThicknessCount),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Root Amplitude and Phase"
        GOTO 9999
     ENDIF
     CAmplitudeandPhaseRoot = CZERO
  END IF

  RIndividualReflectionsRoot = ZERO
  
  IRootArraySize = SIZE(RIndividualReflectionsRoot)
  
  IF(IWriteFLAG.GE.10) THEN
     
     PRINT*,"REDUCING Reflections",my_rank
     
  END IF

  IF(IImageFLAG.LE.1) THEN
     CALL MPI_REDUCE(RIndividualReflections,RIndividualReflectionsRoot,&
          IRootArraySize,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
          MPI_COMM_WORLD,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " In MPI_REDUCE"
        GOTO 9999
     ENDIF   
  ELSE     
     CALL MPI_REDUCE(CAmplitudeandPhase,CAmplitudeandPhaseRoot,&
          IRootArraySize,MPI_DOUBLE_COMPLEX,MPI_SUM,0, &
          MPI_COMM_WORLD,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " In MPI_REDUCE"
        GOTO 9999
     ENDIF   
  END IF

  IF(IWriteFLAG.GE.10) THEN
     PRINT*,"REDUCED Reflections",my_rank
  END IF
  
  IF(IImageFLAG.GE.2) THEN
     DEALLOCATE(&
          CAmplitudeandPhase,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " Deallocating CAmplitudePhase"
        GOTO 9999
     ENDIF   
  END IF
   
  IF(IImageFLAG.LE.1) THEN
     DEALLOCATE( &
          RIndividualReflections,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " Deallocating RIndividualReflections"
        GOTO 9999
     ENDIF   
  END IF

  IF(my_rank.EQ.0) THEN
     ALLOCATE( &
          RFinalMontageImageRoot(MAXVAL(IImageSizeXY),&
          MAXVAL(IImageSizeXY),IThicknessCount),&
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables Root Montage"
        GOTO 9999
     ENDIF
  END IF

  RFinalMontageImageRoot = ZERO

  IF(my_rank.EQ.0.AND.IImageFLAG.GE.2) THEN
     RIndividualReflectionsRoot = &
          CAmplitudeandPhaseRoot * CONJG(CAmplitudeandPhaseRoot)
  END IF

  IF(my_rank.EQ.0) THEN
     DO IThicknessIndex =1,IThicknessCount
        DO ind = 1,2*IPixelCount
           DO jnd = 1,2*IPixelCount
              CALL MakeMontagePixel(ind,jnd,IThicknessIndex,&
                   RFinalMontageImageRoot,&
                   RIndividualReflectionsRoot(ind,jnd,:,IThicknessIndex),IErr)
              IF( IErr.NE.0 ) THEN
                 PRINT*,"main(", my_rank, ") error ", IErr, &
                      " in MakeMontagePixel"
                 GOTO 9999
              ENDIF
           END DO
        END DO
     END DO
  END IF

  !--------------------------------------------------------------------
  ! Write out Images
  !--------------------------------------------------------------------
  
  IF(my_rank.EQ.0) THEN

     ALLOCATE( &
          RImage(2*IPixelCount,2*IPixelCount), &
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables RImage"
        GOTO 9999
     ENDIF
     
     
     IF(IWriteFLAG.GE.10) THEN
        
        PRINT*,"Writing Images"
     END IF
     DO knd = 1,IThicknessCount
        !--------------------------------------------------------
        ! Write Montage
        !--------------------------------------------------------
        
        RThickness = RInitialThickness + (knd-1)*RDeltaThickness 
        IThickness = RInitialThickness + (knd-1)*RDeltaThickness 
        
        WRITE(CThickness,*) IThickness
        WRITE(CThicknessLength,*) SCAN(CThickness,'0123456789',.TRUE.)-SCAN(CThickness,'0123456789')+1
        
        
        IF(IImageFLAG.GE.0) THEN
           WRITE(surname,"(A2,A1,I5.5,A2,I5.5)") &
                "M-","T",IThickness,"-P",MAXVAL(IImageSizeXY)
           CALL OpenReflectionImage(MontageOut,surname,IErr,0,MAXVAL(IImageSizeXY)) 
           IF( IErr.NE.0 ) THEN
              PRINT*,"main(", my_rank, ") error in OpenData()"
              GOTO 9999
           ENDIF
           IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
              
              PRINT*,"main(", my_rank, ") working on RThickness=", RThickness
              
           END IF

           CALL WriteReflectionImage(MontageOut,RFinalMontageImageRoot(:,:,knd), &
                IErr,MAXVAL(IImageSizeXY),MAXVAL(IImageSizeXY))
!!$           DO ind = 1,MAXVAL(IImageSizeXY)
!!$              WRITE(MontageOut,*) &
!!$                   RFinalMontageImageRoot(ind,:,knd) 
!!$           END DO
           CLOSE(MontageOut,ERR=9999)
           
        END IF
        !--------------------------------------------------------
        ! Write Reflections
        !--------------------------------------------------------
        
        IF(IImageFLAG.GE.1) THEN
          
           WRITE(path,"(A2,A1,I1.1,A2,I1.1,A2,I1.1,A2,I4.4,A2,I5.5)") &
                "F-",&
                "S", IScatterFactorMethodFLAG, &
                "_B", ICentralBeamFLAG, &
                "_M", IMaskFLAG, &
                "_P", IPixelCount, &
                "_T", IThickness
           
           call system('mkdir ' // path)
           
           DO ind = 1,IReflectOut
              CALL OpenReflectionImage(IChOutWIImage,path,IErr,ind,2*IPixelCount)
              IF( IErr.NE.0 ) THEN
                 PRINT*,"main(", my_rank, ") error in OpenReflectionImage()"
                 GOTO 9999
              ENDIF

              IF(IImageFLAG.GE.2) THEN
                 
                 CALL OpenReflectionImage(IChOutWFImageReal,path,IErr,ind,2*IPixelCount)
                 IF( IErr.NE.0 ) THEN
                    PRINT*,"main(", my_rank, ") error in OpenAmplitudeImage()"
                    GOTO 9999
                 ENDIF
                 
                 CALL OpenReflectionImage(IChOutWFImagePhase,path,IErr,ind,2*IPixelCount)
                 IF( IErr.NE.0 ) THEN
                    PRINT*,"main(", my_rank, ") error in OpenPhaseImage()"
                    GOTO 9999
                 ENDIF

              END IF

              !-----------------------------------------------------------------------------
              ! Create An Image
              !-----------------------------------------------------------------------------
              IF(IImageFLAG.GE.2) THEN

                 CALL WriteReflectionImage(IChOutWFImageReal,&
                      REAL(CAmplitudeandPhaseRoot(:,:,ind,knd)),IErr,2*IPixelCount,2*IPixelCount)       
                 IF( IErr.NE.0 ) THEN
                    PRINT*,"main(", my_rank, ") error in WriteReflectionImage()"
                    GOTO 9999
                 ENDIF
                 
                 CALL WriteReflectionImage(IChOutWFImagePhase,&
                      AIMAG(CAmplitudeandPhaseRoot(:,:,ind,knd)),IErr,2*IPixelCount,2*IPixelCount)       
                 IF( IErr.NE.0 ) THEN
                    PRINT*,"main(", my_rank, ") error in WriteReflectionImage()"
                    GOTO 9999
                 ENDIF

              END IF
              
              CALL WriteReflectionImage(IChOutWIImage,&
                   RIndividualReflectionsRoot(:,:,ind,knd),IErr,2*IPixelCount,2*IPixelCount)       
              IF( IErr.NE.0 ) THEN
                 PRINT*,"main(", my_rank, ") error in WriteReflectionImage()"
                 GOTO 9999
              ENDIF
              
              IF(IImageFLAG.GE.2) THEN
                 CLOSE(IChOutWFImageReal,ERR=9999)
                 CLOSE(IChOutWFImagePhase,ERR=9999)
              END IF

              CLOSE(IChOutWIImage,ERR=9999)
           END DO
        END IF
     END DO

     DEALLOCATE( &
          RImage,STAT=IErr)       
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error in Deallocation of RImage"
        GOTO 9999
     ENDIF
     

  END IF
  
  !--------------------------------------------------------------------
  ! free memory
  !--------------------------------------------------------------------
  
  !Dellocate Global Variables
  
  DEALLOCATE( &
       RgVecMatT, &
       Rhklpositions, RMask,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in Deallocation of RgVecMatT etc"
     GOTO 9999
  ENDIF
  DEALLOCATE( &
       CUgMat,IPixelLocations, &
       RDevPara,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in Deallocation of CUgmat etc"
     GOTO 9999
  ENDIF
  
  DEALLOCATE( &
       RIndividualReflectionsRoot,STAT=IErr)       
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error in Deallocation of RIndividualReflectionsRoot"
     GOTO 9999
  ENDIF
  
  IF(IImageFLAG.GE.2) THEN
     DEALLOCATE(&
          CAmplitudeandPhaseRoot,STAT=IErr) 
     
     IF( IErr.NE.0 ) THEN
        PRINT*,"main(", my_rank, ") error in Deallocation of CAmplitudeandPhase"
        GOTO 9999
     ENDIF
  END IF
  
  !--------------------------------------------------------------------
  ! finish off
  !--------------------------------------------------------------------
    
  CALL cpu_time(CurrentTime)
  Duration=(CurrentTime-StartTime)
  IHours = FLOOR(Duration/3600.0D0)
  IMinutes = FLOOR(MOD(Duration,3600.0D0)/60.0D0)
  ISeconds = MOD(Duration,3600.0D0)-IMinutes*60.0D0

  PRINT*, "FelixSim(", my_rank, ") ", RStr, ", used time=", IHours, "hrs ",IMinutes,"mins ",ISeconds,"Seconds "

  !--------------------------------------------------------------------
  ! Shut down MPI
  !--------------------------------------------------------------------
9999 &
  CALL MPI_Finalize(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"main(", my_rank, ") error ", IErr, " in MPI_Finalize()"
     STOP
  ENDIF

  ! clean shutdown
  STOP
  
!!$800 PRINT*,"main(", my_rank, "): ERR in CLOSE()"
!!$  IErr= 1
!!$  RETURN

END PROGRAM FelixSim
