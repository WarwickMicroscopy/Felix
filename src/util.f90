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
! $Id: util.f90,v 1.95 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Log: util.f90,v $
! Revision 1.95  2014/04/28 12:26:19  phslaz
! Fixed minor bugs with the new reflection pool request
!
! Revision 1.93  2014/04/24 10:48:50  phslaz
! Resorthkl corrected to sort by gmax and not dot(hkl,hkl), language used in resorthkl needs attention, continuous changing of integers is horrible. Codes makes and runs
!
! Revision 1.90  2014/04/09 13:45:39  phslaz
! cleaned up the write flags also added in some of the amplitude/phase imaging
!
! Revision 1.89  2014/03/27 21:39:14  phslaz
! Added two new flags IImageFlag and IOutputFlag
! IImageFLAG = 0 (Montage) 1 (Montage and Reflections)
! IOutputFLAG = 0 (nothing) 1 (EigenSpectra) 2 (UgMat) 3 (Wavefunctions)
! Have also put many Print statments into IWriteflag control
! code compiles and runs
!
! Revision 1.88  2014/03/27 18:13:43  phslaz
! MPI_Reduce attempt, compiles but fails to fun
!
! Revision 1.87  2014/03/26 17:04:52  phslaz
! Felix now creates images
!
! Revision 1.84  2014/03/25 15:45:31  phslaz
! conflict resolution
!
! Revision 1.83  2014/03/25 15:35:34  phsht
! included consistent start of each source file and GPL statements
!
! Revision 1.82  2014/03/25 14:03:39  phslaz
! fixed allocation issue with FelixCoefficients subroutine, allocated memory wasnt being deallocated within loop
!
! Revision 1.81  2014/03/24 12:50:31  phslaz
! fixed write flag 2 and rewrote lacbed to read file once instead of every thickness
!
! Revision 1.80  2014/03/21 15:55:36  phslaz
! New Lacbed code Working
!
! Revision 1.79  2014/03/13 18:10:32  phslaz
! Seg Fault due to 32bit suspected constraint
!
! Revision 1.78  2014/03/07 10:49:45  phslaz
! Corrected issues with inpcif, should now work with badly structured cifs
!
! Revision 1.76  2014/02/21 15:26:42  phslaz
! collapsed a few parts of main.f90 into subroutines in util.f90
!
! Revision 1.75  2014/02/07 14:33:05  phslaz
! LACBED code now reads eigen spectra output
!
! Revision 1.74  2014/01/23 18:50:55  phslaz
! Absorption and 2Pi convention installed and checked with matlab
!
! Revision 1.73  2014/01/20 18:33:59  phslaz
! Isotropic Debye Waller Factor and Atomic Site Partial Occupancy done
!
! Revision 1.71  2014/01/17 15:43:10  phslaz
! Bug Fix : Periodic Table now compiles
!
! Revision 1.69  2014/01/16 16:12:42  phsht
! work on scattering factors
!
! Revision 1.68  2014/01/13 13:46:20  phslaz
! *** empty log message ***
!
! Revision 1.65  2014/01/07 17:11:52  phslaz
! Bug Fix : Hklmake was being called before the basevector assignment and as such, being called with lots of zeros...
!
! Revision 1.64  2014/01/07 11:49:17  phsht
! transformation to microscope reference included, needs testing
!
! Revision 1.63  2013/12/19 16:30:28  phsht
! new version of HKLMake(), using .cif information
!
! Revision 1.62  2013/12/19 14:58:57  phsht
! symmetry operations now correctly interpreted from .cif;
! structure correctly read in
!
! Revision 1.61  2013/12/17 17:40:53  phsht
! make inpcif.f90 which now seems to work
!
! Revision 1.60  2013/11/27 12:30:21  phsht
! more consistency in MPI output routines and their error handling
!
! Revision 1.59  2013/10/21 15:56:31  phslaz
! Changed Variable names
!
! Revision 1.58  2013/10/03 15:48:57  phsht
! added routine HKLMake() in util.f90 to make BWM compatible;
! checked for 64 pixels with input file of this version

! Revision 1.57  2013/10/03 12:52:25  phsht
! two errors corrected for beam selection parameter
!
! Revision 1.56  2013/10/03 11:15:16  phsht
! new input file structure
!
! Revision 1.55  2013/10/02 20:43:28  phsht
! replaced old input names with new names;
! structure of input file still needs changing
!
! Revision 1.54  2013/09/24 16:29:10  phslaz
! Ugmateffective Fixed Now, images now match matlab :)
!
! Revision 1.53  2013/09/24 15:10:07  phslaz
! Added Near Zero Checks
!
! Revision 1.52  2013/09/19 11:19:33  phsht
! eigenvector calculation now working
!
! Revision 1.51  2013/09/18 16:08:42  phsht
! work with Keiths to include changes done with Richard yesterday and
! also to continue until intensities; reached these, but don't quite agree
! yet. Nevertheless, up to the diagonalization, things do work! And eigenvalues
! are also correct, just eigenvectors not yet.
!
! Revision 1.50  2013/09/10 16:51:54  phsht
! before implementation of UgMatEffective
!
! Revision 1.49  2013/09/09 14:09:46  phsht
! works up to the potential
!
! Revision 1.48  2013/09/09 10:58:20  phsht
! subroutines up to this stage seems to produce same output as matlab code
!
! Revision 1.47  2013/09/05 15:20:55  phslaz
! Bug Fix : Didnt declare kstep
!
! Revision 1.46  2013/09/05 15:19:15  phslaz
! Added calculation of kstep
!
! Revision 1.45  2013/09/05 15:08:23  phslaz
! hklpositions is causing issues and isnt massive important right now so ive commented it out for the time being
!
! Revision 1.44  2013/09/05 14:54:46  phslaz
! Bug Fix : Was getting an Unclassifed Function Error from hklpositions, hopefully this will fix it
!
! Revision 1.43  2013/09/05 14:43:39  phslaz
! Added Calculation of positions of reflections in the final image
!
! Revision 1.42  2013/09/05 14:22:48  phslaz
! Beam Selection Works Perfectly
!
! Revision 1.41  2013/09/05 14:21:24  phslaz
! MIssed a comment
!
! Revision 1.40  2013/09/05 14:16:19  phslaz
! Bug Fix : Forgot to uncomment the error
!
! Revision 1.39  2013/09/05 14:04:20  phslaz
! Debugging : Seems to the Sg Calculation Causing the issue, have commented
!
! Revision 1.38  2013/09/05 14:02:56  phslaz
! Debugging : Added Print Lines
!
! Revision 1.37  2013/09/05 11:10:55  phslaz
! Bug Fix : Getting a Segementatiion Fault from beam selection, trying to track it
!
! Revision 1.36  2013/09/05 11:07:40  phslaz
! Bug Fix : forgot an =
!
! Revision 1.35  2013/09/05 11:05:31  phslaz
! Bug Fix i forgot a THEN, i thknk thats the issue
!
! Revision 1.34  2013/09/05 11:01:38  phslaz
! Added Beam Selection Criteria
!
! Revision 1.33  2013/09/05 09:07:33  phslaz
! Resorthkl being called with size(hkl) = 7101 which is the size of the array not the number of reflections (which is 1/3 that)
!
! Revision 1.32  2013/09/05 08:50:38  phslaz
! There was a second call of ReSortHKL but was calling with iErr not Size(HKL) changed it, we'll see what happens
!
! Revision 1.31  2013/09/05 08:46:21  phslaz
! ResortHKL still doing something funny, have added some print lines
!
! Revision 1.30  2013/09/05 08:32:37  phslaz
! Loop wasnt looping over the whole hkl matrix
!
! Revision 1.29  2013/09/05 08:06:10  phsht
! ReSortHKL() now sorts smaller entries first
!
! Revision 1.28  2013/09/05 00:05:03  phslaz
! Resort is backwards, need to invert
!
! Revision 1.27  2013/09/04 23:50:05  phslaz
! Messing around with variables in ILength and ind
!
! Revision 1.26  2013/09/04 23:44:55  phslaz
! Sorting seems to work backwards :)
!
! Revision 1.25  2013/09/04 19:21:46  phsht
! corrected calling of ReSortHKL(); but not sure that result is ok
!
! Revision 1.24  2013/09/04 16:16:06  phslaz
! Changed HKLarray from (N,THREEDIMENSION) to HKLarray(ILength,THREEDIMENSION)
!
! Revision 1.23  2013/09/04 16:12:55  phslaz
! I think ive broke resortHKl
!
! Revision 1.22  2013/09/04 16:09:50  phslaz
! Yeah my previous update broke everything, ive commented out the printstatement
!
! Revision 1.21  2013/09/04 16:07:49  phslaz
! Added Some Prints in, still trying to figure out how the resortHKL works, should it be called with Size(HKL)?
!
! Revision 1.20  2013/09/04 15:38:47  phslaz
! Bug Fix : Trying to diagnose the ResortHKL subroutine, see if theres anything wrong
!
! Revision 1.19  2013/09/04 15:22:30  phslaz
!
! Added a PRINT line in resorthkl trying to bugfix
!
! Revision 1.18  2013/09/04 15:16:57  phslaz
! Removed the definition of gVecMag here
!
! Revision 1.17  2013/09/04 14:59:35  phslaz
! Attempting to Print out HKL after sort
!
! Revision 1.16  2013/09/04 14:52:13  phslaz
! Added A print for minimumG
!
! Revision 1.15  2013/09/04 14:27:16  phslaz
! i wrote int instead of ind, because im an idiot
!
! Revision 1.14  2013/09/04 14:24:18  phslaz
! Modifed a line to output 8 g vector magnitudes not all of them (for debuging)
!
! Revision 1.13  2013/09/03 16:31:33  phslaz
! Bug Fix : Values looked good but i think its running on the unsorted HKL so ive called it prior to the function
!
! Revision 1.12  2013/09/03 16:23:34  phslaz
! Added gVecMag output (to check functionality)
!
! Revision 1.11  2013/09/03 16:20:03  phslaz
! Bug Fix : HopeFully fixed gVecMag calculation
!
! Revision 1.10  2013/09/03 16:10:49  phslaz
! Uncommented gVecMag calculation and added real variable gvecmag(Size(HKL,DIM=1)) (not sure if this will work)
!
! Revision 1.9  2013/09/03 16:05:07  phslaz
! Bug Fix : End - Endif
!
! Revision 1.8  2013/09/03 16:03:05  phslaz
! Added "minimum G Mag"
!
! Revision 1.7  2013/09/02 15:42:42  phsht
! code checked with matlab version up to and including BraggCentral
!
! Revision 1.6  2013/08/29 11:25:15  phsht
! small change
!
! Revision 1.5  2013/06/11 14:53:08  phsht
! more work in the Diffraction defs part
!
! Revision 1.4  2013/06/10 15:08:02  phsht
! work this morning
!
! Revision 1.3  2013/06/10 08:20:28  phsht
! more work before realizing that MATLAB file is not quite the latest version
!
! Revision 1.2  2013/06/07 07:11:28  phsht
! more ground work
!
! Revision 1.1  2013/04/03 19:35:45  phsht
! first installation of basic Fortran routines/structure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! -----------------------------------------------------------------------
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE MicroscopySettings( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  REAL(RKIND) norm,dummy

  INTEGER IErr

  IF(((IWriteFLAG.GE.0).AND.(my_rank.EQ.0)).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"MicroscopySettings(",my_rank,")"
  END IF

  RElectronVelocity= &
       RSpeedOfLight * &
       SQRT( 1.0D0 - ( &
       (RElectronMass*RSpeedOfLight**2) / &
       (RElectronCharge*RAcceleratingVoltage*1.0D3 + &
        RElectronMass*RSpeedOfLight**2) &
       )**2 )

  RElectronWaveLength= RPlanckConstant / &
       ( SQRT(2.0D0*RElectronMass*RElectronCharge*RAcceleratingVoltage*1.D3) * &
         SQRT(1.0D0 + (RElectronCharge*RAcceleratingVoltage*1.D3) / &
          (2.D0*RElectronMass*RSpeedOfLight**2) &
       )) * RAngstromConversion

  RElectronWaveVectorMagnitude=TWOPI/RElectronWaveLength

  RRelativisticCorrection= &
       1.D0 / SQRT( 1.D0 - (RElectronVelocity/RSpeedOfLight)**2 )

  RRelativisticMass= RRelativisticCorrection*RElectronMass

  IF(((IWriteFLAG.GE.1).AND.(my_rank.EQ.0)).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"MicroscopySettings(",my_rank,") ElectronVelocity =", RElectronVelocity
     PRINT*,"MicroscopySettings(",my_rank,") ElectronWaveLength =", RElectronWaveLength
     PRINT*,"MicroscopySettings(",my_rank,") ElectronWaveVectorMagnitude =", RElectronWaveVectorMagnitude
     PRINT*,"MicroscopySettings(",my_rank,") RelativisticCorrection =", RRelativisticCorrection

  END IF
     
  RETURN

END SUBROUTINE MicroscopySettings

! -----------------------------------------------------------------------
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE Crystallography( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
   USE MyMPI

  IMPLICIT NONE

  REAL(RKIND) norm
  REAL(RKIND), DIMENSION(THREEDIM) :: &
       RXDirO, RYDirO, RZDirO, RYDirC
  REAL(RKIND), DIMENSION(THREEDIM,THREEDIM) :: &
       RTMatC2O,RTMatO2M

  INTEGER IErr, ind,jnd,knd, Ifullind, Iuniind
  LOGICAL Lunique
  REAL(RKIND) :: RTTEst

  PRINT*,"DBG: Crystallography()"

  RYDirC = CROSS(RZDirC,RYDirC)

  ! Setup Crystal Lattice Vectors: orthogonal reference frame in Angstrom units

  RaVecO(1)= RLengthX
  RaVecO(2)= 0.D0
  RaVecO(3)= 0.D0

  RbVecO(1)= RLengthY*COS(RGamma)
  RbVecO(2)= RLengthY*SIN(RGamma)
  RbVecO(3)= 0.D0

  RcVecO(1)= RLengthZ*COS(RBeta)
  RcVecO(2)= RLengthZ*(COS(RAlpha)-COS(RBeta)*COS(RGamma))/SIN(RGamma)
  RcVecO(3)= RLengthZ*( &
       SQRT(1.D0- &
        COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)-COS(RGamma)*COS(RGamma) + &
        2.D0*COS(RAlpha)*COS(RBeta)*COS(RGamma)) / &
       SIN(RGamma))

  IF(IVolumeFLAG .EQ. 0) THEN
     RVolume= RLengthX*RLengthY*RLengthZ* &
          SQRT(1.0D0 - &
          COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)-COS(RGamma)*COS(RGamma) + &
          2.0D0*COS(RAlpha)*COS(RBeta)*COS(RGamma))
  ENDIF

  RTTest = DOT_PRODUCT(RaVecO/DOT_PRODUCT(RaVecO,RaVecO),RbVecO/DOT_PRODUCT(RbVecO,RbVecO))*&
       DOT_PRODUCT(RbVecO/DOT_PRODUCT(RbVecO,RbVecO),RcVecO/DOT_PRODUCT(RcVecO,RcVecO))*&
       DOT_PRODUCT(RcVecO/DOT_PRODUCT(RcVecO,RcVecO),RaVecO/DOT_PRODUCT(RaVecO,RaVecO))

  !uncomment this when the hexagonal selection rules have been installed

!!$  IF((RAlpha.GT.PI.OR.RBeta.GT.PI.OR.RGamma.GT.PI).AND.SCAN(SSpaceGroupName,'rR').NE.0) THEN
!!$     SSpaceGroupName = TRIM(ADJUSTL("H"))
!!$     PRINT*,"Crystal is Hexagonal"
!!$  END IF

  IF(ABS(RTTest).LE.TINY.AND.SCAN(SSpaceGroupName,'rR').NE.0) THEN
     SSpaceGroupName = TRIM(ADJUSTL("V"))
     PRINT*,"Crystal is Obverse"
  END IF

  ! Set up Reciprocal Lattice Vectors: orthogonal reference frame in 1/Angstrom units

  RarVecO= TWOPI*CROSS(RbVecO,RcVecO)/DOT(RbVecO,CROSS(RcVecO,RaVecO))
  RbrVecO= TWOPI*CROSS(RcVecO,RaVecO)/DOT(RcVecO,CROSS(RaVecO,RbVecO))
  RcrVecO= TWOPI*CROSS(RaVecO,RbVecO)/DOT(RaVecO,CROSS(RbVecO,RcVecO))

  DO ind=1,THREEDIM
     IF (abs(RarVecO(ind)).lt.1.D-3) THEN
        RarVecO(ind) = 0.0D0
     ENDIF
     IF (abs(RbrVecO(ind)).lt.1.D-3) THEN
        RbrVecO(ind) = 0.0D0
     ENDIF
     IF (abs(RcrVecO(ind)).lt.1.D-3) THEN
        RcrVecO(ind) = 0.0D0
     ENDIF
  ENDDO

  ! Transform from orthogonal reference frame to microscope reference frame
  
  ! RTmat transforms from crystal reference to orthogonal frame
  RTMatC2O(:,1)= RaVecO(:)
  RTMatC2O(:,2)= RbVecO(:)
  RTMatC2O(:,3)= RcVecO(:)

  IF((IWriteFLAG.GE.3.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: RTMatC2O=", RTMatC2O
  END IF

  ! R?DirO vectors are reference vectors in orthogonal frame
  RXDirO= RXDirC(1)*RarVecO + RXDirC(2)*RbrVecO + RXDirC(3)*RcrVecO
  RXDirO= RXDirO/SQRT(DOT_PRODUCT(RXDirO,RXDirO))

  RZDirO= RZDirC(1)*RaVecO + RZDirC(2)*RbVecO + RZDirC(3)*RcVecO
  RZDirO= RZDirO/SQRT(DOT_PRODUCT(RZDirO,RZDirO))

  RYDirO= CROSS(RZDirO,RXDirO)

  ! RTmatO2M transforms from orthogonal to microscope reference frame
  RTMatO2M(1,:)= RXDirO(:)
  RTMatO2M(2,:)= RYDirO(:)
  RTMatO2M(3,:)= RZDirO(:)
  
  IF((IWriteFLAG.GE.3.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: RTMatO2M=", RTMatO2M
  END IF

  ! now transform from crystal reference frame to orthogonal and then to microscope frame

  RXDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RXDirC))
  RYDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RYDirC))
  RZDirM = MATMUL(RTMatO2M,MatMUL(RTMatC2O,RZDirC))

  ! since R?Vec is already in orthogonal frame, only transform into microscope frame needed
  RaVecM= MATMUL(RTMatO2M,RaVecO)
  RbVecM= MATMUL(RTMatO2M,RbVecO)
  RcVecM= MATMUL(RTMatO2M,RcVecO)
  
  ! create new set of reciprocal lattice vectors

  RarVecM= TWOPI*CROSS(RbVecM,RcVecM)/DOT(RbVecM,CROSS(RcVecM,RaVecM))
  RbrVecM= TWOPI*CROSS(RcVecM,RaVecM)/DOT(RcVecM,CROSS(RaVecM,RbVecM))
  RcrVecM= TWOPI*CROSS(RaVecM,RbVecM)/DOT(RaVecM,CROSS(RbVecM,RcVecM))

  ! Calculate the FULL set of possible fractional atomic positions

  ALLOCATE( &
       RFullAtomicFracCoordVec( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1),&
       THREEDIM), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       SFullAtomicNameVec( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       RFullPartialOccupancy( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       RFullIsotropicDebyeWallerFactor( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE( &
       IFullAtomNumber( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       IFullAnisotropicDWFTensor( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  DO ind=1, SIZE(RSymVec,DIM=1)

     DO jnd=1, SIZE(RAtomSiteFracCoordVec,DIM=1)
        Ifullind= SIZE(RSymVec,1)*(jnd-1) + ind

        RFullAtomicFracCoordVec(Ifullind,:)= &
             MATMUL(RSymMat(ind,:,:),RAtomSiteFracCoordVec(jnd,:)) &
             + RSymVec(ind,:)
        SFullAtomicNameVec(Ifullind) = SAtomName(jnd)
        RFullPartialOccupancy(Ifullind) = RAtomicSitePartialOccupancy(jnd)
        RFullIsotropicDebyeWallerFactor = RIsotropicDebyeWallerFactors(jnd)
        IFullAtomNumber(Ifullind) = IAtomNumber(jnd)
        IFullAnisotropicDWFTensor(Ifullind) = IAnisotropicDWFTensor(jnd)

        ! renormalize such that all values are non-negative
        DO knd=1,THREEDIM
           IF( RFullAtomicFracCoordVec(Ifullind,knd) .LT. ZERO) THEN
              RFullAtomicFracCoordVec(Ifullind,knd)= &
                   RFullAtomicFracCoordVec(Ifullind,knd)+1.D0
           ENDIF
        ENDDO

     ENDDO

  ENDDO

  DO ind = 1,SIZE(RFullAtomicFracCoordVec,DIM=1)
     DO jnd = 1,SIZE(RFullAtomicFracCoordVec,DIM=2)
        IF (RFullAtomicFracCoordVec(ind,jnd).LT.ZERO) THEN
           RFullAtomicFracCoordVec(ind,jnd) = &
                RFullAtomicFracCoordVec(ind,jnd) + ONE
        END IF
        IF (RFullAtomicFracCoordVec(ind,jnd).GE.ONE) THEN
           RFullAtomicFracCoordVec(ind,jnd) = &
                RFullAtomicFracCoordVec(ind,jnd) - ONE
        END IF
     END DO
  END DO

  ! Calculate the set of unique fractional atomic positions
  
  ALLOCATE( &
       MNP(ITotalAtoms,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       SMNP(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       RDWF(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       ROcc(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       IAtoms(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       IAnisoDWFT(ITotalAtoms), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Crystallography(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  MNP(1,:)= RFullAtomicFracCoordVec(1,:)
  SMNP(1)= SFullAtomicNameVec(1)
  RDWF(1) = RFullIsotropicDebyeWallerFactor(1)
  ROcc(1) = RFullPartialOccupancy(1)
  IAtoms(1) = IFullAtomNumber(1)
  IAnisoDWFT(1) = IFullAnisotropicDWFTensor(1)
  Iuniind=1

  DO ind=2,SIZE(RFullAtomicFracCoordVec,1)

     DO jnd=1,Iuniind

        IF ( RFullAtomicFracCoordVec(ind,1) .EQ. MNP(jnd,1) .AND. &
             RFullAtomicFracCoordVec(ind,2) .EQ. MNP(jnd,2) .AND. &
             RFullAtomicFracCoordVec(ind,3) .EQ. MNP(jnd,3) .AND. &
             SFullAtomicNameVec(ind) .EQ. SMNP(jnd) ) THEN
           !this seems NOT a unique coordinate
           Lunique=.FALSE.
           EXIT
        ENDIF
        Lunique=.TRUE.
     ENDDO

     IF(Lunique .EQV. .TRUE.) THEN
        Iuniind=Iuniind+1
        MNP(Iuniind,:)= RFullAtomicFracCoordVec(ind,:)
        SMNP(Iuniind)= SFullAtomicNameVec(ind)
!!$        IF(IWriteFLAG.GE.1) THEN
!!$           PRINT*,"MNP(Iuniind,:) = ",MNP(Iuniind,:),SMNP(Iuniind)
!!$        END IF
        RDWF(Iuniind) = RFullIsotropicDebyeWallerFactor(ind)
        ROcc(Iuniind) = RFullPartialOccupancy(ind)
        IAtoms(Iuniind) = IFullAtomNumber(ind)
        IAnisoDWFT(Iuniind) = IFullAnisotropicDWFTensor(ind)
        
     ENDIF
     
  ENDDO

  PRINT*,"Iuniind = ",Iuniind
  
  
  DO ind=1,Iuniind     
     IF((IWriteFLAG.GE.3.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        IF (ind.EQ.1) THEN
           PRINT*,"DBG : Crystallography : ","MNP ","SMNP ","RDWF ","ROcc "
        END IF
        
        PRINT *,"DBG: Crystallography : ",MNP(ind,:),SMNP(ind),RDWF(ind),ROcc(ind),IAtoms(ind),ind
     END IF
  ENDDO

  !IAnisoDWFT now contains a list of indices referring to the correct Anisotropic Debye Wall Factor Tensor for each atom in the unit cell

  DO ind=1,ITotalAtoms
     
     IF((IWriteFLAG.GE.3.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"DBG: IAnisoDWFT",IAnisoDWFT(ind)
     END IF
     
  END DO
  
  ! Calculate atomic position vectors r from Fractional Coordinates and Lattice Vectors
  ! in microscopy reference frame in Angstrom units
  
  DO ind=1,SIZE(MNP,DIM=1)
     DO jnd=1,THREEDIM
        RrVecMat(ind,jnd)= MNP(ind,1)*RaVecM(jnd) + MNP(ind,2)*RbVecM(jnd) + MNP(ind,3)*RcVecM(jnd)
     ENDDO
  ENDDO
  
  INAtomsUnitCell= SIZE(MNP,DIM=1)
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: main(", my_rank, ") NAtomsUnitCell=", INAtomsUnitCell
  END IF
  
  DO ind = 1,INAtomsUnitCell
     DO jnd = 1,THREEDIM
        IF (ABS(RrVecMat(ind,jnd)).LE.TINY) THEN
           RrVecMat(ind,jnd) = ZERO
        ENDIF
     ENDDO
  ENDDO

  ! select the bases in the reciprocal space

  RETURN

END SUBROUTINE Crystallography

! -----------------------------------------------------------------------
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE HKLMake(Ihklmax,Rhkl0Vec,RAcceptanceAngle,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) Ihklmax
  REAL(RKIND) RAcceptanceAngle
  REAL(RKIND), DIMENSION(THREEDIM) :: Rhkl0Vec

  REAL(RKIND), DIMENSION(THREEDIM) :: Rhkl0UnitVec, RhklDummyVec, RhklDummyUnitVec, &
       RuvwVec

  INTEGER IErr, ind,jnd,knd, INhkl
  
  PRINT*,"DBG: HKLMake()"

  SELECT CASE(SSpaceGroupName)
  CASE("F")
     RInvBaseVec= RInvBaseVecF
  CASE("I")
     RInvBaseVec= RInvBaseVecI
  CASE("A")
     RInvBaseVec= RInvBaseVecA
  CASE("B")
     RInvBaseVec= RInvBaseVecB
  CASE("C")
     RInvBaseVec= RInvBaseVecC
  CASE("R")
     RInvBaseVec= RInvBaseVecR
  CASE("V")
     RInvBaseVec= RInvBaseVecV
  CASE("P")
     RInvBaseVec= RInvBaseVecP
  CASE DEFAULT
     PRINT*,"HKLMake(): unknown space group", SSpaceGroupName, "--- aborting"
     IErr=1
     RETURN
  END SELECT

  Rhkl0UnitVec= Rhkl0Vec/SQRT(DOT_PRODUCT(REAL(Rhkl0Vec,RKIND),REAL(Rhkl0Vec,RKIND)))

  INhkl= 0
  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1

           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)

           RhklDummyUnitVec= RhklDummyVec / &
                SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))

           RuvwVec= MATMUL(RInvBaseVec,RhklDummyVec)

           IF( ABS(RuvwVec(1)-NINT(RuvwVec(1))) .LE. TINY .AND. &
                ABS(RuvwVec(2)-NINT(RuvwVec(2))) .LE. TINY .AND. &
                ABS(RuvwVec(3)-NINT(RuvwVec(3))) .LE. TINY ) THEN

              IF(IZolzFLAG.EQ.1) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                   .LE.SIN(RAcceptanceAngle)) THEN
                 INhkl = INhkl +1       
                 
              ENDIF
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  IF(IWriteFLAG.GE.1) THEN
     PRINT*,"DBG: HKLMake, Nhkl=", INhkl
  END IF

  ALLOCATE(RHKL(INhkl+1,THREEDIM), STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"hklmake(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Rhkl"
     RETURN
  ENDIF
  
  IF(IWriteFLAG.GE.3) THEN
     PRINT*,"DBG: BaseVec",RInvBaseVec
  END IF

  INhkl= 0
  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1

           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)

           RhklDummyUnitVec= RhklDummyVec / &
                SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))

           RuvwVec= MATMUL(RInvBaseVec,RhklDummyVec)

           IF( ABS(RuvwVec(1)-NINT(RuvwVec(1))) .LE. TINY .AND. &
                ABS(RuvwVec(2)-NINT(RuvwVec(2))) .LE. TINY .AND. &
                ABS(RuvwVec(3)-NINT(RuvwVec(3))) .LE. TINY ) THEN

              IF(IZolzFLAG.EQ.1) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                    RHKL(INhkl,:)= RhklDummyVec
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                 INhkl =  INhkl + 1
                 RHKL(INhkl,:) = RhklDummyVec                 
              END IF
           END IF
         ENDDO
     ENDDO
  ENDDO

  RHKL(INhkl+1,:)= (/ 0.0D0,0.0D0,0.0D0 /)

END SUBROUTINE HKLMake

! -----------------------------------------------------------------------
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE DiffractionPatternDefinitions( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE SPara
  USE CPara
  USE IChannels
  USE BlochPara
  
  USE MyMPI

  IMPLICIT NONE

  REAL(RKIND) norm, dummyVec(THREEDIM),dummy

  INTEGER IErr, ind,jnd
  
  PRINT*,"DBG: DiffractionPatternDefinitions()"

  CALL NewHKLMake(IHKLMAXValue,RZDirC,2*PI/180.0D0,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error in NewHKLMake()"
     RETURN
  ENDIF

  CALL ReSortHKL( RHKL, SIZE(RHKL,1))
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefintions(): error in ReSortHKL()"
     RETURN
  ENDIF
  
  ALLOCATE(&
       RgVecMat(SIZE(RHKL,DIM=1),THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RgVecMat(HKL)"
     RETURN
  ENDIF

  ALLOCATE(&
       RgVecMatT(SIZE(RHKL,DIM=1),THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RgVecMatT(HKL)"
     RETURN
  ENDIF

  ALLOCATE(&
       RgVecMag(SIZE(RHKL,DIM=1)), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RgVecMag(HKL)"
     RETURN
  ENDIF

  ALLOCATE(&
       RSg(SIZE(RHKL,DIM=1)), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RgVecMag(HKL)"
     RETURN
  ENDIF

  IF(ICentralBeamFLAG.EQ.0) THEN
     !no central beam, 1st beam is used
     dummy= RHKL(1,1)*RHKL(1,1)+RHKL(1,2)*RHKL(1,2)+RHKL(1,3)*RHKL(1,3)
     dummy= SQRT(dummy)
  ELSE
     !there is a central beam, 2nd beam is used
     dummy= RHKL(2,1)*RHKL(2,1)+RHKL(2,2)*RHKL(2,2)+RHKL(2,3)*RHKL(2,3)
     dummy= SQRT(dummy)
  ENDIF

  RBraggCentral= RElectronWaveLength/(2.0D0*RLengthX)*dummy

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: BraggCentral=", RBraggCentral
  END IF

  DO ind=1,SIZE(RHKL,DIM=1)
     DO jnd=1,THREEDIM
        RgVecMatT(ind,jnd)= &
             RHKL(ind,1)*RarVecM(jnd) + &
             RHKL(ind,2)*RbrVecM(jnd) + &
             RHKL(ind,3)*RcrVecM(jnd)
     ENDDO
  ENDDO
  
  ! G vector magnitudes in 1/Angstrom units

  DO ind=1,SIZE(RHKL,DIM=1)
     RgVecMag(ind)= SQRT(DOT(RgVecMatT(ind,:),RgVecMatT(ind,:)))
  ENDDO

  RBSMaxGVecAmp = RgVecMag(IMinReflectionPool)

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DiffractionPatternDefinitions (",my_rank,") GMax = ",RBSMaxGVecAmp
  END IF

  ! smallest g is gmag(2) IF 000 beam is included !!!add error catch here
  
  IF (ICentralBeamFLAG.EQ.1) THEN
     RMinimumGMag = RgVecMag(2)
  ELSE
     RMinimumGMag = RgVecMag(1)
  ENDIF

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: MinimumGMag = ",RMinimumGMag
  END IF

  ! Calculate the Angles to the centre of all the disks from the central disk
  
  DO ind=1,SIZE(RHKL,DIM=1)
     RSg(ind) = 2*RElectronWaveVectorMagnitude* &
          ((-2*RElectronWaveVectorMagnitude +&
          SQRT((2*RElectronWaveVectorMagnitude)**2 + &
          4*RgVecMag(ind)**2))/2)
  ENDDO
  
  ! Determine the number of Gs within GMax (input variable) 
  ! Furthermore, determine which Gs have Sg < SgMax (Strong Beams)
  
  nReflections = 0
  nStrongBeams = 0
  nWeakBeams = 0
  
  DO ind=1,SIZE(RHKL,DIM=1)
     IF (RgVecMag(ind).LE.RBSMaxGVecAmp) THEN
        nReflections = nReflections + 1
        IF (RSg(ind).LE.RBSMaxDeviationPara) THEN
           nStrongBeams = nStrongBeams + 1
        ENDIF
     ENDIF
  ENDDO
  
  IF (nReflections.LT.IReflectOut) THEN

     IReflectOut = nReflections

  END IF

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"DBG: No. of Reflections = ",nReflections

  END IF
  ! resolution in k space
  RDeltaK = RMinimumGMag*RConvergenceAngle/IPixelCount

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: DeltaK = ", RDeltaK
  END IF

  RETURN

END SUBROUTINE DiffractionPatternDefinitions

! -----------------------------------------------------------------------
!
!	IErr	error code
! ----------------------------------------------------------------------

SUBROUTINE ImageInitialization( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  REAL(RKIND) dummyCA

  INTEGER IErr, ind
  
  PRINT*,"DBG: ImageInitialization()"

  !Determine Positions of reflections in final image (may not need to be here)

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: nReflections,MinGMag=", nReflections, RMinimumGMag
  END IF

  ! positions of the centres of the disks

  DO ind=1,nReflections
     Rhklpositions(ind,1) = &
          RgVecMatT(ind,1)/RMinimumGMag
     Rhklpositions(ind,2) = &
          RgVecMatT(ind,2)/RMinimumGMag
  ENDDO
  
  ! size of final image
  
  IF(RConvergenceAngle .LT. ONE) THEN
     dummyCA=RConvergenceAngle
  ELSE
     dummyCA=0.95D0
  ENDIF
  
DO ind=1,SIZE(Rhklpositions,DIM=2)
     IImageSizeXY(ind)= CEILING(&
          4.0D0*REAL(IPixelCount,RKIND)/dummyCA * &
          (MAXVAL(ABS(Rhklpositions(1:IReflectOut,ind)))+1.0D0) )
  ENDDO
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: IImageSizeXY=", IImageSizeXY
  END IF
    
  RETURN

END SUBROUTINE ImageInitialization

!--------------------------------------------------------------------
!	ReSort:
!
!	sort the Lyapunov eigenvalues s.t. the largest comes first. RESORT()
!	is based on ShellSort from "Numerical Recipes", routine SHELL().
!---------------------------------------------------------------------

SUBROUTINE ReSortHKL( RHKLarray, N )

  USE MyNumbers
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  
  IMPLICIT NONE

  INTEGER N
  REAL(RKIND) RHKLarray(N,THREEDIM)
  REAL(RKIND) RhklarraySearch(THREEDIM), RhklarrayCompare(THREEDIM)
  
  REAL(KIND=RKIND) ALN2I, LocalTINY
  PARAMETER (ALN2I=1.4426950D0, LocalTINY=1.D-5)
  
  INTEGER NN,M,L,K,J,I,LOGNB2, index
  REAL(KIND=RKIND) dummy
  
  PRINT*,"DBG: ReSort()"
  
  LOGNB2=INT(LOG(REAL(N))*ALN2I+LocalTINY)
  M=N
  DO 12 NN=1,LOGNB2
     M=M/2
     K=N-M
     DO 11 J=1,K
        I=J
3       CONTINUE
        L=I+M
        RhklarraySearch = Rhklarray(L,1)*RarVecO + &
             Rhklarray(L,2)*RbrVecO + &
             Rhklarray(L,3)*RcrVecO    
        RhklarrayCompare = Rhklarray(I,1)*RarVecO + &
             Rhklarray(I,2)*RbrVecO + &
             Rhklarray(I,3)*RcrVecO
        IF( &
             DOT_PRODUCT(RHKLarraySearch(:),RHKLarraySearch(:)) .LT. &
             DOT_PRODUCT(RHKLarrayCompare(:),RHKLarrayCompare(:))) THEN
!!$           
!!$             DOT_PRODUCT(RHKLarray(L,:),RHKLarray(L,:)) .LT. &
!!$             DOT_PRODUCT(RHKLarray(I,:),RHKLarray(I,:))) THEN
!!$           
           DO 100 index=1,THREEDIM
              dummy        = RHKLarray(I,index)
              RHKLarray(I,index)= RHKLarray(L,index)
              RHKLarray(L,index)= dummy
100        ENDDO
           
           I=I-M
           IF(I.GE.1) GOTO 3
        ENDIF
11   ENDDO
12 ENDDO
  
  !PRINT*,"Finishing ResortHKL"

  !	PRINT*,"array0(1),array0(N)",array0(1),array0(N)
  RETURN

END SUBROUTINE ReSortHKL

!---------------------------------------------------------------------
SUBROUTINE CONVERTAtomName2Number(name, number, IErr)

  USE IPara
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER IErr, ind, number
  CHARACTER*2 name

  INTEGER, PARAMETER :: NElements=103

  CHARACTER*2 A(NElements)

  DATA A/" H", "He", "Li", "Be", " B", " C", " N", "O", "F", "Ne", &
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", &
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", &
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", &
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",& 
        "Md","No","Lr"/

  DO ind=1,NElements
     IF(TRIM(name)==TRIM(A(ind))) THEN
        number= ind
        IF((IWriteFLAG.EQ.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"DBG: name, number ", name, number
        END IF
           RETURN
     ENDIF
  ENDDO

  PRINT*,"CONVERTAtomName2Number(): could not find index for atom of name ", name
  IErr=1
  RETURN

  PRINT*,"DBG: name, number ", name, number
END SUBROUTINE CONVERTAtomName2Number

!---------------------------------------------------------------------
SUBROUTINE ImageMask (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER ind,jnd, ierr
  REAL(RKIND) :: Rradius, RImageRadius
  
  PRINT*,"DBG: ImageMask()"

  IPixelTotal =0
  SELECT CASE (IMaskFLAG)
  CASE(0) ! circle
     DO ind=1,2*IPixelCount
        DO jnd=1,2*IPixelCount
           Rradius= (ind-(REAL(IPixelCount,RKIND)+0.5))**2 + &
                (jnd-(REAL(IPixelCount,RKIND)+0.5))**2
           Rradius=SQRT(DBLE(Rradius))
           RImageRadius = IPixelCount+0.5
           IF(Rradius.LE.RImageRadius) THEN
              RMask(jnd,ind) = 1
              IPixelTotal= IPixelTotal + 1
              
           ELSE
              RMask(jnd,ind) = 0
           ENDIF
        ENDDO
     ENDDO
  CASE(1) ! square
     RMask = 1
     IPixelTotal = (2*IPixelCount)**2
  END SELECT
  
  ALLOCATE( &
       IPixelLocations(IPixelTotal,2), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Imagemask(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  IPixelTotal = 0
 
  SELECT CASE (IMaskFLAG)
  CASE(0) ! circle
     DO ind=1,2*IPixelCount
        DO jnd=1,2*IPixelCount
           Rradius= (ind-(REAL(IPixelCount,RKIND)+0.5))**2 + &
                (jnd-(REAL(IPixelCount,RKIND)+0.5))**2
           Rradius=SQRT(DBLE(Rradius))
           RImageRadius = IPixelCount+0.5
           IF(Rradius.LE.RImageRadius) THEN
              !RMask(jnd,ind) = 1
              IPixelTotal= IPixelTotal + 1
              IPixelLocations(IPixelTotal,1) = ind
              IPixelLocations(IPixelTotal,2) = jnd
           ELSE
              !RMask(jnd,ind) = 0
           ENDIF
        ENDDO
     ENDDO
  CASE(1) ! square
     IPixelTotal = 0
     DO ind = 1,2*IPixelCount
        DO jnd = 1,2*IPixelCount
           IPixelTotal = IPixelTotal+1
           IPixelLocations(IPixelTotal,1) = ind
           IPixelLocations(IPixelTotal,2) = jnd
        END DO
     END DO
  END SELECT
  
END SUBROUTINE ImageMask

!---------------------------------------------------------------------
SUBROUTINE GMatrixInitialisation (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  
  IMPLICIT NONE
  
  INTEGER ind,jnd,ierr

  PRINT*,"DBG: GMatrixInitialisation()"
  
  DO ind=1,nReflections
     DO jnd=1,nReflections
        
        RgMatMat(ind,jnd,:)= RgVecMatT(ind,:)-RgVecMatT(jnd,:)
        RgMatMag(ind,jnd)= SQRT(DOT(RgMatMat(ind,jnd,:),RgMatMat(ind,jnd,:)))
        
     ENDDO
  ENDDO
  
  RgMatMag = RgMatMag/TWOPI
  
END SUBROUTINE GMatrixInitialisation

!---------------------------------------------------------------------
SUBROUTINE UgCalculation (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,ierr, currentatom, iAtom
  COMPLEX(CKIND) CVgij
  REAL(RKIND) RAtomicFormFactor

  PRINT*,"DBG: UgCalculation()"
  
  DO ind=1,nReflections
     DO jnd=1,nReflections
        
        CVgij= 0.0D0
        
        DO iAtom=1, INAtomsUnitCell
           currentatom = IAtoms(iAtom)
           ! calculate f_e(q) as in Eq. (C.15) of Kirkland, "Advanced Computing in EM"
           
           ! for siplicity, we only use Si data from the structure factor data
           ! needs to be made more general later
           RAtomicFormFactor = &
                ! 3 Lorentzians
                RScattFactors(currentatom,1) / &
                (RgMatMag(ind,jnd)**2 + RScattFactors(currentatom,2)) + &
                RScattFactors(currentatom,3) / &
                (RgMatMag(ind,jnd)**2 + RScattFactors(currentatom,4)) + &
                RScattFactors(currentatom,5) / &
                (RgMatMag(ind,jnd)**2 + RScattFactors(currentatom,6)) + &
             ! 3 Gaussians
                RScattFactors(currentatom,7) * &
                EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(currentatom,8)) + &
                RScattFactors(currentatom,9) * &
                EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(currentatom,10)) + &
                RScattFactors(currentatom,11) * &
                EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(currentatom,12))
            
          ! initialize potential as in Eq. (6.10) of Kirkland

           RAtomicFormFactor = RAtomicFormFactor*ROcc(iAtom)

           IF (IAnisoDebyeWallerFactorFlag.EQ.0) THEN
              
              IF(RDWF(iAtom).GT.1.OR.RDWF(iAtom).LT.0) THEN
                 RDWF(iAtom) = RDebyeWallerConstant/(8*PI**2)
              END IF
              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-RgMatMag(ind,jnd)**2*RDWF(iAtom)/3)
              
           ELSE
              
              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-TWOPI*DOT_PRODUCT(RgMatMat(ind,jnd,:), &
                   MATMUL( RAnisotropicDebyeWallerFactorTensor( &
                   IAnisoDWFT(iAtom),:,:), &
                   RgMatMat(ind,jnd,:))))
              
           END IF
           
           CVgij = CVgij + &
                RAtomicFormFactor * &
                EXP(-CIMAGONE* &
                DOT_PRODUCT(RgMatMat(ind,jnd,:), RrVecMat(iAtom,:)) &
                )
        ENDDO
        
        CUgMat(ind,jnd)=(((TWOPI**2)* RRelativisticCorrection) / &
             (PI * RVolume)) * CVgij

        IF (IAbsorbFlag.EQ.1) THEN
           CUgMat(ind,jnd) = &
                CUgMat(ind,jnd) + &
                ABS(CUgMat(ind,jnd))*(RAbsorptionPercentage/100.D0)*CONE       
        END IF
        
     ENDDO
  ENDDO

  RMeanInnerCrystalPotential= REAL(CUgMat(1,1))
  
  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: RMeanInnerCrystalPotential = ",RMeanInnerCrystalPotential
  END IF
  IF(IZolzFLAG.EQ.1) THEN
     DO ind=1,nReflections
        CUgMat(ind,ind)=CUgMat(ind,ind)-RMeanInnerCrystalPotential
     ENDDO
  END IF

END SUBROUTINE UgCalculation

!---------------------------------------------------------------------
SUBROUTINE BlochCoefficientCalculation(ind,jnd,IErr)

  USE MyNumbers
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,hnd,knd,&
       ierr,IThickness, &
       IThicknessIndex, ILowerLimit, &
       IUpperLimit
       
  REAL(RKIND) Rx0,Ry0, RThickness

  CHARACTER*40 surname
  
  COMPLEX(CKIND), DIMENSION(:,:), ALLOCATABLE :: &
       CGeneralSolutionMatrix, CGeneralEigenVectors
  COMPLEX(CKIND),DIMENSION(:),ALLOCATABLE :: &
       CGeneralEigenValues

  Rx0=(ind-IPixelCount-0.5D0)*RDeltaK ! x-position in the disk
  
  Ry0=(jnd-IPixelCount-0.5D0)*RDeltaK ! y-position in the disk
    
  ! we are inside the mask
  IPixelComputed= IPixelComputed + 1

  !--------------------------------------------------------------------
  ! protocol progress
  !--------------------------------------------------------------------
  
  IF(IWriteFLAG.GE.10) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, "): working on pixel (", ind, ",", jnd,") of (", &
          2*IPixelCount, ",", 2*IPixelCount, ") in total."
  ENDIF
       
  !--------------------------------------------------------------------
  ! calculate deviation parameter Sg for the tilted Ewald spheres
  !--------------------------------------------------------------------
  
  ! TiltedK used to be called Kprime2
  ! the vector of the incoming tilted beam

  CALL CalculateKVectors(Rx0,Ry0,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " In Calculation of KVectors"
     RETURN
  ENDIF

  ! Compute the deviation parameter for ALL reflections
  ! within RBSMaxGVecAmp

  CALL DeviationParameterCalculation(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in Calculation of Deviation Parameter"
     RETURN
  ENDIF

  ! select only those beams where the Ewald sphere is close to the
  ! reciprocal lattice, i.e. within RBSMaxDeviationPara

  CALL DetermineStrongAndWeakBeams(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in Determination of Strong and Weak beams"
     RETURN
  ENDIF
 
  ! select the highest reflection that corresponds to a strong beam
  nBeams= IStrongBeamIndex

  !--------------------------------------------------------------------
  ! ALLOCATE memory for eigen problem
  !--------------------------------------------------------------------
  
  !Eigen Problem Solving
  ALLOCATE( &
       CBeamProjectionMatrix(nBeams,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CBeamProjectionMatrix"
     RETURN
  ENDIF
  ALLOCATE( &
       CDummyBeamMatrix(nBeams,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CDummyBeamMatrix"
     RETURN
  ENDIF
  ALLOCATE( &
       CUgMatEffective(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CUgMatEffective"
     RETURN
  ENDIF
  
  !Allocate General Solution Specific Arrays
  
  IF(IZolzFLAG.EQ.0) THEN
     
     ALLOCATE( &
          CGeneralSolutionMatrix(2*nBeams,2*nBeams), & 
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables General Solution Matrix"
        PRINT*,"Failure Occured at Thickness,ChunkPixel,nBeams = ",IPixelCountTotal,nBeams
        RETURN
     ENDIF
     ALLOCATE( &
          CGeneralEigenVectors(2*nBeams,2*nBeams), &
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables CEigenVectors"
        RETURN
     ENDIF
     ALLOCATE(&
          CGeneralEigenValues(2*nBeams), &
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables CEigenVectors"
        RETURN
     ENDIF
  END IF
  
  
  ALLOCATE( & 
       CEigenVectors(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CEigenVectors"
     RETURN
  ENDIF

  ALLOCATE( &
       CEigenValues(nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CEigenValues"
     RETURN
  ENDIF
  ALLOCATE( &
       CInvertedEigenVectors(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CInvertedEigenVectors"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       CAlphaWeightingCoefficients(nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CAlphaWeightingCoefficients"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       CEigenValueDependentTerms(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CEigenValueDependentTerms"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       CWaveFunctions(nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CWaveFunctions"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       RWaveIntensity(nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RWaveIntensity"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       CPsi0(nBeams), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CPsi0"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     RETURN
  ENDIF

  !--------------------------------------------------------------------
  ! construct the effective UgMat (for strong beams only at the moment)
  !--------------------------------------------------------------------
  
  IF(IWriteFLAG.GE.10) THEN 
     PRINT*,"BlochCoefficientCalculation(", my_rank, &
          ") using n(Strong)Beams= ", nBeams, " beams", &
          " with nWeakBeams=", IWeakBeamIndex
  ENDIF
  
  !--------------------------------------------------------------------
  ! back to eigen problem solution
  !--------------------------------------------------------------------
  
  ! compute the effective Ug matrix by selecting only those beams
  ! for which IStrongBeamList has an entry
  
  CBeamProjectionMatrix= CZERO
  DO knd=1,nBeams
     CBeamProjectionMatrix(knd,IStrongBeamList(knd))=CONE
  ENDDO
  
  CUgMatEffective= &
       MATMUL( &
       CBeamProjectionMatrix, &
       MATMUL(CUgMat,TRANSPOSE(CBeamProjectionMatrix)) &
       )
  
  IF (IZolzFLAG.EQ.0) THEN
  
     
     ! General Solution from page 457 of 
     ! "Diffraction of Electrons by Perfect Crystals"
     ! by A.J.F.Metherell
     
     !FILL 0 sub-Maxtrix 
     
     CGeneralSolutionMatrix = CZERO

     CGeneralSolutionMatrix(1:nBeams,(nBeams+1):(nBeams*2)) = CUgMatEffective

     DO hnd = 1,nBeams

        ! FILL I sub-Matrix
        
        CGeneralSolutionMatrix(hnd+nBeams,hnd) = CONE

        ! Fill B sub-Matrix
        CGeneralSolutionMatrix(hnd+nBeams,hnd+nBeams) = &
            -2*RgVecMatT(IStrongBeamList(hnd),3) !2*gz Terms
        
        ! Calculate Beta Values for D sub-Matrix
        
        CGeneralSolutionMatrix(hnd,hnd+nBeams) = &
             (RBigK**2 - & !K^2
             ( &
             (RTiltedK(1))**2 + & !kx^2
             (RTiltedK(2))**2 + & !ky^2
             2*(RTiltedK(1))*RgVecMatT(IStrongBeamList(hnd),1) + & !2*kx*gx
             2*(RTiltedK(2))*RgVecMatT(IStrongBeamList(hnd),2) + & !2*ky*gy
             RgVecMatT(IStrongBeamList(hnd),1)**2 + &  !gx^2
             RgVecMatT(IStrongBeamList(hnd),2)**2 + &  !gx^2
             RgVecMatT(IStrongBeamList(hnd),3)**2 & !gx^2
             ))
     END DO
  ELSE
     
     CUgMatEffective = CUgMatEffective/(TWO*RBigK)
     
      ! set the diagonal parts of the matrix to be equal to 
     ! strong beam deviation parameters (*2 BigK) 
     DO hnd=1,nBeams
        CUgMatEffective(hnd,hnd) = RDevPara(IStrongBeamList(hnd))
     ENDDO
     
  END IF

  !PRINT*,"SIZE of CUgMatEffective = ",SIZE(CUgMatEffective,DIM=1),SIZE(CUgMatEffective,DIM=2)
   
  !--------------------------------------------------------------------
  ! diagonalize the UgMatEffective
  !--------------------------------------------------------------------
   
  IF (IZolzFLAG.EQ.0) THEN
     CALL EigenSpectrum(2*nBeams, &
          CGeneralSolutionMatrix, &
          CGeneralEigenValues(:), CGeneralEigenVectors(:,:), &
          IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in EigenSpectrum()"
        RETURN
     ENDIF

     CEigenValues = CGeneralEigenValues((nBeams+1):(2*nBeams))
     CEigenVectors = CGeneralEigenVectors((nBeams+1):(nBeams*2),1:nBeams)

     DEALLOCATE(&
          CGeneralEigenVectors, &
          CGeneralEigenValues, &
          CGeneralSolutionMatrix)
  ELSE
     CALL EigenSpectrum(nBeams, &
          CUgMatEffective, &
          CEigenValues(:), CEigenVectors(:,:), &
          IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in EigenSpectrum()"
        RETURN
     ENDIF
  END IF
 
  DO IThicknessIndex=1,IThicknessCount,1
     
     RThickness = RInitialThickness + (IThicknessIndex-1)*RDeltaThickness 
     IThickness = RInitialThickness + (IThicknessIndex-1)*RDeltaThickness 
     
     CALL CreateWaveFunctions(rthickness,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in CreateWavefunction()"
        RETURN
     ENDIF

     IF (IOutputFLAG.GE.3) THEN
        IF(IPixelComputed.EQ.1) THEN
           ! wave functions
           CALL OpenData_MPI(IChOutWF_MPI, "WF", surname, IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in OpenDataMPI()"
              RETURN
           ENDIF
        ELSE
           
           ! wave functions
           CALL OpenDataForAppend_MPI(IChOutWF_MPI, "WF", surname, IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in OpenDataForAppend_ MPI()"
              RETURN
           ENDIF
        END IF
        CALL WriteDataC_MPI(IChOutWF_MPI, ind,jnd, &
             CFullWaveFunctions(:), &
             nReflections, 1, IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in WriteDataC_MPI of IChOutWF()"
           RETURN
        ENDIF
        CALL MPI_FILE_CLOSE(IChOutWF_MPI, IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in MPI_FILE_CLOSE of IChOutWF()"
           RETURN
        ENDIF
        
     END IF
     
     !Collection Wave Intensities from all thickness for later writing

     IF(IImageFLAG.LE.1) THEN
        RIndividualReflections(ind,jnd,1:IReflectOut,IThicknessIndex) = &
             RFullWaveIntensity(1:IReflectOut)
     ELSE
        CAmplitudeandPhase(ind,jnd,1:IReflectOut,IThicknessIndex) = &
             CFullWavefunctions(1:IReflectOut)
     END IF

  END DO

  
  
  !--------------------------------------------------------------------
  ! OUTPUT EIGENsystem data for given pixel
  !--------------------------------------------------------------------
  
  IMAXCBuffer = 2*13*SIZE(CEigenVectors)+2*13*SIZE(CEigenValues)+3*SIZE(IStrongBeamList)+3*6*ADD_OUT_INFO
  
  IF(IOutputFLAG.GE.1) THEN
     CALL WriteEigenSystem_MPI(IChOutES_MPI, ind,jnd,nReflections,nBeams, &
          CEigenValues,CEigenVectors, IStrongBeamList,nBeams,nBeams, 1, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in WriteEigenSystem_ MPI()"
        RETURN
     ENDIF
  ENDIF
  
  IMAXCBuffer = 2*14*SIZE(CUgMatEffective)+7*6*ADD_OUT_INFO
  
  ! UgMatEffective
  IF(IOutputFLAG.GE.2) THEN
     CALL WriteDataC_MPI(IChOutUM_MPI, ind,jnd, &
          CUgMatEffective(:,:), nBeams*nBeams, 1, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in WriteDataC_ MPI() of IChOutUM"
        RETURN
     ENDIF
  ENDIF
  
  
  !--------------------------------------------------------------------
  ! DEALLOCATE eigen problem memory
  !--------------------------------------------------------------------
  
  DEALLOCATE( &
       CUgMatEffective,CPsi0,&
       CInvertedEigenVectors, CAlphaWeightingCoefficients, &
       CEigenValues,CEigenVectors,CEigenValueDependentTerms, &
       CBeamProjectionMatrix, CDummyBeamMatrix,CWavefunctions, &
       RWaveIntensity,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in Deallocation()"
     RETURN
  ENDIF
  
END SUBROUTINE BlochCoefficientCalculation

!---------------------------------------------------------------------
SUBROUTINE CountTotalAtoms(IErr)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER ind,jnd,knd,hnd,ierr, ifullind, iuniind
  LOGICAL Lunique

  ALLOCATE( &
       RFullAtomicFracCoordVec( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1),&
       THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE( &
       SFullAtomicNameVec( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  
  RFullAtomicFracCoordVec = ZERO
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"SIZE OF RFULLATOMICFRACCOORDVEC = ",SIZE(RFullAtomicFracCoordVec,1)
  END IF

  DO ind=1, SIZE(RSymVec,DIM=1)
     
     DO jnd=1, SIZE(RAtomSiteFracCoordVec,DIM=1)
       
        Ifullind= SIZE(RSymVec,1)*(jnd-1) + ind
        
        RFullAtomicFracCoordVec(Ifullind,:)= &
             MATMUL(RSymMat(ind,:,:),RAtomSiteFracCoordVec(jnd,:)) &
             + RSymVec(ind,:)
        SFullAtomicNameVec(Ifullind) = SAtomName(jnd)
        
        ! renormalize such that all values are non-negative
        DO knd=1,THREEDIM
           IF( RFullAtomicFracCoordVec(Ifullind,knd) .LT. ZERO) THEN
              RFullAtomicFracCoordVec(Ifullind,knd)= &
                   RFullAtomicFracCoordVec(Ifullind,knd)+1.D0
           ENDIF
        ENDDO
        
     ENDDO
     
  ENDDO

  DO ind = 1,SIZE(RFullAtomicFracCoordVec,DIM=1)
     DO jnd = 1,SIZE(RFullAtomicFracCoordVec,DIM=2)
        IF (RFullAtomicFracCoordVec(ind,jnd).LT.ZERO) THEN
           RFullAtomicFracCoordVec(ind,jnd) = &
                RFullAtomicFracCoordVec(ind,jnd) + ONE
        END IF
        IF (RFullAtomicFracCoordVec(ind,jnd).GE.ONE) THEN
           RFullAtomicFracCoordVec(ind,jnd) = &
                RFullAtomicFracCoordVec(ind,jnd) - ONE
        END IF
     END DO
  END DO

  
  ! Calculate the set of unique fractional atomic positions
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"CountTotalAtoms(",my_rank,") ITotalAtoms = ",ITotalAtoms
  END IF

  ALLOCATE( &
       MNP(1000,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       SMNP(1000), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  
  MNP = ZERO
  
  MNP(1,:)= RFullAtomicFracCoordVec(1,:)
  SMNP(1)= SFullAtomicNameVec(1)
  
  Iuniind = 1
  
  IF(ITotalAtoms.EQ.0)THEN
     
     DO ind=2,SIZE(RFullAtomicFracCoordVec,1)
        
        DO jnd=1,Iuniind
           
           IF ( RFullAtomicFracCoordVec(ind,1) .EQ. MNP(jnd,1) .AND. &
                RFullAtomicFracCoordVec(ind,2) .EQ. MNP(jnd,2) .AND. &
                RFullAtomicFracCoordVec(ind,3) .EQ. MNP(jnd,3) .AND. &
                SFullAtomicNameVec(ind) .EQ. SMNP(jnd) ) THEN
              !this seems NOT a unique coordinate
              Lunique=.FALSE.
              EXIT
           ENDIF
           Lunique=.TRUE.
        ENDDO
        
        IF(Lunique .EQV. .TRUE.) THEN
           Iuniind=Iuniind+1
           MNP(Iuniind,:)= RFullAtomicFracCoordVec(ind,:)
           SMNP(Iuniind)= SFullAtomicNameVec(ind)
        ENDIF
        
     ENDDO
     
     ITotalAtoms = Iuniind
     
  END IF
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"CountTotalAtoms(",my_rank,") ITotalAtoms = ",ITotalAtoms
  END IF

  DEALLOCATE( &
       MNP,SMNP, &
       RFullAtomicFracCoordVec, &
       SFullAtomicNameVec,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
  
END SUBROUTINE CountTotalAtoms

SUBROUTINE CreateWavefunctions(rthickness,IErr)

 USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER ind,jnd,knd,hnd,ierr, ifullind, iuniind,gnd,ichnk
  REAL(RKIND) rthickness 
   
  !--------------------------------------------------------------------
  ! calculate wavefunctions
  !--------------------------------------------------------------------
  
  CPsi0= CZERO
  IF(nBeams .GE. 0) CPsi0(1) = CONE
  
  ! Invert the EigenVector matrix
  CInvertedEigenVectors= CONJG(TRANSPOSE(CEigenVectors(:,:)))
  
  !From EQ 6.32 in Kirkland Advance Computing in EM
  CAlphaWeightingCoefficients = MATMUL(CInvertedEigenVectors(1:nBeams,1:nBeams),CPsi0) 
  
  CEigenValueDependentTerms= ZERO
  
  DO hnd=1,nBeams !IReflectOut 
     
     ! This needs to be a diagonal matrix
     CEigenValueDependentTerms(hnd,hnd) = &
          EXP(CIMAGONE*RThickness*CEigenValues(hnd)) 
     
  ENDDO
  
  ! EQ 6.35 in Kirkland Advance Computing in EM
  ! C-1*C*alpha 
  
  CWaveFunctions(:) = &
       MATMUL( &
       MATMUL(CEigenVectors(1:nBeams,1:nBeams),CEigenValueDependentTerms), & 
       CAlphaWeightingCoefficients(:) &
       )
  
  DO hnd=1,nBeams
     RWaveIntensity(hnd)= &
          CONJG(CWaveFunctions(hnd)) * CWaveFunctions(hnd)
  ENDDO
  
  
  !PRINT*,"This Far"
  
  !--------------------------------------------------------------------
  ! rePADDing of wave function and intensities with zero's 
  !--------------------------------------------------------------------
  
  CFullWaveFunctions=CZERO
  RFullWaveIntensity=ZERO
  
  DO knd=1,nBeams
     CFullWaveFunctions(IStrongBeamList(knd))=CWaveFunctions(knd)
     RFullWaveIntensity(IStrongBeamList(knd))=RWaveIntensity(knd)
  ENDDO
  
END SUBROUTINE CreateWavefunctions

SUBROUTINE MakeMontagePixel(ind,jnd,ithicknessindex,RMontageImage,RIntensity,Ierr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  integer hnd,ind,jnd,knd,Ierr
  INTEGER(IKIND) IThicknessindex, IXpos, IYpos
  REAL(RKIND), DIMENSION(MAXVAL(IImageSizeXY),MAXVAL(IImageSizeXY),IThicknessCount),INTENT(OUT) :: &
       RMontageImage
  REAL(RKIND), DIMENSION(IReflectOut) :: RIntensity
 
  DO hnd = 1,IReflectOut
     
     IF (RConvergenceAngle.LT.ONE) THEN
        IXpos = NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(hnd,2))
        IYpos = NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/RConvergenceAngle)*RhklPositions(hnd,1))
     ELSE
        
        !If the Convergence angle is > 1 causing disk overlap in experimental pattern, then plot as if convergence angle was 0.95 (non-physical but makes a pretty picture)
        IXpos = NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/0.95D0)*RhklPositions(hnd,2))
        IYpos = NINT(MAXVAL(IImageSizeXY)/TWO+TWO*(IPixelCount/0.95D0)*RhklPositions(hnd,1))
     ENDIF
     !PRINT*,"IXpos,IYpos,IPixelCount,jnd,ind = ",IXpos,IYpos,IPixelCount,jnd,ind
     RMontageImage(IXpos-IPixelCount+jnd,IYpos-IPixelCount+ind,IThicknessIndex) = &
          RMontageImage(IXpos-IPixelCount+jnd,IYpos-IPixelCount+ind,IThicknessIndex) + &
          RIntensity(hnd)
  END DO

END SUBROUTINE MakeMontagePixel

SUBROUTINE CalculateKVectors(Rx0,Ry0,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  REAL(RKIND) Rx0,Ry0
  INTEGER(IKIND) :: IErr
  
  RTiltedK(1)= Rx0
  RTiltedK(2)= Ry0
  RTiltedK(3)= SQRT(RBigK**2 - Rx0**2 - Ry0**2)
  
END SUBROUTINE CalculateKVectors

SUBROUTINE DeviationParameterCalculation(IErr)

USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  INTEGER(IKIND) knd,IErr
  
  DO knd=1,nReflections
     ! DevPara used to be called Sg in the book
     
     RDevPara(knd)= &
          -( RBigK + DOT(RgVecMatT(knd,:),RTiltedK(:)) /RBigK) + &
          SQRT( &
          ( RBigK**2 + DOT(RgVecMatT(knd,:),RTiltedK(:)) )**2 /RBigK**2 - &
          (RgVecMag(knd)**2 + &
          2.0D0* DOT(RgVecMatT(knd,:),RTiltedK(:))) &
          )
  END DO

END SUBROUTINE DeviationParameterCalculation

SUBROUTINE DetermineStrongAndWeakBeams(IErr)

  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI
  
  INTEGER(IKIND) ind,knd,IErr,IMinimum,IMaximum,ICheck,jnd
  REAL(RKIND) RDummySg(nReflections)

  !Determine RBSMaxDeviationPara

  RDummySg = ABS(RDevPara)

  DO ind=1,IMinStrongBeams
     IMinimum = MINLOC(RDummySg,1)
     IF(ind.EQ.IMinStrongBeams) THEN
        RBSMaxDeviationPara = ABS(RDummySg(IMinimum))
     ELSE
        RDummySg(IMinimum) = 1000000 !Large number
     END IF
  END DO
  
  IStrongBeamIndex=0
  IWeakBeamIndex=0
  DO knd=1,nReflections
     IF( ABS(RDevPara(knd)) .LE. RBSMaxDeviationPara ) THEN
        IStrongBeamIndex= IStrongBeamIndex +1
        IStrongBeamList(IStrongBeamIndex)= knd
     ENDIF
  ENDDO
  
  RDummySg = ABS(RMeanInnerCrystalPotential/RDevPara)
  
  jnd = 0

  !Determine RBSBethePara

  DO ind=1,nReflections
     ICheck = 0
     IMaximum = MAXLOC(RDummySg,1)

     DO knd = 1,IStrongBeamIndex
        IF(IMaximum.EQ.IStrongBeamList(knd)) THEN
           ICheck = 1
           EXIT
        END IF
     END DO

     IF(ICheck.EQ.0) THEN
        jnd = jnd+1
     END IF

     IF(jnd.EQ.IMinWeakBeams) THEN
        RBSBethePara = (RDummySg(IMaximum))
     ELSE
        RDummySg(IMaximum) = 0.D0 !Large number
     END IF
  END DO

  IWeakBeamIndex=0
  DO knd=1,nReflections
     IF( (ABS(RDevPara(knd)) .GT. RBSMaxDeviationPara).AND. &
          (ABS(RMeanInnerCrystalPotential/RDevPara(knd)) .GE. RBSBethePara) ) THEN
        IWeakBeamIndex= IWeakBeamIndex +1
        IWeakBeamList(IWeakBeamIndex)= knd
     ENDIF
  ENDDO

END SUBROUTINE DetermineStrongAndWeakBeams

SUBROUTINE NewHKLmake(Ihklmax,Rhkl0Vec,RAcceptanceAngle,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) IErr, Ihklmax,ind,jnd,knd,INhkl
  REAL(RKIND) RAcceptanceAngle
  REAL(RKIND), DIMENSION(THREEDIM) :: Rhkl0Vec,RhklDummyUnitVec,RhklDummyVec,Rhkl0UnitVec

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"NewHKLmake(",my_rank,")"
  END IF
  INhkl = 0

  Rhkl0UnitVec= Rhkl0Vec/SQRT(DOT_PRODUCT(REAL(Rhkl0Vec,RKIND),REAL(Rhkl0Vec,RKIND)))

  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1
           
           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)

           RhklDummyUnitVec= RhklDummyVec / &
                SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))
!-NINT(MOD(RhklDummyVec(1)+RhklDummyVec(2),2.D0))
           SELECT CASE(SSpaceGroupName)
           CASE("F") !Face Centred
              IF(((ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),2.D0)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),2.D0)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),2.D0)).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),2.D0))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),2.D0))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),2.D0))).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),2.D0))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),2.D0))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),2.D0))).GT.TINY))) THEN
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                    
                 ENDIF
                 ! INhkl = INhkl + 1
              END IF              
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),2.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),2.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                    
                 ENDIF
              END IF
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),2.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                    
                 ENDIF
              END IF
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),2.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),3.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),3.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
           CASE("P")! Primitive
              !INhkl = INhkl + 1
              
              IF(IZolzFLAG.EQ.1) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                   .LE.SIN(RAcceptanceAngle)) THEN
                 INhkl = INhkl +1       
              ENDIF
           CASE DEFAULT
              PRINT*,"HKLMake(): unknown space group", SSpaceGroupName, "--- aborting"
              IErr=1
           END SELECT

        END DO
     END DO
  END DO
  
  Allocate(&
       RHKL((INhkl+1),THREEDIM),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"hklMake(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables Rhkl"
     RETURN
  ENDIF
  
  INhkl = 0

  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1
           
           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)

           RhklDummyUnitVec= RhklDummyVec / &
                SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))
           
           SELECT CASE(SSpaceGroupName)
           CASE("F") !Face Centred
              IF(((ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),2.D0)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),2.D0)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),2.D0)).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),2.D0))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),2.D0))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),2.D0))).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),2.D0))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),2.D0))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),2.D0))).GT.TINY))) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),2.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),2.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),2.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),2.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),3.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),3.D0)).LE.TINY) THEN
                 !INhkl = INhkl + 1
                 !RHKL(INhkl,:) = RhklDummyVec
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       RHKL(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    RHKL(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
           CASE("P")! Primitive
              !INhkl = INhkl + 1
              !RHKL(INhkl,:) = RhklDummyVec
              IF(IZolzFLAG.EQ.1) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                    RHKL(INhkl,:)= RhklDummyVec
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RAcceptanceAngle)) THEN
                 INhkl =  INhkl + 1
                 RHKL(INhkl,:) = RhklDummyVec                 
              END IF
           CASE DEFAULT
              PRINT*,"HKLMake(): unknown space group", SSpaceGroupName, "--- aborting"
              IErr=1
              RETURN
           END SELECT

        END DO
     END DO
  END DO

  
  RHKL(INhkl+1,:)= (/ 0.0D0,0.0D0,0.0D0 /)

  IF(IWriteFLAG.GE.4) THEN
     DO ind=1,INhkl
        PRINT*,RHKL(ind,:)
     END DO
  END IF

END SUBROUTINE NewHKLmake
