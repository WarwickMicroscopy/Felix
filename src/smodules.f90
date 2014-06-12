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
! $Id: smodules.f90,v 1.63 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!$Log: smodules.f90,v $
!Revision 1.63  2014/04/28 12:26:19  phslaz
!Fixed minor bugs with the new reflection pool request
!
!Revision 1.61  2014/04/23 17:18:00  phslaz
!Improved Error checking, all subroutines now include ierr and return to main (in felixsim) or lacbed (in felixdraw) upon ierr.ne.0 and call MPI_FINALISE
!
!Revision 1.60  2014/04/14 16:51:12  phslaz
!Seemingly fixed the rhombahedral problem, turns out theres was a mistake in inpcif where the 3rd angle was being read in incorrectly, have also written a new hklmake which is more understandable and applies selection rules directly rather than mathematically
!
!Revision 1.59  2014/04/09 13:45:39  phslaz
!cleaned up the write flags also added in some of the amplitude/phase imaging
!
!Revision 1.58  2014/03/27 21:39:14  phslaz
!Added two new flags IImageFlag and IOutputFlag
!IImageFLAG = 0 (Montage) 1 (Montage and Reflections)
!IOutputFLAG = 0 (nothing) 1 (EigenSpectra) 2 (UgMat) 3 (Wavefunctions)
!Have also put many Print statments into IWriteflag control
!code compiles and runs
!
!Revision 1.57  2014/03/27 18:13:43  phslaz
!MPI_Reduce attempt, compiles but fails to fun
!
!Revision 1.56  2014/03/26 17:04:52  phslaz
!Felix now creates images
!
!Revision 1.52  2014/03/25 15:45:31  phslaz
!conflict resolution
!
!Revision 1.51  2014/03/25 15:35:34  phsht
!included consistent start of each source file and GPL statements
!
!Revision 1.50  2014/03/21 15:55:36  phslaz
!New Lacbed code Working
!
!Revision 1.48  2014/03/07 10:49:45  phslaz
!Corrected issues with inpcif, should now work with badly structured cifs
!
!Revision 1.46  2014/02/21 15:26:42  phslaz
!collapsed a few parts of main.f90 into subroutines in util.f90
!
!Revision 1.45  2014/02/20 13:17:31  phsht
!removed WF/WI/EX/EV/etc outputs from BER-BLOCH
!combined EX+EV output into ES output and made MPI compatible
!restructured ES file
!
!Revision 1.44  2014/02/20 10:15:23  phslaz
!Working towards improved cif read in, also lacbed now creates montages
!
!Revision 1.43  2014/02/17 14:08:56  phslaz
!Lacbed now works
!
!Revision 1.42  2014/02/07 14:33:05  phslaz
!LACBED code now reads eigen spectra output
!
!Revision 1.41  2014/02/04 15:19:52  phsht
!added more output channels
!
!Revision 1.40  2014/01/20 18:33:59  phslaz
!Isotropic Debye Waller Factor and Atomic Site Partial Occupancy done
!
!Revision 1.39  2014/01/20 15:58:50  phslaz
!Isotropic Debye Waller Factor and Partial Occupancy input from cif
!
!Revision 1.38  2014/01/17 16:57:27  phslaz
!InpCif now reads in isotropic debye waller factors but there are not currently used
!
!Revision 1.37  2014/01/16 16:12:42  phsht
!work on scattering factors
!
!Revision 1.36  2014/01/13 13:46:20  phslaz
!*** empty log message ***
!
!Revision 1.33  2014/01/07 11:49:17  phsht
!transformation to microscope reference included, needs testing
!
!Revision 1.32  2013/12/19 16:30:28  phsht
!new version of HKLMake(), using .cif information
!
!Revision 1.31  2013/12/19 14:58:57  phsht
!symmetry operations now correctly interpreted from .cif;
!structure correctly read in
!
!Revision 1.30  2013/12/17 17:40:53  phsht
!make inpcif.f90 which now seems to work
!
!Revision 1.29  2013/11/29 16:46:20  phsht
!removed thickness loop again, was unwieldly
!
!Revision 1.28  2013/11/27 12:30:21  phsht
!more consistency in MPI output routines and their error handling
!
!Revision 1.27  2013/11/26 22:25:56  phsht
!MPIIO for the rest of the files, but something not right yet
!
!Revision 1.26  2013/11/26 18:27:18  phsht
!MPI file IO now working for WI (intensities)
!
!Revision 1.25  2013/11/25 18:26:33  phsht
!progress on the MPI version
!
!Revision 1.24  2013/11/13 17:48:02  phsht
!1st attempt at LACBED code
!
!Revision 1.23  2013/10/30 13:00:44  phslaz
!Rearranged Variable declarations, allocations and deallocations
!
!Revision 1.22  2013/10/21 15:56:31  phslaz
!Changed Variable names
!
!Revision 1.21  2013/10/03 13:41:41  phsht
!typo in inout.f90 let to wrong thickness values being used;
!now works again
!
!Revision 1.20  2013/10/03 11:15:16  phsht
!new input file structure
!
!Revision 1.19  2013/10/02 20:43:28  phsht
!replaced old input names with new names;
!structure of input file still needs changing
!
!Revision 1.18  2013/09/23 16:52:29  phsht
!OPEN, write and CLOSE of data files now implemented
!
!Revision 1.17  2013/09/10 16:51:54  phsht
!before implementation of UgMatEffective
!
!Revision 1.16  2013/09/09 14:09:46  phsht
!works up to the potential
!
!Revision 1.15  2013/09/09 10:58:33  phsht
!subroutines now works
!
!Revision 1.14  2013/09/06 09:21:01  phslaz
!Added the creation of the circular mask (needed later) made allocatable in smodules, allocated in main
!
!Revision 1.13  2013/09/05 11:01:38  phslaz
!Added Beam Selection Criteria
!
!Revision 1.12  2013/09/04 15:31:15  phslaz
!Bug Fix : I think im making a right hash of this trying to Allocate gVecMag
!
!Revision 1.11  2013/09/04 15:15:45  phslaz
!Defined gVecMag here as opposed to in util
!
!Revision 1.10  2013/09/04 07:22:44  phsht
!added declarations for RBSMaxDeviationPara, BSsgMax, RBSBethePara
!
!Revision 1.9  2013/09/02 15:42:42  phsht
!code checked with matlab version up to and including BraggCentral
!
!Revision 1.8  2013/09/02 15:13:14  phsht
!added two more flags to the input file and inout.f90
!
!Revision 1.7  2013/09/02 13:48:36  phslaz
!"beta" -> "RBeta"
!
!Revision 1.6  2013/06/11 14:53:08  phsht
!more work in the Diffraction defs part
!
!Revision 1.5  2013/06/10 15:08:02  phsht
!work this morning
!
!Revision 1.4  2013/06/10 08:20:28  phsht
!more work before realizing that MATLAB file is not quite the latest version
!
!Revision 1.3  2013/06/07 07:57:32  phsht
!some constants and parameters defined
!
!Revision 1.2  2013/06/07 07:11:28  phsht
!more ground work
!
!Revision 1.1  2013/04/03 19:35:45  phsht
!first installation of basic Fortran routines/structure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!--------------------------------------------------------------------
MODULE CConst
  CHARACTER*18, PARAMETER :: RStr= "$Revision: 1.63 $ "
  CHARACTER*30, PARAMETER :: DStr= "$Date: 2014/04/28 12:26:19 $ "
  CHARACTER*16, PARAMETER :: AStr= "$Author: phslaz $"

  CHARACTER*8 CSpaceGrp(230)

  DATA CSpaceGrp/"P1","P-1","P2","P21","C2","Pm","Pc","Cm",&
       "Cc","P2/m","P21/m","C2/m","P2/c","P21/c","C2/c", &
       "P222","P2221","P21212","P212121","C2221","C222","F222", &
       "I222","I212121","Pmm2","Pmc21","Pcc2","Pma2","Pca21", &
       "Pnc2","Pmm21","Pba2","Pna21","Pnn2","Cmm2","Cmc21","Ccc2", &
       "Amm2","Aem2","Ama2","Aea2","Fmm2","Fdd2","Imm2","Iba2", &
       "Ima2","Pmmm","Pnnn","Pccm","Pban","Pmma","Pnna","Pmna","Pcca", &
       "Pbam","Pccn","Pbcm","Pnnm","Pmmn","Pbcn","Pbca","Pnma","Cmcm", &
       "Cmce","Cmmm","Cccm","Cmme","Ccce","Fmmm","Fddd","Immm","Ibam", &
       "Ibca","Imma","P4","P41","P42","P43","I4","I41","P-4","I-4", &
       "P4/m","P42/m","P4/n","P42/n","I4/m","I41/a","P422","P4212", &
       "P4122","P41212","P4222","P42212","P4322","P43212","I422", &
       "I4122","P4mm","P4bm","P42cm","P42nm","P4cc","P4nc","P42mc", &
       "P42bc","I4mm","I4cm","I41md","I41cd","P-42m","P-42c","P-421m", &
       "P-421c","P-4m2","P-4c2","P-4b2","P-4n2","I-4m2","I-4c2", &
       "I-42m","I-42d","P4/mmm","P4/mcc","P4/nbm","P4/nnc","P4/mbm", &
       "P4/mnc","P4/nmm","P4/ncc","P42/mmc","P42/mcm","P42/nbc","P42/nnm", &
       "P42/mbc","P42/mnm","P42/mnc","P42/ncm","I4/mmm","I4/mcm","I41/amd", &
       "I41/acd","P3","P31","P32","R3","P-3","R-3","P312","P321","P3112", &
       "P3121","P3212","P3221","R32","P3m1","P31m","P3c1","P31c","R3m", &
       "R3c","P-31m","P-31c","P-3m1","P-3c1","R-3m","R-3c","P6","P61", &
       "P65","P62","P64","P63","P-6","P6/m","P63/m","P622","P6122", &
       "P6522","P6222","P6422","P6322","P6mm","P6cc","P63cm","P63mc", &
       "P-6m2","P-6c2","P-62m","P-62c","P6/mmm","P6/mcc","P63/mcm", &
       "P63/mmc","P23","F23","I23","P213","I213","Pm-3","Pn-3","Fm-3", &
       "Fd-3","Im-3","Pa-3","Ia-3","P432","P4232","F432","F4132","I432", &
       "P4332","P4132","I4132","P-43m","F-43m","I-43m","P-43n","F-43c", &
       "I-43d","Pm-3m","Pn-3n","Pm-3n","Pn-3m","Fm-3m","Fm-3c","Fd-3m", &
       "Fd-3c","Im-3m","Ia-3d"/
       

END MODULE CConst

!--------------------------------------------------------------------
MODULE IConst
  USE MyNumbers
  INTEGER(IKIND), PARAMETER :: &
       MAXWriteFLAG= 10, &
       THREEDIM= 3, &
       ADD_OUT_INFO=6, &
       IParallelFLAG=0
END MODULE IConst

!--------------------------------------------------------------------
MODULE RConst
  USE MyNumbers
  
  REAL(RKIND), PARAMETER :: &
       RSpeedOfLight=2.99762458D+8, &
       RElectronMass=9.10938291D-31, &
       RElectronMassMeV=0.510998928, &
       RPlanckConstant=6.62606957D-34, &
       RElectronCharge=1.602176565D-19, &
       RAngstromConversion=1.D10
    
END MODULE RConst

!--------------------------------------------------------------------
MODULE IPara
  USE MyNumbers
  
  !Write Out
  
  INTEGER(IKIND) :: &
       IMAXRBuffer,  IMAXCBuffer     
  
  !Input Flags

  INTEGER(IKIND) :: &
       IWriteFLAG, IScatterFactorMethodFLAG, &
       ICentralBeamFLAG, IMaskFLAG, IVolumeFLAG, &
       IZolzFLAG,IAbsorbFLAG, IAnisoDebyeWallerFactorFlag, &
       IImageFLAG,IOutputFLAG,IBeamConvergenceFLAG,  &
       IPseudoCubicFLAG,IXDirectionFLAG,IBinorTextFLAG

  !Minimum Reflections etc
  INTEGER(IKIND) :: &
       IMinReflectionPool,IMinStrongBeams,IMinWeakBeams

  !OtherFLAGS

  INTEGER(IKIND) :: &
       IDiffractionFLAG=0

  !Disk Radius

  INTEGER(IKIND) :: &
       IPixelCount
  
  !Crystal Settings

  INTEGER(IKIND) :: &
       ITotalAtoms

  ! Name2Atom index
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       IAtomNumber,IAtoms

  !Microscope Settings

  INTEGER(IKIND) :: &
       IIncidentBeamDirectionX, IIncidentBeamDirectionY, &
       IIncidentBeamDirectionZ, &
       IXDirectionX, IXDirectionY, IXDirectionZ, &
       INormalDirectionX,INormalDirectionY,INormalDirectionZ

  !LACBED Input

  INTEGER(IKIND) :: &
       IReflectOut

  !Beams from selection criteria

  INTEGER(IKIND) :: &  
       nReflections,nStrongBeams,nWeakBeams,nBeams,IHKLMAXValue
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       IAnisotropicDWFTensor, IAnisoDWFT

  !Main

  INTEGER(IKIND) :: &
       IPixelTotal, INAtomsUnitCell,IPixelComputed

  INTEGER, DIMENSION(2) :: & 
       IImageSizeXY
  
  !LACBED

  INTEGER(IKIND),DIMENSION(:,:), ALLOCATABLE :: &
       ILACBEDStrongBeamList, IPixelLocation, ISymmetryRelations
  INTEGER(IKIND),DIMENSION(:), ALLOCATABLE :: &
       InBeams,IStrongBeamList

  !inpcif

  INTEGER(IKIND) :: &
      ISymCount
  
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       IFullAtomNumber, IFullAnisotropicDWFTensor

  INTEGER(IKIND) :: &
       IPixelCountTotal

  !LACBED Writing

  INTEGER(IKIND) :: &
       ISeperateFolderFlag
  ! Thickness loop Variables

  INTEGER(IKIND) :: &
       IThicknessCount

  INTEGER(IKIND),DIMENSION(:,:),ALLOCATABLE :: &
       IPixelLocations

  !Refine Parameters

  INTEGER(IKIND), DIMENSION(2) :: &
       IOffset
   
END MODULE IPara

!--------------------------------------------------------------------
MODULE RPara
  USE MyNumbers
  USE RConst
  USE IConst

  !INPUT Section
  
  !Beam Selection Criteria
  
  REAL(RKIND) :: &
       RBSMaxDeviationPara, RBSMaxGVecAmp, RBSBethePara, &
       RConvergenceTolerance,RBSBmax, RBSPMax
  
  !Crystal Settings
  
  REAL(RKIND) :: &
       RLengthX, RLengthY, RLengthZ, RVolume, &
       RAlpha, RBeta, RGamma, &
       RDebyeWallerConstant,RAbsorptionPercentage
  
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: &
       RIsotropicDebyeWallerFactors, RAtomicSitePartialOccupancy, RDWF, ROcc
  
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
       RSymVec,RAtomSiteFracCoordVec, MNP,&
       RUniqueKey
  
  REAL(RKIND), DIMENSION(:,:,:), ALLOCATABLE :: &
       RSymMat

  !Microscope Parameters

  REAL(RKIND) :: &
       RConvergenceAngle, RAcceleratingVoltage

  !LACBED Input

  REAL(RKIND) :: &
       RInitialThickness, &
       RFinalThickness, &
       RDeltaThickness       

  !Debye Waller Factor not sure if we use this 
  REAL(RKIND) :: & 
  !     RGVectorUsePercentage, &
       RMeanSquaredDisplacement


  !HKL indices 
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
       RHKL

  ! scattering factors
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
      RScattFactors 

  ! Microscopy Settings
  REAL(RKIND) :: &
       RElectronVelocity, RElectronWaveLength, &
       RElectronWaveVectorMagnitude, RRelativisticCorrection, &
       RRelativisticMass, RBraggCentral

  ! Crystallography 
  ! Real Space and Reciprocal Lattice Vectors in Orthogonal and Microscope
  ! reference framce
  REAL(RKIND), DIMENSION(THREEDIM) :: &
       RXDirM, RYDirM, RZDirM,& 
       RaVecO, RbVecO, RcVecO, &
       RaVecM, RbVecM, RcVecM, &
       RarVecO, RbrVecO, RcrVecO, &
       RarVecM, RbrVecM, RcrVecM, &
       RXDirC, RZDirC, &
       RNormDirC,RNormDirM
  
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
       RrVecMat

  REAL(RKIND) :: RBaseVec(THREEDIM,THREEDIM), &
       RInvBaseVec(THREEDIM,THREEDIM)
  
  REAL(RKIND), DIMENSION(:,:,:), ALLOCATABLE :: &
       RAnisotropicDebyeWallerFactorTensor
  
  !Diffraction Pattern Definitions
  
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: &
       RgVecMag, RSg
  
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
       RgVecMat, RgVecMatT
  
  REAL(RKIND), DIMENSION(THREEDIM,THREEDIM) :: &
       RTMat

  REAL(RKIND) :: &
       RDeltaK, RMinimumGMag
  
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       RGn

  !Image Initialisation
  
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: & 
       Rhklpositions
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RFinalMontageImage

  !Main Program
  
  REAL(RKIND) :: &
       RMeanInnerCrystalPotential
  
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: & 
       RMask, RgMatMag
  
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: & 
       RgMatMat
  
  !LACBED Program
  
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
       RFullAtomicFracCoordVec
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: &
       RFullPartialOccupancy, RFullIsotropicDebyeWallerFactor

  !WaveFunction Arrays

  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       RWaveIntensity,RFullWaveIntensity
    
  REAL(RKIND), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       RIndividualReflections

  !Refinement Variables

  REAL(RKIND), DIMENSION(:,:),ALLOCATABLE :: &
       RImageIn

  REAL(RKIND) :: &
       RCrossCorrelation

END MODULE RPara

MODULE CPara

  USE MyNumbers

  COMPLEX(CKIND), DIMENSION(:,:), ALLOCATABLE :: &
       CUgMat,CUgMatPrime, CUgMatEffective,CEigenValuesChunk
  COMPLEX(CKIND), DIMENSION(:,:,:), ALLOCATABLE :: &
       CEigenVectorsChunk
  COMPLEX(CKIND),DIMENSION(:),ALLOCATABLE :: &
       CAlphaWeightingCoefficients, CPsi0
  COMPLEX(CKIND),DIMENSION(:,:), ALLOCATABLE :: &
       CEigenValueDependentTerms,CInvertedEigenVectors, &
       CBeamProjectionMatrix,CDummyBeamMatrix
  COMPLEX(CKIND),DIMENSION(:),ALLOCATABLE :: &
       CEigenValues,CGammaValues, CWaveFunctions,CFullWaveFunctions
  COMPLEX(CKIND),DIMENSION(:,:),ALLOCATABLE :: &
       CEigenVectors
  COMPLEX(CKIND), DIMENSION(:,:,:,:), ALLOCATABLE :: &
       CAmplitudeandPhase

END MODULE CPara

!--------------------------------------------------------------------
MODULE SPara
  USE MyNumbers
  
  CHARACTER*1 SSpaceGroupName
  CHARACTER*2, DIMENSION(:), ALLOCATABLE :: &
       SFullAtomicNameVec
  
  CHARACTER*2, DIMENSION(:), ALLOCATABLE :: &
       SAtomName, SMNP
  
END MODULE SPara

!--------------------------------------------------------------------
! Input- and Outputchannels
MODULE IChannels
  INTEGER, PARAMETER :: &
       IChInp= 40, &
       IChOutWF= 41, IChOutWI= 42, &
       IChOutEV= 43, IChOutEX= 44, &
       IChOutUM= 45, IChOut=46, &
       IChInImage = 51
  INTEGER :: &
       IChOutWF_MPI, IChOutWI_MPI, &
       IChOutES_MPI, IChOutUM_MPI, &
       IChOut_MPI 
  INTEGER, PARAMETER :: &
       IChOutWFImageReal= 47, IChOutWFImagePhase= 48, &
       IChOutWIImage= 49, MontageOut = 50
END MODULE IChannels

MODULE BlochPara
  
  USE MyNumbers
  USE IConst
  USE MPI
  USE MyMPI
  
  !--------------------------------------------------------------------
  ! eigen problem variables
  !--------------------------------------------------------------------
  
  REAL(RKIND) RBigK
  
  INTEGER(IKIND) IStrongBeamIndex, IWeakBeamIndex
  INTEGER(IKIND),DIMENSION(:), ALLOCATABLE :: &
       IWeakBeamList
  REAL(RKIND),DIMENSION(:), ALLOCATABLE :: &
       RDevPara
  REAL(RKIND), DIMENSION(THREEDIM) :: &
       RTiltedK
  REAL(8), DIMENSION(:), ALLOCATABLE :: &
       RROutArray, RIOutArray
  COMPLEX(CKIND),DIMENSION(:,:), ALLOCATABLE :: &
       CEigenSaveTemp
END MODULE BlochPara
