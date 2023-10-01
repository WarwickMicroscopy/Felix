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
! Time:
! Status:  
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

! All modules & procedures conatained in this file:
! SConst
! IConst
! RConst
! IPara
! RPara
! CPara
! SPara
! IChannels
! BlochPara
! Refinement

!>
!! Module-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
MODULE SConst

! the following three lines are now automatically generated via git tags
!  CHARACTER(50), PARAMETER :: RStr= "Version: 1.2"
!  CHARACTER(50), PARAMETER :: DStr= "Date: 30-08-2022"
!  CHARACTER(50), PARAMETER :: AStr= "Refinements B,C,D,F,H and S working, no HOLZ" 

  CHARACTER(8) SAllSpaceGrp(230)
!NB needs more work here, does not have non-standard settings or modern versions
!with the letter e
  DATA SAllSpaceGrp/"P1","P-1","P2","P21","C2","Pm","Pc","Cm",&
       "Cc","P2/m","P21/m","C2/m","P2/c","P21/c","C2/c", &
       "P222","P2221","P21212","P212121","C2221","C222","F222", &
       "I222","I212121","Pmm2","Pmc21","Pcc2","Pma2","Pca21", &
       "Pnc2","Pmm21","Pba2","Pna21","Pnn2","Cmm2","Cmc21","Ccc2", &
       "Amm2","Aem2","Ama2","Aea2","Fmm2","Fdd2","Imm2","Iba2", &
       "Ima2","Pmmm","Pnnn","Pccm","Pban","Pmma","Pnna","Pmna","Pcca", &
       "Pbam","Pccn","Pbcm","Pnnm","Pmmn","Pbcn","Pbca","Pnma","Cmcm", &
       "Cmca","Cmmm","Cccm","Cmme","Ccca","Fmmm","Fddd","Immm","Ibam", &
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

  CHARACTER(2) :: SElementSymbolMatrix(110)!N.B. Number must equal INElements
  DATA SElementSymbolMatrix/"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", &
       "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", &
       "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", &
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", &
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",& 
        "Md","No","Lr","D","Ja","Jb","Jc","Jd","Je","Jf"/!Note element 'Q', 'Ja'... etc. added to END of list

  CHARACTER(8) :: SAlphabet(26)
  DATA SAlphabet/"Aa","Bb","Cc","Dd","Ee","Ff","Gg","Hh","Ii","Jj","Kk","Ll",&
       "Mm","Nn","Oo","Pp","Qq","Rr","Ss","Tt","Uu","Vv","Ww","Xx","Yy","Zz"/

END MODULE SConst
!--------------------------------------------------------------------

!>
!! Module-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
MODULE IConst
  USE MyNumbers
  INTEGER(IKIND), PARAMETER :: &
       MAXWriteFLAG= 10, &
       ADD_OUT_INFO=6, &
       IParallelFLAG=0,&
       IRandomFLAG = 1, &
       IFixedSeed = 123456787

END MODULE IConst
!--------------------------------------------------------------------


!>
!! Module-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
MODULE RConst
  USE MyNumbers
  
  REAL(RKIND), PARAMETER :: &
       RSpeedOfLight=REAL(2.99762458D+8,RKIND), &! in m/s
       RElectronMass=REAL(9.10938291D-31,RKIND), &!in kg
       RPlanckConstant=REAL(6.62606957D-34,RKIND), &! in kg m^2 /s
       RElectronCharge=REAL(1.602176565D-19,RKIND), &!in C
       RAngstromConversion=REAL(1.D10,RKIND)! So [1A (in m)] * RAngstromConversion = 1
  REAL(RKIND), PARAMETER :: RTolerance =REAL(1E-7,RKIND)
    
END MODULE RConst
!--------------------------------------------------------------------

!>
!! Module-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
MODULE IPara
  USE MyNumbers
  USE IConst 
     
  !Write Out
  INTEGER(IKIND) :: IMAXRBuffer,  IMAXCBuffer     
  !Input Flags
  INTEGER(IKIND) :: IWriteFLAG,IDebugFLAG,IScatterFactorMethodFLAG, &
       IVolumeFLAG,IHolzFLAG,IAbsorbFLAG, &
       IBeamConvergenceFLAG,IDevFLAG, &
       IRefineModeFLAG,IOutputFLAG,IPrint,&
       IWeightingFLAG,IRefineMethodFLAG,ICorrelationFLAG,IImageProcessingFLAG,&
       IByteSize
  !Minimum Reflections etc
  INTEGER(IKIND) :: IMinStrongBeams,IMinWeakBeams
  !Simulation size
  INTEGER(IKIND) :: ISizeX, ISizeY, IPixelCount, INFrames
  !Crystal Settings
  INTEGER(IKIND) :: IMaxPossibleNAtomsUnitCell
  INTEGER(IKIND),DIMENSION(:,:), ALLOCATABLE :: IhklLattice
  !Name2Atom index
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: IBasisAtomicNumber,IAtomicNumber
  !Microscope Settings
  INTEGER(IKIND) :: IIncidentBeamDirectionX, IIncidentBeamDirectionY, &
       IIncidentBeamDirectionZ, &
       IXDirectionX, IXDirectionY, IXDirectionZ, &
       INormalDirectionX,INormalDirectionY,INormalDirectionZ
  !Iterative Ug
  INTEGER(IKIND) :: INoofUgs
  !LACBED Input
  INTEGER(IKIND) :: INoOfHKLsAll, INoOfHKLsFrame
  !Beams from selection criteria
  INTEGER(IKIND) :: INhkl,nStrongBeams,nWeakBeams,nBeams,IHKLMAXValue
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: IBasisAnisoDW,IStrongBeamList, IAnisoDW
  !Main
  INTEGER(IKIND) :: IPixelTotal, INAtomsUnitCell
  INTEGER, DIMENSION(2) :: IImageSizeXY
  !Refinement FLAGS
  INTEGER(IKIND) :: IImageOutputFLAG
  !LACBED
  INTEGER(IKIND),DIMENSION(:,:), ALLOCATABLE :: IPixelLocation, ISymmetryRelations, IgOutList
  INTEGER(IKIND),DIMENSION(:), ALLOCATABLE :: InBeams,IEquivalentUgKey
  !inpcif
  INTEGER(IKIND) :: ISymCount,ISpaceGrp,ILN
  INTEGER(IKIND) :: IPixelCountTotal
  !Thickness loop Variables
  INTEGER(IKIND) :: IThicknessCount
  !Tracking reflections
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: IhklsAll,IhklsFrame,ILiveList,ILACBEDList
  INTEGER(IKIND), DIMENSION(50) :: ILACBEDFlag
  !Ug Calculation
  INTEGER(IKIND) :: ICurrentZ,IPsize
  !MPI pixel tracking
  INTEGER(IKIND) :: ILocalPixelCountMin,ILocalPixelCountMax
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: ILocalPixelOffset,ILocalNPix
END MODULE IPara
!--------------------------------------------------------------------

!>
!! Module-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
MODULE RPara
  USE MyNumbers
  USE RConst
  USE IConst

  !INPUT Section 
  !Crystallography
  REAL(RKIND) :: RCellA,RCellB,RCellC,RVolume,RAlpha,RBeta,RGamma, &
       RDebyeWallerConstant,RAbsorptionPercentage
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: RBasisIsoDW, RBasisOccupancy, RIsoDW,&
         ROccupancy,RLatMag
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RSymVec,RBasisAtomPosition, RBasisAtomDelta,&
       RAtomPosition,RUniqueKey,RgLattice
  REAL(RKIND), DIMENSION(:,:,:), ALLOCATABLE :: RSymMat
  !Microscope Parameters
  REAL(RKIND) :: RConvergenceAngle,RAcceleratingVoltage
  REAL(RKIND) :: RElectronVelocity,RElectronWaveLength, &
       RElectronWaveVectorMagnitude,RRelativisticCorrection, &
       RRelativisticMass,RBraggCentral,RFrameAngle
  !LACBED
  REAL(RKIND) :: RInitialThickness,RFinalThickness,RDeltaThickness, &
       RInitialDebyeWallerFactor,RFinalDebyeWallerFactor,&
       RDeltaDebyeWallerFactor
  !Iterative Ugs
  REAL(RKIND) :: RPercentageUgChange
  !Debye Waller Constant, g-vector magnitude, dummy [s'x s'y] for absorption calc
  REAL(RKIND) :: RCurrentB,RCurrentGMagnitude,RSprimeY,RPScale
  !HKL indices 
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: Rhkl 
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RInputHKLs
  ! scattering factors
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RScattFactors!,RPseudoAtom
  ! Crystallography 
  ! Real Space and Reciprocal Lattice Vectors in Orthogonal and Microscope
  ! reference frames
  REAL(RKIND), DIMENSION(ITHREE) :: RXDirM,RYDirM,RZDirM,& 
       RaVecO, RbVecO, RcVecO, &
       RaVec_0, RbVec_0, RcVec_0, &
       RaVecM, RbVecM, RcVecM, &
       RarVecO, RbrVecO, RcrVecO, &
       RarVecM, RbrVecM, RcrVecM, &
       RXDirC_0, RZDirC_0, RXDirC, RZDirC, &
       RXDirO, RYDirO, RZDirO, RNormDirC,RNormDirM
  REAL(RKIND) :: RarMag, RbrMag, RcrMag
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RAtomCoordinate
  REAL(RKIND) :: RBaseVec(ITHREE,ITHREE), &
       RInvBaseVec(ITHREE,ITHREE)
  REAL(RKIND), DIMENSION(:,:,:), ALLOCATABLE :: RAnisotropicDebyeWallerFactorTensor
  !Diffraction Pattern Definitions
  REAL(RKIND) RBigK, RDevLimit, RMeanInnerPotential, RScattFacToVolts, RDeltaK, RMinimumGMag, RGVectorMagnitude
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: RgPoolMag
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RgPool
  REAL(RKIND), DIMENSION(ITHREE,ITHREE) :: RTMat
  REAL(RKIND),DIMENSION(ITHREE) :: RGVector
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RgDotNorm
  REAL(RKIND),DIMENSION(:), ALLOCATABLE :: RDevPara  ! deviation parameter for each g, for a given pixel
  REAL(RKIND),DIMENSION(:), ALLOCATABLE :: RDevC  ! deviation parameter for each g at the image centre
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: RgMatrixMagnitude, RgSumMat
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RgMatrix, RgPoolList
  !WaveFunction Arrays
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RWaveIntensity,RFullWaveIntensity,RSumIntensity
  REAL(RKIND), DIMENSION(:,:,:), ALLOCATABLE :: RIndividualReflections
  !Refinement Variables
  ! Simulated Images as a 1D array (no. of patterns, no of thicknesses, totalpixels)
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RSimulatedPatterns,RTempImage
  ! Simulated Images as images (width,height, no.of patterns, no of thicknesses)
  REAL(RKIND),DIMENSION(:,:,:,:),ALLOCATABLE :: RImageSimi
  ! Output images - up to 20 per frame (may need increasing in the future!)
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_1
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_2
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_3
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_4
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_5
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_6
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_7
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_8
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_9
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_10
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_11
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_12
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_13
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_14
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_15
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_16
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_17
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_18
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_19
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_20
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_21
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_22
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_23
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_24
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_25
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_26
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_27
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_28
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_29
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_30
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_31
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_32
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_33
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_34
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_35
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_36
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_37
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_38
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_39
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_40
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_41
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_42
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_43
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_44
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_45
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_46
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_47
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_48
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_49
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_50
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_51
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_52
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_53
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_54
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_55
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_56
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_57
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_58
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_59
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_60
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_61
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_62
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_63
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_64
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_65
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_66
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_67
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_68
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_69
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_70
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_71
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_72
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_73
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_74
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_75
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_76
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_77
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_78
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_79
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_80
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_81
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_82
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_83
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_84
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_85
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_86
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_87
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_88
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_89
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_90
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_91
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_92
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_93
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_94
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_95
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_96
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_97
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_98
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_99
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: RLACBED_100
  !Gaussian blur radius in pixels
  REAL(RKIND) :: RBlurRadius
END MODULE RPara
!--------------------------------------------------------------------

!>
!! Module-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
MODULE CPara
  USE MyNumbers

  COMPLEX(CKIND),DIMENSION(:),ALLOCATABLE :: CAlphaWeightingCoefficients, CPsi0,CUniqueUg,CEigenValues,&
                CGammaValues, CWaveFunctions,CFullWaveFunctions
  COMPLEX(CKIND), DIMENSION(:,:), ALLOCATABLE :: CUgMatNoAbs,CUgMatPrime,CUgMat,CUgSgMatrix,CEigenValuesChunk,&
                CEigenVectors,CEigenValueDependentTerms,CInvertedEigenVectors,CBeamProjectionMatrix,&
                CDummyBeamMatrix
  COMPLEX(CKIND), DIMENSION(:,:,:), ALLOCATABLE :: CEigenVectorsChunk,CAmplitudeandPhase,CPseudoAtom,CPseudoScatt

END MODULE CPara
!--------------------------------------------------------------------

!>
!! Module-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
MODULE SPara
  USE MyNumbers
  
  CHARACTER(500) :: SPrintString
  CHARACTER(40) :: SChemicalFormula
  CHARACTER(1) :: SSpaceGroupName
  CHARACTER(10) :: SSpaceGrp
  CHARACTER(30), DIMENSION(:), ALLOCATABLE :: SSymString
  CHARACTER(5), DIMENSION(:), ALLOCATABLE :: SBasisAtomLabel,SAtomLabel
  CHARACTER(2), DIMENSION(:), ALLOCATABLE :: SBasisAtomName, SAtomName
  CHARACTER(1), DIMENSION(:), ALLOCATABLE :: SWyckoffSymbol
  
END MODULE SPara
!--------------------------------------------------------------------

!>
!! Module-description: Input and output channels
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
MODULE IChannels
  INTEGER, PARAMETER :: IChInp=40, IChOutIM=41, IChOutRC=42,IChOutIhkl=43
END MODULE IChannels
!--------------------------------------------------------------------

