!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
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
!  This file is part of felixsim.
!
!  felixsim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  felixsim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with felixsim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: out.f90,v 1.59 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!---------------------------------------------------------------------
!This file contains all the output subroutines
!---------------------------------------------------------------------



! --------------------------------------------------------------------
! OpenData
! --------------------------------------------------------------------


SUBROUTINE OpenData(IChOutWrite, prefix, surname, IErr)

  USE MyNumbers
  USE WriteToScreen

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  CHARACTER*27 :: surname, surnamelength
  CHARACTER*2 :: prefix,postfix
  INTEGER(IKIND) :: IChOutWrite, IErr
  CHARACTER*34 :: filename
  INTEGER(IKIND) :: index

 ! CALL Message("OpenData",IMust,IErr)
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
!this is redundant
  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(KIND=IKIND) IChOutWrite, IErr

  CHARACTER*34 filename

  IF((IWriteFLAG.GE.10.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"OpenImageForReadIn()"
  END IF
  IF((IWriteFLAG.GE.10.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
    PRINT*,filename
  END IF

  OPEN(UNIT= IChInImage, ERR= 10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),FORM='UNFORMATTED',&
       ACCESS='DIRECT',IOSTAT=Ierr,RECL=2*IPixelCount*8)
  RETURN

  ! error in OPEN detected
10 PRINT*,"OpenImageForReadIn(): ERR in OPEN()"
  IErr= 1
  RETURN
  
END SUBROUTINE OpenImageForReadIn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE ReadImageForRefinement(IErr)
!this is redundant
  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind

  IF((IWriteFLAG.GE.10.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ReadImageForRefinement()"
  END IF

  DO ind=1,2*IPixelCount
     READ(IChInImage,rec=ind,ERR=10) RImageIn(ind,:)
  END DO
  RETURN

10 IErr=1
  PRINT*,"ReadImageForRefinement (", my_rank, ") error in READ()",IErr
  RETURN
  
END SUBROUTINE ReadImageForRefinement

! --------------------------------------------------------------------
! Open Reflection Image
! --------------------------------------------------------------------
SUBROUTINE OpenReflectionImage(IChOutWrite, surname, IErr,IReflectWriting,IImageSizeX,ind)

  USE MyNumbers
  USE WriteToScreen

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  USE IChannels
  IMPLICIT NONE

  CHARACTER(*) :: surname
  CHARACTER*20 :: prefix,postfix,h,k,l
  INTEGER(IKIND) :: IChOutWrite, IErr,IReflectWriting,IImageSizeX
  CHARACTER*250 filename
  CHARACTER*40 fileext
  CHARACTER*60 Simagesize
  INTEGER index,ind

  !!$  Only Prints out this message once when iterating (i.e. when in 1st iteration)

  IF (IMessageCounter.LT.1) THEN
     CALL Message("OpenReflectionImage",IMust,IErr)
      CALL Message("OpenReflectionImage",IMust+IDebug,IErr,&
          MessageString = "is looping. Dependent on ImageFLAG also (called more than once while looping)")
     IMessageCounter = IMessageCounter +1
  END IF

  SELECT CASE(IChOutWrite)
  CASE(MontageOut)
  CASE DEFAULT
     IF(IHKLSelectFLAG.EQ.0) THEn
        WRITE(h,*)  NINT(Rhkl(IReflectWriting,1))
        WRITE(k,*)  NINT(Rhkl(IReflectWriting,2))
        WRITE(l,*)  NINT(Rhkl(IReflectWriting,3))
     ELSE
        WRITE(h,*)  NINT(Rhkl(IOutPutReflections(IReflectWriting),1))
        WRITE(k,*)  NINT(Rhkl(IOutPutReflections(IReflectWriting),2))
        WRITE(l,*)  NINT(Rhkl(IOutPutReflections(IReflectWriting),3))
     END IF
  END SELECT

  WRITE(Simagesize,"(A2,I3.3,A2,I3.3)") &
       "_",IImageSizeX,&
       "x",IImageSizeX

  WRITE(fileext,*) TRIM(ADJUSTL(".bin")) 
  
  SELECT CASE(IChOutWrite)
  
  CASE(IChOutWFImageReal)        
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE FUNCTION REAL PART (WR*.txt)"
     END IF
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"/Real-",&
     TRIM(ADJUSTL(h)),&
     TRIM(ADJUSTL(k)),&
     TRIM(ADJUSTL(l)),&
     TRIM(ADJUSTL(fileext))
	 
  CASE(IChOutWFImagePhase)        
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE FUNCTION PHASE PART (WP*.txt)"
     END IF
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"/Imag-",&
     TRIM(ADJUSTL(h)),&
     TRIM(ADJUSTL(k)),&
     TRIM(ADJUSTL(l)),&
     TRIM(ADJUSTL(fileext))
	 
  CASE(IChOutWIImage) 
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for INTENSITIES"
     END IF
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"/",&
     TRIM(ADJUSTL(h)),&
     TRIM(ADJUSTL(k)),&
     TRIM(ADJUSTL(l)),&
     TRIM(ADJUSTL(fileext))
	 
  CASE(MontageOut)        
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"Montage_",&
          TRIM(ADJUSTL(fileext))
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE INTENSITIES"
     END IF
	 
  CASE DEFAULT
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening UNKNOWN channel ", IChOutWrite
     END IF
	 
  END SELECT
  
  CALL Message("OpenReflectionImage",IInfo,IErr, MessageVariable = "filename", &
       MessageString = filename)


  OPEN(UNIT=IChOutWrite, ERR=10, STATUS= 'UNKNOWN', FILE=TRIM(ADJUSTL(filename)),FORM='UNFORMATTED',&
       ACCESS='DIRECT',IOSTAT=Ierr,RECL=IImageSizeX*8)
  RETURN
   
  ! error in OPEN detected
10 PRINT*,"OpenReflectionImage(): ERR in OPEN()",IErr
  IErr= 1
  RETURN
  
END SUBROUTINE OpenReflectionImage

!-----------------------------------------------------------------
! Write Reflection Images
!-----------------------------------------------------------------

SUBROUTINE WriteReflectionImage( IChOutWrite, data, IErr,IImageSizeX,IImageSizeY,knd)
!this is now redundant
  USE MyNumbers
  USE WriteToScreen

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(KIND=IKIND) IErr,IImageSizeY,IImageSizeX
  REAL(KIND=RKIND),DIMENSION(IImageSizeX,IImageSizeY) :: data
  CHARACTER*100 CSizeofData
  INTEGER ind,knd, IChOutWrite
  CHARACTER*100 SFormatString


  IF (IMessageCounter.LT.2) THEN
     CALL Message("WriteReflectionImage",IMust,IErr)
     CALL Message("WriteReflectionImage",IMust+IDebug,IErr, &
          MessageString = "is looping. Dependent on ImageFLAG also (called more than once while looping)")
     IMessageCounter = IMessageCounter +1
  END IF
     
  DO ind = 1,(IImageSizeY)
     WRITE(IChOutWrite,rec=ind) data(ind,:)
  END DO


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

SUBROUTINE WriteCif(IErr)

  USE MyNumbers
  USE WriteToScreen

  USE IConst
  USE RConst
  
  USE IPara
  USE RPara

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INCLUDE       'ciftbx-f90.cmn'

  INTEGER(IKIND) :: IErr
  REAL(RKIND),DIMENSION(SIZE(RBasisAtomPosition,DIM=1),SIZE(RBasisAtomPosition,DIM=2)) :: &
       ROutputData
  REAL(RKIND),DIMENSION(2,3) :: RUnitCellParameters
  LOGICAL :: f1

  IF(.NOT.dict_('cif_core.dic','valid')) THEN
     PRINT*,"Requested Core Dictionary not Present"
  END IF

  IF(.NOT.pfile_('felixoutput.cif')) THEN
     PRINT*,"Cif file already exists"
  END IF

  f1 = pdata_('DataBlock') !Open a Data Block

  call close_

END SUBROUTINE WriteCif
  
!Write out the sample input file, when none provided
SUBROUTINE WriteOutInputFile (IErr)
  
  USE WriteToScreen
  USE MyNumbers
  
  USE IConst
  USE RConst
  
  USE IPara
  USE RPara
  USE CPara
  USE SPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND):: IErr

!!$  IF(ISoftwareMode.LT.2) THEN
     CALL Message("WriteOutInputFile",IMust,IErr)
     
     OPEN(UNIT= IChInp,FILE= "felix.inp.sample",&
       STATUS= 'UNKNOWN')
     
     CALL WriteToScreenandFile(ADJUSTL("# Input file for felixsim/draw/refine version :VERSION: Build :BUILD:"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# ------------------------------------"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# ------------------------------------"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# control flags"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IWriteFLAG                = 3"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IImageFLAG                = 01"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IScatterFactorMethodFLAG  = 0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IMaskFLAG                 = 1"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IHolzFLAG                 = 0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IAbsorbFLAG               = 1"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IAnisoDebyeWallerFlag     = 0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# radius of the beam in pixels"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IPixelCount               = 64"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# beam selection criteria"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IMinReflectionPool        = 100"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IMinStrongBeams           = 20"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IMinWeakBeams             = 0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# crystal settings"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RDebyeWallerConstant      = 0.4668"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RAbsorptionPer            = 2.9"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# microscope settings"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("ROuterConvergenceAngle    = 6.0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IIncidentBeamDirection    = [0,1,1]"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("IXDirection               = [1,0,0]"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("INormalDirection          = [0,1,1]"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RAcceleratingVoltage (kV) = 200.0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RAcceptanceAngle          = 0.0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("# Image Output Options"),IErr)
     CALL WriteToScreenandFile(ADJUSTL(""),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RInitialThickness        = 400.0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RFinalThickness          = 700.0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("RDeltaThickness          = 10.0"),IErr)
     CALL WriteToScreenandFile(ADJUSTL("INoOfLacbedPatterns              = 7"),IErr)

     IF(ISoftwareMode.EQ.2) THEN
        CALL WriteToScreenandFile(ADJUSTL(""),IErr)
        CALL WriteToScreenandFile(ADJUSTL("#Refinement Specific Flags"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("IRefineModeFLAG          = B"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("IWeightingFLAG           = 0"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("IContinueFLAG            = 0"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("ICorrelationFLAG         = 0"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("IImageProcessingFLAG     = 0"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("RBlurRadius              = 1.45"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("INoofUgs                 = 10"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("IAtomicSites             = (1,2,6)"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("IPrint                   = 10"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("RSimplexLengthScale      = 5.0"),IErr)
        CALL WriteToScreenandFile(ADJUSTL("RExitCriteria            = 0.0001"),IErr)
     END IF
        CLOSE(UNIT=IChInp)
!!$END IF
        
END SUBROUTINE WriteOutInputFile

SUBROUTINE WriteToScreenandFile(SStringtoWrite,IErr)
  !This is a pointless subroutine
  USE WriteToScreen
  USE MyNumbers
  
  USE IConst
  USE RConst
  
  USE IPara
  USE RPara
  USE CPara
  USE SPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND):: IErr
  CHARACTER(*) :: SStringtoWrite

  PRINT*,TRIM(ADJUSTL(SStringtoWrite))
  WRITE(UNIT=IChInp,FMT='(A)') TRIM(ADJUSTL(SStringtoWrite))

 END SUBROUTINE WriteToScreenandFile
