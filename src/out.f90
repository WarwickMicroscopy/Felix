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

  CHARACTER*27 :: &
       surname, surnamelength
  CHARACTER*2 :: &
       prefix,postfix
  INTEGER(IKIND) :: &
       IChOutWrite, IErr
  CHARACTER*34 :: &
       filename
  INTEGER(IKIND) :: &
       index

 ! CALL Message("OpenData",IMust,IErr)
  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN

     PRINT*,"OpenData()"

  END IF
  
  WRITE(surnamelength,*) LEN_TRIM(surname)

  WRITE(filename,"(A2,A2,A1,A"//TRIM(surnamelength)//",A4)") "F-",prefix,"-",surname,".txt"

!  CALL Message("OpenData",IAllInfo,IErr,MessageVariable = "filename", & 
!       MessageString = filename)

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
       ACCESS='DIRECT',IOSTAT=Ierr,RECL=2*IPixelCount*8)

  RETURN

  ! error in OPEN detected
10 PRINT*,"OpenImageForReadIn(): ERR in OPEN()"
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

  DO ind=1,2*IPixelCount
     READ(IChInImage,rec=ind) RImageIn(ind,:)
  END DO
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadImageForRefinement (", my_rank, ") error in READ()",IErr
     RETURN
  ENDIF
  
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

  CHARACTER(*) :: &
       surname
  CHARACTER*20 :: &
       prefix,postfix,h,k,l
  INTEGER(IKIND) :: &
       IChOutWrite, IErr,IReflectWriting,IImageSizeX
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
        WRITE(h,*)  NINT(RHKL(IReflectWriting,1))
        WRITE(k,*)  NINT(RHKL(IReflectWriting,2))
        WRITE(l,*)  NINT(RHKL(IReflectWriting,3))
     ELSE
        
        WRITE(h,*)  NINT(RHKL(IOutPutReflections(IReflectWriting),1))
        WRITE(k,*)  NINT(RHKL(IOutPutReflections(IReflectWriting),2))
        WRITE(l,*)  NINT(RHKL(IOutPutReflections(IReflectWriting),3))
     END IF
  END SELECT

  WRITE(Simagesize,"(A2,I5.5,A2,I5.5)") &
       "-P",IImageSizeX,&
       "-P",IImageSizeX

  SELECT CASE (IBinorTextFLAG)
  CASE(0)
     WRITE(fileext,*) TRIM(ADJUSTL(".bin")) 
  CASE(1)
     WRITE(fileext,*) TRIM(ADJUSTL(".txt"))
  END SELECT
  
  SELECT CASE(IChOutWrite)
  CASE(IChOutWFImageReal)        
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"/f-WF-A-hkl-",&
          TRIM(ADJUSTL(h)),&
          TRIM(ADJUSTL(k)),&
          TRIM(ADJUSTL(l)),&
          TRIM(ADJUSTL(Simagesize)),&
          TRIM(ADJUSTL(fileext))
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE FUNCTION REAL PART (WR*.txt)"
     END IF
  CASE(IChOutWFImagePhase)        
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"/f-WF-P-hkl-",&
          TRIM(ADJUSTL(h)),&
          TRIM(ADJUSTL(k)),&
          TRIM(ADJUSTL(l)),&
          TRIM(ADJUSTL(Simagesize)),&
          TRIM(ADJUSTL(fileext))
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE FUNCTION PHASE PART (WP*.txt)"
     END IF
  CASE(IChOutWIImage) 
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"/f-WI-hkl-",&
          TRIM(ADJUSTL(h)),&
          TRIM(ADJUSTL(k)),&
          TRIM(ADJUSTL(l)),&
          TRIM(ADJUSTL(Simagesize)),&
          TRIM(ADJUSTL(fileext))     
     IF (IWriteFLAG.GE.10) THEN
        PRINT*, "OpenImage: opening image for WAVE INTENSITIES"
     END IF
  CASE(MontageOut)        
     WRITE(filename,*) TRIM(ADJUSTL(surname)),"-WI-M",&
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

SUBROUTINE WriteReflectionImage( IChOutWrite, data, IErr,IImageSizeX,IImageSizeY,knd)

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

SUBROUTINE WriteCifFile(IErr)

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

  INTEGER(IKIND) :: &
       IErr
  REAL(RKIND),DIMENSION(SIZE(RAtomSiteFracCoordVec,DIM=1),SIZE(RAtomSiteFracCoordVec,DIM=2)) :: &
       ROutputData
  REAL(RKIND),DIMENSION(2,3) :: &
       RUnitCellParameters
  LOGICAL :: &
       f1

  IF(.NOT.dict_('cif_core.dic','valid')) THEN
     PRINT*,"Requested Core Dictionary not Present"
  END IF

  IF(.NOT.pfile_('felixoutput.cif')) THEN
     PRINT*,"Cif file already exists"
  END IF

  f1 = pdata_('DataBlock') !Open a Data Block

  call close_

END SUBROUTINE WriteCifFile
  
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

  IF(ISoftwareMode.LT.2) THEN

     CALL Message("WriteOutInputFile",IMust,IErr)
     
     OPEN(UNIT= IChInp,FILE= "felix.inp.simdraw_sample",&
       STATUS= 'UNKNOWN')
  
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Input file for felixsim/draw/refine version :VERSION: Build :BUILD:")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# ------------------------------------")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# ------------------------------------")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# felixsim input")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# control flags")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IWriteFLAG                = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IImageFLAG                = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IOutputFLAG               = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IBinorTextFLAG            = 0") 
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IScatterFactorMethodFLAG  = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("ICentralBeamFLAG          = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMaskFLAG                 = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IZolzFLAG                 = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IAbsorbFLAG               = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IAnisoDebyeWallerFlag     = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IBeamConvergenceFLAG      = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IPseudoCubicFLAG          = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionFLAG           = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# radius of the beam in pixels")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IPixelCount               = 64")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# beam selection criteria")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMinReflectionPool        = 600")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMinStrongBeams           = 125")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMinWeakBeams             = 20")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RBSBMax                   = 0.1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RBSPMax                   = 0.1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RConvergenceTolerance (%) = 1.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# crystal settings")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RDebyeWallerConstant      = 0.4668")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RAbsorptionPer            = 2.9")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# microscope settings")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("ROuterConvergenceAngle    = 6.0") 
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RInnerConvergenceAngle    = 0.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IIncidentBeamDirectionX   = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IIncidentBeamDirectionY   = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IIncidentBeamDirectionZ   = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionX              = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionY              = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionZ              = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INormalDirectionX         = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INormalDirectionY         = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INormalDirectionZ         = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RAcceleratingVoltage (kV) = 200.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Image Output Options")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RInitialThickness        = 300.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RFinalThickness          = 1300.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RDeltaThickness          = 10.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IReflectOut              = 49")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
  ELSE
     
     OPEN(UNIT= IChInp,FILE= "felix.inp.refine_sample",&
       STATUS= 'UNKNOWN')
  
  
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Input file for felixsim/draw/refine version :VERSION: Build :BUILD:")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# ------------------------------------")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# ------------------------------------")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# felixsim input")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# control flags")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IWriteFLAG                = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IImageFLAG                = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IOutputFLAG               = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IBinorTextFLAG            = 0") 
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IScatterFactorMethodFLAG  = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("ICentralBeamFLAG          = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMaskFLAG                 = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IZolzFLAG                 = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IAbsorbFLAG               = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IAnisoDebyeWallerFlag     = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IBeamConvergenceFLAG      = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IPseudoCubicFLAG          = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionFLAG           = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# radius of the beam in pixels")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IPixelCount               = 16")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# beam selection criteria")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMinReflectionPool        = 600")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMinStrongBeams           = 200")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IMinWeakBeams             = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RBSBMax                   = 0.1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RBSPMax                   = 0.1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RConvergenceTolerance (%) = 1.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# crystal settings")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RDebyeWallerConstant      = 0.4668")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RAbsorptionPer            = 2.9")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# microscope settings")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("ROuterConvergenceAngle    = 6.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RInnerConvergenceAngle    = 0.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IIncidentBeamDirectionX   = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IIncidentBeamDirectionY   = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IIncidentBeamDirectionZ   = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionX              = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionY              = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IXDirectionZ              = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INormalDirectionX         = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INormalDirectionY         = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INormalDirectionZ         = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RAcceleratingVoltage (kV) = 200.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Image Output Options")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RInitialThickness        = 400.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RFinalThickness          = 700.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RDeltaThickness          = 10.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IReflectOut              = 1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# felixrefine Input")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("#Refinement Specific Flags")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IRefineModeFLAG          = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IWeightingFLAG           = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IContinueFLAG            = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("ICorrelationFLAG         = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IImageProcessingFLAG     = 0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Debye Waller Factor Iteration")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RInitialDebyeWallerFactor = 0.1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RFinalDebyeWallerFactor = 1.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RDeltaDebyeWallerFactor = 0.1")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IElementsforDWFchange = {0}")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Ug Iteration")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("INoofUgs                  = 10")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RLowerBoundUgChange       = 50.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RUpperBoundUgChange       = 50.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RDeltaUgChange            = 50.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Structural Refinement")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IAtomicSites              = (1,2,6)")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Refinement Output")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("IPrint                    = 10")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("# Simplex Initialisation")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("RSimplexLengthScale       = 5.0")
     WRITE(UNIT= IChInp,FMT='(A)') ADJUSTL("")
  CLOSE(UNIT=IChInp)
END IF

END SUBROUTINE WriteOutInputFile
