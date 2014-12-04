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

  USE IChannels
  USE MPI
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

  USE IChannels
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

  USE IChannels
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

  USE MPI
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

  USE IChannels

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

  USE IChannels

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

  USE IChannels

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
