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
! $Id: FelixSim.f90,v 1.89 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE MontageSetup(IThicknessIndex,knd,ind,jnd,RFinalMontageImageRoot,RIndividualReflectionsRoot,IErr)
  
  USE MyNumbers
  
  USE MPI
  USE MyMPI
  
  USE IPara; USE RPara
  
  IMPLICIT NONE

  INTEGER(IKIND):: IErr, &
       IThicknessIndex,knd,ind,jnd
  REAL(RKIND),DIMENSION(MAXVAL(IImageSizeXY),&
       MAXVAL(IImageSizeXY),IThicknessCount):: RFinalMontageImageRoot
  REAL(RKIND),DIMENSION(IReflectOut,IThicknessCount,IPixelTotal):: &
       RIndividualReflectionsRoot
  
     DO IThicknessIndex =1,IThicknessCount
        DO knd = 1,IPixelTotal
           ind = IPixelLocations(knd,1)
           jnd = IPixelLocations(knd,2)
           CALL MontageInitialisation(ind,jnd,IThicknessIndex,&
                RFinalMontageImageRoot,&
                RIndividualReflectionsRoot(:,IThicknessIndex,knd),IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"FelixSim(", my_rank, ") error ", IErr, &
                   " in MakeMontagePixel"
              RETURN
           ENDIF
        END DO
     END DO
  

END SUBROUTINE MontageSetup
