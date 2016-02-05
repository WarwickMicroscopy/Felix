!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all rights reserved
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

SUBROUTINE ImageSetup (IErr)

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$%
!!$%      Determines size of final images and creates image shape from IMaskFLAG
!!$%
!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  USE MyNumbers
  USE WriteToScreen

  USE IPara; USE RPara

  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr

  
  ! Image initialisation 
  CALL Message("ImageSetup",IMust,IErr)

  CALL ImageInitialisation( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"ImageSetup(", my_rank, ") error", IErr, &
	  "in ImageInitialistion()"
     RETURN
  END IF
  
  !--------------------------------------------------------------------
  ! define image masks
  !--------------------------------------------------------------------
      
  ALLOCATE(RMask(2*IPixelCount,2*IPixelCount),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ImageSetup(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variable RMask"
     RETURN
  END IF

  !Calls subroutine that sets up masking image

  CALL ImageMaskInitialisation(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ImageSetup(", my_rank, ") error ", IErr, &
          " in ImageMaskInitialisation"
     RETURN
  END IF

  CALL Message("ImageSetup",IInfo,IErr, &
       MessageVariable = "IPixelTotal", IVariable = IPixelTotal)

END SUBROUTINE ImageSetup
