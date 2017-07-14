!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-17, all rights reserved
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: specimentsetup.f90,v 1.59 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Calls the subroutines which set up the experimental part of felix
SUBROUTINE ExperimentalSetup (IErr)

  USE WriteToScreen
  USE MyNumbers
  USE IConst
  
  USE IPara; USE RPara; USE SPara; USE CPara
 
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr

  !--------------------------------------------------------------------
  ! crystallography settings
  CALL CrystallographyInitialisation( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"ExperimentalSetup(",my_rank,")error",IErr,"in CrystallographyInitialisation()"
     RETURN
  ENDIF

  !--------------------------------------------------------------------
  ! diffraction initialization
  !--------------------------------------------------------------------
  CALL DiffractionPatternInitialisation( IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"ExperimentalSetup(",my_rank,")error",IErr,"in DiffractionPatternInitialisation()"
     RETURN
  ENDIF

END SUBROUTINE ExperimentalSetup


