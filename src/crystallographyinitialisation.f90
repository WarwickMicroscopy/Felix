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

SUBROUTINE CrystallographyInitialisation( IErr )

  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr, ind

  !Removed PseudoCubic translation

  CALL ReciprocalLattice(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
             " in ReciprocalLattice"
     RETURN
  ENDIF

  CALL UniqueAtomPositions(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
             " in UniqueAtomPositions "
     RETURN
  ENDIF

  !zz CALL CrystalUniqueFractionalAtomicPostitionsCalculation(IErr)
!zz    IF( IErr.NE.0 ) THEN
!zz      PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
!zz              " in CrystalUniqueFractionalAtomicPostitionsCalculation "
!zz      RETURN
!zz   ENDIF

END SUBROUTINE CrystallographyInitialisation
