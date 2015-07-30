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
! $Id: Inputmodules.f90,v 1.11 2014/03/25 15:37:30 phsht Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE InputModules
  IMPLICIT NONE
CONTAINS
!!$Function that reads in a 3D vector from a file
!!$has to have the form [N,N,N] where N is an Integer of any length
!!$Possibly need to make flexible for any type of bracket
  INTEGER(IKIND), DIMENSION(3) FUNCTION &
       ThreeDimVectorReadIn(SUnformattedVectorString,SOpenBracketDummy, &
       SCloseBracket)
  
  IMPLICIT NONE
  
  CHARACTER* INTENT(IN) :: &
       SUnformattedVector,SOpenBracketDummy,SCloseBracketDummy

  CHARACTER*1 :: &
       SComma=',',SOpenBracket,SCloseBracket
  
  LOGICAL :: &
       LBACK=.TRUE.
  
  INTEGER(IKIND) :: &
       IOpenBracketPosition, ICloseBracketPosition, IFirstCommaPosition, &
       ILastCommaPosition
  
!!$Trim and adjustL bracket to ensure one character string
  SOpenBracket=TRIM(ADJUSTL(SOpenBracketDummy)
  SCloseBracket=TRIM(ADJUSTL(SCloseBracketDummy)
  
!!$Read in string to for the incident beam direction, convert these
!!$to x,y,z coordinate integers with string manipulation
  
  IOpenBracketPosition=INDEX(SIncidentBeamDirection,SOpenBracket)
  ICloseBracketPosition=INDEX(SIncidentBeamDirection,SCloseBracket)
  IFirstCommaPosition=INDEX(SIncidentBeamDirection,SComma)
  ILastCommaPosition=INDEX(SIncidentBeamDirection,SComma,LBACK)
  
  
  READ(ADJUSTL(TRIM(SUnformattedVector(IOpenBracketPosition+1:IFirstCommaPosition-1))),*) &
       ThreeDimVectorReadIn(1)
  READ(ADJUSTL(TRIM(SUnformattedVector(IFirstCommaPosition+1:ILastCommaPosition-1))),*) &
       ThreeDimVectorReadIn(2)
  READ(ADJUSTL(TRIM(SUnformattedVector(ILastCommaPosition+1:ICloseBracketPosition-1))),*) &
       ThreeDimVectorReadIn(3)
  
END FUNCTION 3DVectorReadIn
END MODULE InputModules
