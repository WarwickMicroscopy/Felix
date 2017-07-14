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

MODULE InputModules
 USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

CONTAINS
!!$Function that reads in a 3D vector from a file
!!$has to have the form [N,N,N] where N is an Integer of any length
!!$Possibly need to make flexible for any type of bracket
 SUBROUTINE ThreeDimVectorReadIn(SUnformattedVector,SOpenBracketDummy, &
       SCloseBracketDummy,RFormattedVector)
  
  IMPLICIT NONE
  
  CHARACTER(*), INTENT(IN) :: &
       SUnformattedVector,SOpenBracketDummy,SCloseBracketDummy
  REAL(RKIND),INTENT(OUT),DIMENSION(3) :: &
       RFormattedVector

  CHARACTER*1 :: &
       SComma=',',SOpenBracket,SCloseBracket
  CHARACTER*100 :: &
       SFormattedVectorX,SFormattedVectorY,SFormattedVectorZ
 
  LOGICAL :: &
       LBACK=.TRUE.
  
  INTEGER(IKIND) :: &
       IOpenBracketPosition, ICloseBracketPosition, IFirstCommaPosition, &
       ILastCommaPosition
  
!!$Trim and adjustL bracket to ensure one character string
  SOpenBracket=TRIM(ADJUSTL(SOpenBracketDummy))
  SCloseBracket=TRIM(ADJUSTL(SCloseBracketDummy))
  
!!$Read in string to for the incident beam direction, convert these
!!$to x,y,z coordinate integers with string manipulation
  
  IOpenBracketPosition=INDEX(SUnformattedVector,SOpenBracket)
  ICloseBracketPosition=INDEX(SUnformattedVector,SCloseBracket)
  IFirstCommaPosition=INDEX(SUnformattedVector,SComma)
  ILastCommaPosition=INDEX(SUnformattedVector,SComma,LBACK)

!!$  Separate the Unformatted Vector String into its three X,Y,Z components
  SFormattedVectorX=TRIM(ADJUSTL(SUnformattedVector(IOpenBracketPosition+1:IFirstCommaPosition-1)))
  SFormattedVectorY=TRIM(ADJUSTL(SUnformattedVector(IFirstCommaPosition+1:ILastCommaPosition-1)))
  SFormattedVectorZ=TRIM(ADJUSTL(SUnformattedVector(ILastCommaPosition+1:ICloseBracketPosition-1)))

!!$Read in each string component into its associated position in real variable RFormattedVector
  READ(SFormattedVectorX,*) RFormattedVector(1)
  READ(SFormattedVectorY,*) RFormattedVector(2)
  READ(SFormattedVectorZ,*) RFormattedVector(3)
  
END SUBROUTINE ThreeDimVectorReadIn
END MODULE InputModules
