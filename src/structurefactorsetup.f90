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

SUBROUTINE StructureFactorSetup(IErr)

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$%
!!$%     Calculate g-vector matrix (All inter g vectors) and from them
!!$%     the Structure factors which will be available
!!$%
!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  USE WriteToScreen
  USE MyNumbers
  USE IConst

  USE IPara; USE RPara ; USE CPara
  USE BlochPara

  USE MPI
  USE MyMPI


  IMPLICIT NONE

  INTEGER(IKIND) :: ind,jnd,IErr

  CALL Message("StructureFactorSetup",IMust,IErr)

  !--------------------------------------------------------------------
  ! Calculate Reflection Matrix
  !--------------------------------------------------------------------
  !Allocation--------------------------------------------------------
  ALLOCATE(RgMatMat(nReflections,nReflections,THREEDIM),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(",my_rank,")error allocating RgMatMat"
     RETURN
  END IF 
  ALLOCATE(RgMatMag(nReflections,nReflections),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(",my_rank,")error allocating RgMatMag"
     RETURN
  END IF

  !RgMatMat is a matrix of g-vectors that corresponds to the CUgMatNoAbs matrix
  !RgMatMag is a matrix of their magnitides
  !RgPool is a list of g-vectors in the microscope ref frame, units of 1/A, multiplied by 2 pi
  DO ind=1,nReflections
     DO jnd=1,nReflections
        RgMatMat(ind,jnd,:)= RgPoolT(ind,:)-RgPoolT(jnd,:)
        RgMatMag(ind,jnd)= SQRT(DOT_PRODUCT(RgMatMat(ind,jnd,:),RgMatMat(ind,jnd,:)))
     ENDDO
  ENDDO
  !Ug take the 2 pi back out of the magnitude...   
  RgMatMag = RgMatMag/TWOPI

  CALL StructureFactorInitialisation (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(",my_rank,")error in StructureFactorInitialisation"
     RETURN
  END IF
  
  !For symmetry determination, only in Ug refinement
  IF((IRefineModeSelectionArray(1).EQ.1)) THEN
    RgSumMat = SUM(ABS(RgMatMat),3)+RgMatMag+ABS(REAL(CUgMatNoAbs))+ABS(AIMAG(CUgMatNoAbs))
  END IF

  !Dellocation-------------------------------------------------------- 
  DEALLOCATE(RgMatMag,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(",my_rank,")error deallocating RgMatMag"
     RETURN
  END IF
  DEALLOCATE(RgMatMat,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(",my_rank,")error deallocating RgMatMat"
     RETURN
  END IF

END SUBROUTINE StructureFactorSetup
