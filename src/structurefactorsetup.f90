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
  CHARACTER*200 :: SPrintString
  
  CALL Message("StructureFactorSetup",IMust,IErr)

  !--------------------------------------------------------------------
  ! Calculate Reflection Matrix
  !--------------------------------------------------------------------
  !Allocation--------------------------------------------------------
  ALLOCATE(RgMatrix(nReflections,nReflections,ITHREE),STAT=IErr)
  ALLOCATE(RgMatrixMagnitude(nReflections,nReflections),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(",my_rank,")error in allocations"
     RETURN
  END IF

  !RgMatrix is a matrix of g-vectors that corresponds to the CUgMatNoAbs matrix
  !RgMatrixMagnitude is a matrix of their magnitides
  !RgPool is a list of 2pi*g-vectors in the microscope ref frame, units of 1/A (NB exp(i*q.r), optical convention)
  DO ind=1,nReflections
     DO jnd=1,nReflections
        RgMatrix(ind,jnd,:)= RgPool(ind,:)-RgPool(jnd,:)
        RgMatrixMagnitude(ind,jnd)= SQRT(DOT_PRODUCT(RgMatrix(ind,jnd,:),RgMatrix(ind,jnd,:)))
     ENDDO
  ENDDO
  IF(IWriteFLAG.EQ.3.AND.my_rank.EQ.0) THEN
	PRINT*,"g-vector magnitude matrix (2pi/A)"
	DO ind =1,8
     WRITE(SPrintString,FMT='(16(1X,F5.2))') RgMatrixMagnitude(ind,1:8)
     PRINT*,TRIM(ADJUSTL(SPrintString))
    END DO
  END IF
  
  CALL StructureFactorInitialisation (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(",my_rank,")error in StructureFactorInitialisation"
     RETURN
  END IF
  
  !For symmetry determination, only in Ug refinement
  IF (IRefineMode(1).EQ.1 .OR. IRefineMode(12).EQ.1) THEN
    RgSumMat = SUM(ABS(RgMatrix),3)+RgMatrixMagnitude+ABS(REAL(CUgMatNoAbs))+ABS(AIMAG(CUgMatNoAbs))
  END IF

  !Dellocation-------------------------------------------------------- 
  DEALLOCATE(RgMatrixMagnitude,STAT=IErr)
  DEALLOCATE(RgMatrix,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"StructureFactorSetup(",my_rank,")error deallocations"
     RETURN
  END IF

END SUBROUTINE StructureFactorSetup
