!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all rights reserved
!
! Version: $VERSION$
! Date:    $DATE$
! Time:    $TIME$
! Status:  $RLSTATUS$
! Build:   $BUILD$
! Author:  $AUTHOR$
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

SUBROUTINE CrystallographyInitialisation( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND)IErr, ind
 

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"CrystallographyInitialisation(",my_rank,")"
  END IF
    
  IF(IDiffractionFLAG.EQ.1) THEN
     DEALLOCATE(&
          RFullAtomicFracCoordVec,SFullAtomicNameVec,&
          RFullPartialOccupancy,RFullIsotropicDebyeWallerFactor, &
          IFullAtomNumber, IFullAnisotropicDWFTensor, &
          MNP,SMNP,RDWF,ROcc,IAtoms,IAnisoDWFT,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
             " in DEALLOCATE()"
        RETURN
     ENDIF
  END IF
  
  IF(IPseudoCubicFLAG.EQ.1.AND.IDiffractionFLAG.EQ.0) THEN
     RZDirC(1) = (IIncidentBeamDirectionX/3.0D0)+&
          (IIncidentBeamDirectionY/3.0D0)-&
          (2.0D0*(IIncidentBeamDirectionZ/3.0D0))
     RZDirC(2) = (-IIncidentBeamDirectionX/3.0D0)+&
          (2.0D0*(IIncidentBeamDirectionY/3.0D0))-&
          (IIncidentBeamDirectionZ/3.0D0)
     RZDirC(3) = (IIncidentBeamDirectionX/3.0D0)+&
          (IIncidentBeamDirectionY/3.0D0)+&
          (IIncidentBeamDirectionZ/3.0D0)
     RXDirC(1) = (IXDirectionX/3.0D0)+&
          (IXDirectionY/3.0D0)-&
          (2.0D0*(IXDirectionZ/3.0D0))
     RXDirC(2) = (-IXDirectionX/3.0D0)+&
          (2.0D0*(IXDirectionY/3.0D0))-&
          (IXDirectionZ/3.0D0)
     RXDirC(3) = (IXDirectionX/3.0D0)+&
          (IXDirectionY/3.0D0)+&
          (IXDirectionZ/3.0D0)
     RNormDirC(1) = (INormalDirectionX/3.0D0)+&
          (INormalDirectionY/3.0D0)-&
          (2.0D0*(INormalDirectionZ/3.0D0))
     RNormDirC(2) = (-INormalDirectionX/3.0D0)+&
          (2.0D0*(INormalDirectionY/3.0D0))-&
          (INormalDirectionZ/3.0D0)
     RNormDirC(3) = (INormalDirectionX/3.0D0)+&
          (INormalDirectionY/3.0D0)+&
          (INormAlDirectionZ/3.0D0)

     DO ind =1,3
        IF(ABS(RZDirC(ind)).LE.TINY) THEN
           RZDirC(ind) = 100000000.0D0 ! A large number
        END IF
        IF(ABS(RXDirC(ind)).LE.TINY) THEN
           RXDirC(ind) = 100000000.0D0 ! A large number
        END IF
        IF(ABS(RNormDirC(ind)).LE.TINY) THEN
           RNormDirC(ind) = 100000000.0D0 ! A large number
        END IF
     END DO
        RZDirC = RZDirC/MINVAL(ABS(RZDirC))
        RXDirC = RXDirC/MINVAL(ABS(RXDirC))

     DO ind =1,3
        IF(RZDirC(ind).GT.10000000.0D0) THEN
           RZDirC(ind) = ZERO ! A large number
        END IF
        IF(RXDirC(ind).GT.10000000.0D0) THEN
           RXDirC(ind) = ZERO ! A large number
        END IF
        IF(RNormDirC(ind).GT.10000000.0D0) THEN
           RNormDirC(ind) = ZERO ! A large number
        END IF
     END DO
     
     RZDirC = REAL(NINT(RZDirC))
     RXDirC = REAL(NINT(RXDirC))
     RNormDirC = REAL(NINT(RNormDirC))
     
  END IF

!!$  IF(my_rank.EQ.0) THEN
!!$     PRINT*,RZDirC,RXDirC
!!$  END IF

  !RYDirC = CROSS(RZDirC,RYDirC)

  CALL CrystalLatticeVectorDetermination(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
             " in CrystalLatticeVectorDetermination"
     RETURN
  ENDIF

  CALL CrystalFullFractionalAtomicPostitionsCalculation(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
             " in CrystalFullFractionalAtomicPostitionsCalculation "
     RETURN
  ENDIF

  CALL CrystalUniqueFractionalAtomicPostitionsCalculation(IErr)
   IF( IErr.NE.0 ) THEN
     PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
             " in CrystalUniqueFractionalAtomicPostitionsCalculation "
     RETURN
  ENDIF

END SUBROUTINE CrystallographyInitialisation
