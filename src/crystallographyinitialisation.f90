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
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr, ind
 
  CALL Message("CrystallographyInitialisation",IMust,IErr)
    
  IF(IDiffractionFLAG.EQ.1) THEN!RB what is this about
     DEALLOCATE(RFullAtomicFracCoordVec,SFullAtomicNameVec,&
          RFullPartialOccupancy,RFullIsotropicDebyeWallerFactor, &
          IFullAtomicNumber, IFullAnisotropicDWFTensor, &
          MNP,SMNP,RDWF,ROcc,IAtoms,IAnisoDWFT,STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
             " in DEALLOCATE()"
        RETURN
     ENDIF
  END IF
  
  IF(IPseudoCubicFLAG.EQ.1.AND.IDiffractionFLAG.EQ.0) THEN
     RZDirC(1) = (IIncidentBeamDirectionX/THREE)+&
          (IIncidentBeamDirectionY/THREE)-&
          (TWO*(IIncidentBeamDirectionZ/THREE))
     RZDirC(2) = (-IIncidentBeamDirectionX/THREE)+&
          (TWO*(IIncidentBeamDirectionY/THREE))-&
          (IIncidentBeamDirectionZ/THREE)
     RZDirC(3) = (IIncidentBeamDirectionX/THREE)+&
          (IIncidentBeamDirectionY/THREE)+&
          (IIncidentBeamDirectionZ/THREE)
     RXDirC(1) = (IXDirectionX/THREE)+&
          (IXDirectionY/THREE)-&
          (TWO*(IXDirectionZ/THREE))
     RXDirC(2) = (-IXDirectionX/THREE)+&
          (TWO*(IXDirectionY/THREE))-&
          (IXDirectionZ/THREE)
     RXDirC(3) = (IXDirectionX/THREE)+&
          (IXDirectionY/THREE)+&
          (IXDirectionZ/THREE)
     RNormDirC(1) = (INormalDirectionX/THREE)+&
          (INormalDirectionY/THREE)-&
          (TWO*(INormalDirectionZ/THREE))
     RNormDirC(2) = (-INormalDirectionX/THREE)+&
          (TWO*(INormalDirectionY/THREE))-&
          (INormalDirectionZ/THREE)
     RNormDirC(3) = (INormalDirectionX/THREE)+&
          (INormalDirectionY/THREE)+&
          (INormAlDirectionZ/THREE)

     DO ind =1,3
        IF(ABS(RZDirC(ind)).LE.TINY) THEN
           RZDirC(ind) = REAL(100000000.0,RKIND) ! A large number
        END IF
        IF(ABS(RXDirC(ind)).LE.TINY) THEN
           RXDirC(ind) = REAL(100000000.0,RKIND) ! A large number
        END IF
        IF(ABS(RNormDirC(ind)).LE.TINY) THEN
           RNormDirC(ind) = REAL(100000000.0,RKIND) ! A large number
        END IF
     END DO
        RZDirC = RZDirC/MINVAL(ABS(RZDirC))
        RXDirC = RXDirC/MINVAL(ABS(RXDirC))

     DO ind =1,3
        IF(RZDirC(ind).GT.REAL(10000000.0,RKIND)) THEN
           RZDirC(ind) = ZERO ! A large number
        END IF
        IF(RXDirC(ind).GT.REAL(10000000.0,RKIND)) THEN
           RXDirC(ind) = ZERO ! A large number
        END IF
        IF(RNormDirC(ind).GT.REAL(10000000.0,RKIND)) THEN
           RNormDirC(ind) = ZERO ! A large number
        END IF
     END DO
     
     RZDirC = REAL(NINT(RZDirC))
     RXDirC = REAL(NINT(RXDirC))
     RNormDirC = REAL(NINT(RNormDirC))
     
  END IF

  CALL CrystalLatticeVectorDetermination(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
             " in CrystalLatticeVectorDetermination"
     RETURN
  ENDIF

  CALL AllAtomPositions(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
             " in AllAtomPositions "
     RETURN
  ENDIF

  !zz CALL CrystalUniqueFractionalAtomicPostitionsCalculation(IErr)
!zz    IF( IErr.NE.0 ) THEN
!zz      PRINT*,"CrystallographyInitialisation(", my_rank, ") error ", IErr, &
!zz              " in CrystalUniqueFractionalAtomicPostitionsCalculation "
!zz      RETURN
!zz   ENDIF

END SUBROUTINE CrystallographyInitialisation
