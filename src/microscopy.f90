!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
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

SUBROUTINE MicroscopySettings( IErr )

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  REAL(RKIND) norm,dummy

  INTEGER IErr

  IF(((IWriteFLAG.GE.0).AND.(my_rank.EQ.0)).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"MicroscopySettings(",my_rank,")"
  END IF

  RElectronVelocity= &
       RSpeedOfLight * &
       SQRT( 1.0D0 - ( &
       (RElectronMass*RSpeedOfLight**2) / &
       (RElectronCharge*RAcceleratingVoltage*1.0D3 + &
        RElectronMass*RSpeedOfLight**2) &
       )**2 )

  RElectronWaveLength= RPlanckConstant / &
       ( SQRT(2.0D0*RElectronMass*RElectronCharge*RAcceleratingVoltage*1.D3) * &
         SQRT(1.0D0 + (RElectronCharge*RAcceleratingVoltage*1.D3) / &
          (2.D0*RElectronMass*RSpeedOfLight**2) &
       )) * RAngstromConversion

  RElectronWaveVectorMagnitude=TWOPI/RElectronWaveLength

  RRelativisticCorrection= &
       1.D0 / SQRT( 1.D0 - (RElectronVelocity/RSpeedOfLight)**2 )

  RRelativisticMass= RRelativisticCorrection*RElectronMass

  IF(((IWriteFLAG.GE.1).AND.(my_rank.EQ.0)).OR.IWriteFLAG.GE.10) THEN
     
     PRINT*,"MicroscopySettings(",my_rank,") ElectronVelocity =", RElectronVelocity
     PRINT*,"MicroscopySettings(",my_rank,") ElectronWaveLength =", RElectronWaveLength
     PRINT*,"MicroscopySettings(",my_rank,") ElectronWaveVectorMagnitude =", RElectronWaveVectorMagnitude
     PRINT*,"MicroscopySettings(",my_rank,") RelativisticCorrection =", RRelativisticCorrection

  END IF
     
  RETURN

END SUBROUTINE MicroscopySettings
