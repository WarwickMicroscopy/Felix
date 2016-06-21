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

SUBROUTINE MicroscopySettings( IErr )

  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  REAL(RKIND)::norm,dummy,ROneThousand
  INTEGER :: IErr

  CALL Message("MicroscopySettings",IMust,IErr)

  ROneThousand = 1000.0_RKIND

  !Electron velocity in metres per second
  RElectronVelocity= RSpeedOfLight*SQRT(ONE-((RElectronMass*RSpeedOfLight**2)/ &
       (RElectronCharge*RAcceleratingVoltage*ROneThousand + &
        RElectronMass*RSpeedOfLight**2) )**2 )
  
  !Electron wavelength in Angstroms
  RElectronWaveLength= RPlanckConstant / &
       ( SQRT(TWO*RElectronMass*RElectronCharge*RAcceleratingVoltage*ROneThousand) * &
         SQRT(ONE + (RElectronCharge*RAcceleratingVoltage*ROneThousand) / &
          (TWO*RElectronMass*RSpeedOfLight**2) ))* RAngstromConversion
  
  !(NB k=1/lamda and exp(2*pi*i*q.r), optical convention)
  RElectronWaveVectorMagnitude=ONE/RElectronWaveLength
  RRelativisticCorrection= ONE/SQRT( ONE - (RElectronVelocity/RSpeedOfLight)**2 )
  RRelativisticMass= RRelativisticCorrection*RElectronMass
   
  CALL Message("MicroscopySettings",IInfo,IErr, &
       MessageVariable = "ElectronVelocity",RVariable = RElectronVelocity)
  CALL Message("MicroscopySettings",IInfo,IErr, &
       MessageVariable = "ElectronWavelength",RVariable = RElectronWaveLength)
  CALL Message("MicroscopySettings",IInfo,IErr, &
       MessageVariable = "ElectronWaveVectorMagnitude",RVariable = RElectronWaveVectorMagnitude)
  CALL Message("MicroscopySettings",IInfo,IErr, &
       MessageVariable = "RelativisticCorrection",RVariable = RRelativisticCorrection)
     
  RETURN

END SUBROUTINE MicroscopySettings
