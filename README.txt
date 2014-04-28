!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! BlochSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
!
! (C) 2013/14, all rights reserved
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This file is part of BlochSim.
!
!  BlochSim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  BlochSim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with BlochSim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

README.txt file for FelixSim bloch wave method diffraction pattern simulation software.

Input File:

IWriteFlag: Controls the amount of information printed to the screen or
logfile during simulation, values 0-11 increase amount with 0 resulting in no
information being printed appart from start and stop messages.  More output
allows for easier debugging but will slow the execution.

IImageFLAG: Determines what images will be produced by the software :

0 Montages or Diffraction patterns only, one per thickness in .txt file format
1 Montages and Individual reflections saved in .txt file format, Reflections
will be bundled into seperate folders for each thickness.
2 Montages, reflections and amplitude/phase images will be saved, phase and
amplitude images will be saved individually for each reflection and labelled
-P or -A respectively

IOutputFLAG: Determines the amount of calculated Variables which are saved
for later use by FelixDraw.

0 Nothing is saved (Fastest)
1 Eigenspectra will be saved in binary format
2 Eigenspectra and Ug Matrix is saved
3 Eigenspectra, Ug Matrix and Wavefunctions saved

IScatteringFactorMethod: Determines which method by which to calculate potentials
0 Kirkands Method

IZolzFLAG: Choose to limit the simulation to the zeroth order laue zone.  

0 No (Includes HOLZ, slower) 
1 Yes (ZOLZ only, faster)

ICentralBeamFlag: Exclude the [000] beam from the final images which improves
relative intensity

0 No central beam
1 Central beam included

IMaskFLAG: Chooses between a circular or square input beam

0 Square
1 Circular

IAbsorbFLAG: Choose to include absorption in simulation

0 no
1 yes

IAnisotropicDebyeWallerFLAG: Choose to use anisotropic debye waller factors if
available

0 no
1 yes

IPixelCount: Pixel Radius of images, simulation scales as the square of this
  number but primary parallelisation is over pixels (more pixels, more cores
  can be used effectively) 64 is a good for images, 128+ better for
  quantitative analysis

RBSMaxGVecAmp: Maximum gvector magnitude (measured in reciprocal angstroms)
  allowed in the simulation, this defines the number "nreflections" and
  controls the size of the reflection pool accessible during diagonalisation
  (must have a decimal point)

RBSMaxDeviationPara: Maximum magnitude of the deviation parameter (distance
  between the reciprocal lattive point and ewald sphere) this value is
  unsigned and allows points within and outide the ewald sphere to be
  accepted. This paramter controls the number of beams overwhich the
  diagonisation is preformed, increasing it will result in greater accuracy
  but will slow simultion

RBSBethePara: Minimum Bethe potential allowed, this controls which
  nreflections (which been rejected as strong) are considered weak and used to
  perturb the potentials priorto diagonisations, a smaller number increases
  the number of weak beams

RDebyeWallerFactorConstant: If no Debye waller factor is found in the .cif
  file this value will be used for all atomic species missing the factor, note
  this is the B factor not U

RAbsorptionPer: Defines the percentage of absorption applied to the
  potentials, this values is used for all atomic species

RConvergenceAngle: Defines the convergence angle of the beam in units of half
  the minimum gvector magnitude, at a value of 1 all beams will touch at
  tjtheir edges, values greater than 1 would (experimentally) cause beams to
  overlap, in the simulation this will not occur

IIncidentBeamDirectionX: X Component of the incident beam direction (Zone
  axis) expressed in the crystal reference frame in real space

IIncidentBeamDirectionY: Y Component of the incident beam direction (Zone
  axis) expressed in the crystal reference frame in real space

IIncidentBeamDirectionZ: Z Component of the incident beam direction (Zone
  axis) expressed in the crystal reference frame in real space

IXDirectionX: X component of the chosen X-axis expressed in the crystal
  reference frame in reciprocal space

IXDirectionY: Y component of the chosen X-axis expressed in the crystal
  reference frame in reciprocal space

IXDirectionZ: Z component of the chosen X-axis expressed in the crystal
  reference frame in reciprocal space

RAccelerationVoltage: Acceleration voltage of the microscope expressed in KV

RInitialThickness: Lower bound thickness to be applied (Angstroms)

RFinalThickness: Upper Bound Thickness to be Appliced (Angstroms)

RDeltaThickness: Step between thickness (Angstroms)

IReflectOut: The number of the reflections to be included in the final
  image(s)
