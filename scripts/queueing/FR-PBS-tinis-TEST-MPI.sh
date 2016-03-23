#!/bin/bash
#PBS -l nodes=2:ppn=12
#PBS -l pmem=2012mb
#PBS -l walltime=00:30:00
#PBS -M $USER@warwick.ac.uk
#PBS -m a
#PBS -r y
#PBS -V
#PBS -k oe
#PBS -j oe

module load intel impi imkl FFTW

pwd
cd ~/Projects/D-LACBED/Felix/samples/FR-GaAs/
pwd
ls
srun -n 24 ../../src/felixrefine
