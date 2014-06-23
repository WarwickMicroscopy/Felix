#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Felix
#
# Richard Beanland, Keith Evans and Rudolf A Roemer
#
# (C) 2013/14, all rights reserved
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#  This file is part of FelixSim+Draw.
#
#  BlochSim is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  BlochSim is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with BlochSim.  If not, see <http://www.gnu.org/licenses/>.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# $Id: makefile-minerva.GF,v 1.5 2014/03/27 15:27:24 phsht Exp $
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#module swap oldmodules minerva-2.0
#module load goolf-1.7.10

F90 = mpif90
FC = gfortran -ffixed-form -std=legacy
F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none # debug
#F90FLAGS = -m64 -fno-range-check -pg -ffree-line-length-none # profiling
F77FLAGS = -fbounds-check -fno-range-check # debug

#INCLUDEDIR = /hpcwarwick/mathlib/gnu/lib/ # CSC machines
#INCLUDEDIR = /usr/lib/ # RAR Linux

#LIBFLAGS= -L/hpcwarwick/mathlib/gnu/lib/ -lgoto2_nehalem-r1.13 #-llapack -lblas
#LIBFLAGS= -lmkl
LIBFLAGS= \
-L/csc/minerva/software/gompi/1.7.10/ScaLAPACK/2.0.2-openblas-0.2.8-lapack-3.4.2/lib \
-L/csc/minerva/software/gompi/1.7.10/OpenBLAS/0.2.8-lapack-3.4.2/lib \
-lopenblas \
-lfftw3

#LD_LIBRARY_PATH=/hpcwarwick/mathlib/gnu/lib/:$LD_LIBRARY_PATH

SIMOBJECTFILES = gmodules.o smodules.o main.o util.o inout.o inpcif.o ciftbx.o hash_funcs.o eigen.o invert.o microscopy.o diffraction.o crystallography.o bloch.o Ug.o image.o
DRAWOBJECTFILES = gmodules.o smodules.o lacbed.o util.o inout.o \
eigen.o invert.o microscopy.o diffraction.o crystallography.o bloch.o Ug.o image.o inpcif.o ciftbx.o hash_funcs.o 

REFINEOBJECTFILES = gmodules.o smodules.o refinemain.o util.o inout.o inpcif.o ciftbx.o hash_funcs.o eigen.o invert.o microscopy.o diffraction.o crystallography.o bloch.o Ug.o image.o


%.o: %.f90
	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -c $< $(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F90
	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.f
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

all:	FelixSim FelixDraw FelixRefine

FelixSim: $(SIMOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(SIMOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(SIMOBJECTFILES) $(LIBFLAGS)

FelixDraw: $(DRAWOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(DRAWOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(DRAWOBJECTFILES) $(LIBFLAGS)

FelixRefine: $(REFINEOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(REFINEOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(REFINEOBJECTFILES) $(LIBFLAGS)

clean:	
	rm -f core *.mod *.o
