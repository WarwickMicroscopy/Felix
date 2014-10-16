#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# felixsim
#
# Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
#
# (C) 2013/14, all rights reserved
#
# Version: :VERSION:
# Date:    :DATE:
# Time:    :TIME:
# Status:  :RLSTATUS:
# Build:   :BUILD:
# Author:  :AUTHOR:
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#  This file is part of felixsim.
#
#  felixsim is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  felixsim is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with felixsim.  If not, see <http://www.gnu.org/licenses/>.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# $Id: makefile-CSC.GF,v 1.2 2014/03/27 18:13:30 phsht Exp $
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F90 = mpif90
FC = gfortran -ffixed-form -std=legacy
F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none -fcray-pointer # debug
#F90FLAGS = -pg -fbounds-check -fno-range-check -fimplicit-none # profiling
F77FLAGS = -fbounds-check -fno-range-check # debug

INCLUDEDIR = /warwick/mathlib/gnu/x86_64/ # CSC machines
#INCLUDEDIR = /usr/lib/ # RAR Linux

LIBFLAGS= -llapack -lblas -lfftw3

# NOTE if linking against $(NAG_SHARED) you must set the environment variable
# LD_LIBRARY_PATH=/warwick/mathlib/nag/ifort-x86_64/lib:$LD_LIBRARY_PATH
# in order to run the code

SIMOBJECTFILES = gmodules.o smodules.o main.o util.o inout.o inpcif.o ciftbx.o hash_funcs.o eigen.o invert.o microscopy.o diffraction.o crystallography.o bloch.o dqng.o xerror.o d1mach.o Ug.o image.o 
DRAWOBJECTFILES = gmodules.o smodules.o lacbed.o util.o inout.o \
eigen.o invert.o microscopy.o diffraction.o crystallography.o bloch.o Ug.o image.o inpcif.o ciftbx.o hash_funcs.o dqng.o xerror.o d1mach.o

REFINEOBJECTFILES = gmodules.o smodules.o refinemain.o util.o inout.o inpcif.o ciftbx.o hash_funcs.o eigen.o invert.o microscopy.o diffraction.o crystallography.o bloch.o Ug.o image.o refineutils.o dqng.o xerror.o d1mach.o

%.o: %.f90
	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -c $< $(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F90
	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.f
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

all:	felixsim felixdraw felixrefine

felixsim: $(SIMOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(SIMOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(SIMOBJECTFILES) $(LIBFLAGS)

felixdraw: $(DRAWOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(DRAWOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(DRAWOBJECTFILES) $(LIBFLAGS)

felixrefine: $(REFINEOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(REFINEOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(REFINEOBJECTFILES) $(LIBFLAGS)


print:	
	a2ps -E -o print.ps gmodules.f90 smodules.f90 random.f90 main.f90 util.f90 inout.f90 inpcif.f90 eigen.f90 invert.f90 microscopy.f90 diffraction.f90 crystallography.f90 bloch.f90 Ug.f90 image.f90 refinemain.f90 dqng.f xerror.f d1mach.f makefile; convert -density 150 print.ps print.pdf

clean:	
	rm -f core *.mod *.o
