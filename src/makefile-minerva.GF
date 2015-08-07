#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# felix
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
#  This file is part of felixsim+draw.
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
# $Id: makefile-minerva.GF,v 1.5 2014/03/27 15:27:24 phsht Exp $
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#module swap oldmodules minerva-2.0
#module load goolf-1.7.10

F90 = mpif90
FC = gfortran -ffixed-form -std=legacy
F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none # debug
#F90FLAGS = -m64 -fno-range-check -pg -ffree-line-length-none # profiling
F77FLAGS = -fbounds-check -fno-range-check # debug

LIBFLAGS= \
-L/csc/minerva/software/gompi/1.7.10/ScaLAPACK/2.0.2-openblas-0.2.8-lapack-3.4.2/lib \
-L/csc/minerva/software/gompi/1.7.10/OpenBLAS/0.2.8-lapack-3.4.2/lib \
-lopenblas \
-lfftw3

#LD_LIBRARY_PATH=/hpcwarwick/mathlib/gnu/lib/:$LD_LIBRARY_PATH

SIMOBJECTFILES = gmodules.o smodules.o message.o errorcodes.o felixsim.o util.o in.o out.o readinput.o \
readciffile.o experimentalsetup.o imagesetup.o montagesetup.o refineutils.o \
structurefactorsetup.o ciftbx.o hash_funcs.o eigen.o invert.o microscopy.o \
diffractionpatterndefinitions.o diffractionpatterninitialisation.o crystallographyinitialisation.o \
crystallography.o bloch.o Ug.o image.o errorchecks.o writeoutput.o dqng.o xerror.o d1mach.o

REFINEOBJECTFILES = gmodules.o smodules.o inputmodules.o message.o errorcodes.o felixrefine.o util.o in.o out.o readciffile.o \
experimentalsetup.o imagesetup.o montagesetup.o structurefactorsetup.o readinput.o\
eigen.o invert.o microscopy.o diffractionpatterndefinitions.o diffractionpatterninitialisation.o \
crystallographyinitialisation.o crystallography.o bloch.o Ug.o image.o ciftbx.o simplex.o\
hash_funcs.o RefineWriteOut.o felixfunction.o refineutils.o writeoutput.o dqng.o xerror.o d1mach.o errorchecks.o\
symmetry.o

%.o: %.f90
	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -c $< $(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F90
	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.f
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

all:	felixsim felixrefine

felixsim: $(SIMOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(SIMOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(SIMOBJECTFILES) $(LIBFLAGS)

felixrefine: $(REFINEOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(REFINEOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(REFINEOBJECTFILES) $(LIBFLAGS)

print:	
	a2ps -E -o print.ps gmodules.f90 smodules.f90 inputmodule.f90 message.f90 errorcodes.f90 felixsim.f90 util.f90 \
in.f90 out.f90 readinput.f90 readciffile.f90 experimentalsetup.f90 imagesetup.f90 \
montagesetup.f90 structurefactorsetup.f90 eigen.f90 invert.f90 microscopy.f90 \
diffraction.f90 diffractionpatterninitialisation.f90 crystallographyinitialisation.f90 \
crystallography.f90 bloch.f90 Ug.f90 image.f90 ciftbx.f hash_funcs.f \
writeoutput.f90 dqng.f xerror.f d1mach.f felixdraw.f90 felixrefine.f90 refineutils.f90 \
makefile-minerva.GF; convert -density 150 print.ps print.pdf

clean:	
	rm -f core *.mod *.o

#clean all:	
#	rm -f core *.mod *.o *.exe
