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
# $Id: makefile.GF,v 1.2 2014/03/27 18:13:30 phsht Exp $
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unixOS := $(shell uname)

ifeq ($(unixOS),Linux)

 linux := $(shell lsb_release -si)

 ifeq ($(linux),SUSE LINUX) #CSC machines

  modul := $(bash -c "module load ompi/1.6.4/gnu/4.3.4")

  INCLUDEDIR = /warwick/mathlib/gnu/x86_64/ # CSC machines

  LIBFLAGS= -llapack -lblas -lfftw3 

  ifeq ($(compiler),ifort) #change ifort flags and test
   F90 = mpif90
   FC = ifort -ffixed-form -std=legacy
   F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none -fcray-pointer -fexternal-blas
   F77FLAGS = -fbounds-check -fno-range-check

  else
   F90 = mpif90
   FC = gfortran -ffixed-form -std=legacy 
   F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none -fcray-pointer -fexternal-blas # debug
   #F90FLAGS = -pg -fbounds-check -fno-range-check -fimplicit-none # profiling
   F77FLAGS = -fbounds-check -fno-range-check # debug
  endif
 endif

 ifeq ($(linux),Ubuntu) #Ubuntu machine

  INCLUDEDIR = /usr/share/doc

  LIBFLAGS = -llapack -lblas -lfftw3
 
  F90 = mpif90
  F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none -fcray-pointer -fexternal-blas
  FC = gfortran -ffixed-form -std=legacy
  F77FLAGS = -fbounds-check -fno-range-check
 endif

 ifeq ($(linux),openSUSE project) #openSUSE machine

  LIBFLAGS = -llapack -lblas -lfftw3

  INCLUDEDIR = /usr/lib64

  F90 = mpif90
  F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none -fcray-pointer -fexternal-blas
  FC = gfortran -ffixed-form -std=legacy
  F77FLAGS = -fbounds-check -fno-range-check
 endif

 ifeq ($(linux),CentOS) #CentOS machine
 
  LIBFLAGS = -llapack -lblas -lfftw3

  INCLUDEDIR = /usr/lib64/

  F90 = mpif90
  F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none -fcray-pointer -fexternal-blas
  FC = gfortran -ffixed-form -std=legacy
  F77FLAGS = -fbounds-check -fno-range-check
 endif
endif 

ifeq ($(unixOS),Darwin) #MAC

 F90 = mpifort
 FC = gfortran -ffixed-form -std=legacy
 F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none
 #F90FLAGS = -pg -fbounds-check -fno-range-check -fimplicit-none
 F77FLAGS = -fbounds-check -fno-range-check # debug

 INCLUDEDIR = /usr/local/lib/ # RAR Linux

 LIBFLAGS= -llapack -lblas -lfftw3
endif

ifeq ($(unixOS),Cygwin) #needs to be tested
 F90 = mpif90
 FC = gfortran -ffixed-form -std=legacy
 F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none # debug
 #F90FLAGS = -pg -fbounds-check -fno-range-check -fimplicit-none # profiling
 F77FLAGS = -fbounds-check -fno-range-check # debug

 INCLUDEDIR = /lib/ # cygwin default

 LIBFLAGS= -llapack -lblas -lfftw3
endif

ifeq ($(unixOS),minerva) #needs to be tested

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
endif

ifeq ($(unixOS),ARCHER) #needs to be tested

 F90 = ftn
 FC = ftn
 F90FLAGS = # release
 #F90FLAGS = -fbounds-check -fno-range-check -fimplicit-none -fcray-pointer # debug
 #F90FLAGS = pg -fbounds-check -fno-range-check -fimplicit-none # profiling
 F77FLAGS = # debug
 #F77FLAGS = -fbounds-check -fno-range-check # debug

 #INCLUDEDIR = # ARCHER
 #INCLUDEDIR = /warwick/mathlib/gnu/x86_64/ # CSC machines
 #INCLUDEDIR = /usr/lib/ # RAR Linux

 #LIBFLAGS= -llapack -lblas -lfftw3 # CSC
 #LIBFLAGS= -lcray-libsci
endif

# NOTE if linking against $(NAG_SHARED) you must set the environment variable
# LD_LIBRARY_PATH=/warwick/mathlib/nag/ifort-x86_64/lib:$LD_LIBRARY_PATH
# in order to run the code

#object files for felixsim
SIMOBJECTFILES = gmodules.o smodules.o message.o errorcodes.o felixsim.o util.o in.o out.o readinput.o \
readciffile.o experimentalsetup.o imagesetup.o montagesetup.o refineutils.o\
structurefactorsetup.o ciftbx.o hash_funcs.o eigen.o invert.o microscopy.o \
diffractionpatterndefinitions.o diffractionpatterninitialisation.o crystallographyinitialisation.o \
crystallography.o bloch.o Ug.o image.o errorchecks.o writeoutput.o dqng.o xerror.o d1mach.o

#object files for felix refine
REFINEOBJECTFILES = gmodules.o smodules.o message.o errorcodes.o felixrefine.o util.o in.o out.o readciffile.o \
experimentalsetup.o imagesetup.o montagesetup.o structurefactorsetup.o readinput.o\
eigen.o invert.o microscopy.o diffractionpatterndefinitions.o diffractionpatterninitialisation.o \
crystallographyinitialisation.o crystallography.o bloch.o Ug.o image.o ciftbx.o simplex.o\
hash_funcs.o felixfunction.o refineutils.o writeoutput.o dqng.o xerror.o d1mach.o errorchecks.o\
symmetry.o

#for the object files which depend on .f90/.F90 files compile with mpif90
%.o: %.f90
	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -c $< $(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F90
	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)
#for object files which depend on .f/.F files compile with gfortran and the F77FLAGS
%.o: %.f
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

all:	felixsim felixrefine
#compiled with mpif90 
felixsim: $(SIMOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(SIMOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(SIMOBJECTFILES) $(LIBFLAGS)
#compiled with mpif90
felixrefine: $(REFINEOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(REFINEOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(REFINEOBJECTFILES) $(LIBFLAGS)


print:	
	a2ps -E -o print.ps gmodules.f90 smodules.f90 message.f90 felixsim.f90 util.f90 \
in.f90 out.f90 readinput.f90 readciffile.f90 experimentalsetup.f90 imagesetup.f90 \
montagesetup.f90 structurefactorsetup.f90 eigen.f90 invert.f90 microscopy.f90 \
diffraction.f90 simplex.f diffractionpatterninitialisation.f90 crystallographyinitialisation.f90 \
crystallography.f90 bloch.f90 Ug.f90 image.f90 ciftbx.f hash_funcs.f felixfunction.f90 \
writeoutput.f90 dqng.f xerror.f d1mach.f felixdraw.f90 felixrefine.f90 refineutils.f90 symmetry.f90\
makefile-CSC.GF; convert -density 150 print.ps print.pdf
	     	   	  
clean:	
	rm -f core *.mod *.o

#clean all:	
#	rm -f core *.mod *.o *.exe
