#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Felix
#
# Richard Beanland, Keith Evans & Rudolf A Roemer
#
# (C) 2013-17, all rights reserved
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
#  Felix is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  Felix is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with Felix.  If not, see <http://www.gnu.org/licenses/>.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# $Id: makefile.GF,v 1.2 2014/03/27 18:13:30 phsht Exp $
#
# cygwin needs: gcc, gfortran, make, git, openmpi, lapack, fftw3
# CSC modules needed: ompi/2.0.1/gnu/4.9.2
# minerva modules needed: minerva-2.0, goolf-1.7.10
# ARCHER modules needed: cray-libsci/13.0.1, fftw, cray-mpich/7.1.1
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unixOS := $(shell uname -s)
cyg := $(shell expr substr $(unixOS) 1 6)

GLOBALgfortranF90FLAGS = -ffree-line-length-none -fbounds-check -fno-range-check -fimplicit-none # debug
GLOBALgfortranF77FLAGS = -fbounds-check -fno-range-check # debug


ifeq ($(unixOS),Linux)

 linux := $(shell (lsb_release -si))
 tiniscut =$(subst /, , $(MODULEPATH))
 tinis :=$(word 2, $(tiniscut))

  ifeq ($(linux),SUSE LINUX) #CSC machines
   HPC := $(shell expr substr $(shell uname -n) 1 3)
   ARCHER := $(shell expr substr $(shell uname -n) 1 2)
   
   ifeq ($(ARCHER),es) #ARCHER
   $(info Compiling on ARCHER)
   $(info Please ensure to run "module load fftw")

   F90 = ftn
   FC = ftn
   F90FLAGS = $(GLOBALgfortranF90FLAGS) -fcray-pointer
   F77FLAGS = $(GLOBALgfortranF90FLAGS)

   INCLUDEDIR = /opt # ARCHER

   #LIBFLAGS= -llapack -lblas -lfftw3 # CSC
   #LIBFLAGS= -lcray-libsci

   else #CSC

   $(info Compiling on CSC, $(linux))
   $(info Please ensure to run "module load ompi/2.0.1/gnu/4.9.2" )

   INCLUDEDIR = /warwick/mathlib/gnu/x86_64/

   LIBFLAGS= -llapack -lblas -lfftw3 

    ifeq ($(FOO),BAR) #ifort compiler
     F90 = mpif90
     FC = ifort 
     F90FLAGS =
     F77FLAGS =

    else
     F90 = mpif90
     FC = gfortran -ffixed-form -std=legacy 
     F90FLAGS = $(GLOBALgfortranF90FLAGS) -fcray-pointer -fexternal-blas 
     F77FLAGS = $(GLOBALgfortranF77FLAGS) 
    endif	     
   endif
 
 else ifeq ($(linux),Ubuntu) #Ubuntu machine
  $(info Compiling on $(linux), $(unixOS))

  INCLUDEDIR = /usr/share/doc
  LIBFLAGS = -llapack -lblas -lfftw3
  F90 = mpif90
  FC = gfortran -ffixed-form -std=legacy	
  F90FLAGS = $(GLOBALgfortranF90FLAGS) -fcray-pointer -fexternal-blas
  F77FLAGS = $(GLOBALgfortranF90FLAGS)

 else ifeq ($(linux),Debian) #Debian
  $(info Compiling on $(linux), $(unixOS))

  INCLUDEDIR = /usr/share/doc
  LIBFLAGS = -llapack -lblas -lfftw3
  F90 = mpif90
  FC = gfortran -ffixed-form -std=legacy
  F90FLAGS = $(GLOBALgfortranF90FLAGS) -fcray-pointer -fexternal-blas
  F77FLAGS = $(GLOBALgfortranF90FLAGS)

 else ifeq ($(linux),openSUSE project) #openSUSE machine
  $(info Compiling on $(linux), $(unixOS))

  LIBFLAGS = -llapack -lblas -lfftw3
  INCLUDEDIR = /usr/lib64
  F90 = mpif90
  FC = gfortran -ffixed-form -std=legacy
  F90FLAGS = $(GLOBALgfortranF90FLAGS) -fcray-pointer -fexternal-blas
  F77FLAGS = $(GLOBALgfortranF90FLAGS)

 else ifeq ($(linux),CentOS) #CentOS machine
  $(info Compiling on $(linux), $(unixOS))
 
  LIBFLAGS = -llapack -lblas -lfftw3
  INCLUDEDIR = /usr/lib64/
  F90 = /usr/lib64/openmpi/bin/mpif90
  FC = gfortran -ffixed-form -std=legacy
  F90FLAGS = $(GLOBALgfortranF90FLAGS) -fcray-pointer -fexternal-blas
  F77FLAGS = $(GLOBALgfortranF90FLAGS)

 else ifeq ($(linux),LinuxMint)
  $(info Compiling on $(linux))
   
  LIBFLAGS = -llapack -lblas -lfftw3
  INCLUDEDIR = /usr/lib/
  F90 = mpif90
  FC = gfortran -ffixed-form -std=legacy
  F90FLAGS = $(GLOBALgfortranF90FLAGS) -fcray-pointer -fexternal-blas
  F77FLAGS = $(GLOBALgfortranF90FLAGS)

 else ifeq ($(linux),Scientific) #Scientific Linux
  $(info Compiling on Scientific Linux)

  LIBFLAGS = -llapack -lblas -lfftw3
  INCLUDEDIR = /usr/lib64/
  F90 = /usr/lib64/openmpi/bin/mpif90
  FC = gfortran -ffixed-form -std=legacy
  F90FLAGS = $(GLOBALgfortranF90FLAGS) -fcray-pointer -fexternal-blas
  F77FLAGS = $(GLOBALgfortranF90FLAGS)

 else #default or tinis

  ifeq ($(tinis),tinis) #tinis hpc
   $(info Compiling on tinis)
   $(info please run: module load intel impi imkl FFTW)

   F90 = mpiifort
   FC = ifort -fixed
   F90FLAGS = -132 -r8 -u -C -g -traceback -check noarg_temp_created #debug
   #F90FLAGS = -132 -r8 -u -g -p #profiling
   F77FLAGS = -r8 #debug

   LIBFLAGS= \
   -mkl -mkl=cluster -qopenmp -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -lfftw3

   INCLUDEDIR= #/hpcwarwick/mathlib/gnu/lib

   else #default

    $(info Compiling on default)
  
   INCLUDEDIR = /warwick/mathlib/gnu/x86_64/

   LIBFLAGS = -llapack -lblas -lfftw3

    F90 = mpif90
    FC = gfortran -ffixed-form -std=legacy
    F90FLAGS = $(GLOBALgfortranF90FLAGS) -fcray-pointer -fexternal-blas
    F77FLAGS = $(GLOBALgfortranF90FLAGS)
   endif
endif

else ifeq ($(unixOS),Darwin) #MAC
 $(info Compiling on $(unixOS), OS X)

 F90 = mpifort-openmpi-mp
 FC = gfortran-mp-4.5 -ffixed-form -std=legacy
 F90FLAGS = $(GLOBALgfortranF90FLAGS) -fcray-pointer -fexternal-blas
 F77FLAGS = $(GLOBALgfortranF90FLAGS)

 INCLUDEDIR = /usr/local/lib/ # RAR Linux

 LIBFLAGS= -llapack -lblas -lfftw3

else ifeq ($(cyg),CYGWIN)
 $(info Compiling on $(cyg)) 

 F90 = mpif90
 FC = gfortran -ffixed-form -std=legacy
 F90FLAGS = $(GLOBALgfortranF90FLAGS) -fcray-pointer -fexternal-blas
 F77FLAGS = $(GLOBALgfortranF90FLAGS)

 INCLUDEDIR = /lib/ # cygwin default

 LIBFLAGS= -llapack -lblas -lfftw3

else #default
 $(info Compiling on default settings)
 
 INCLUDEDIR = /warwick/mathlib/gnu/x86_64

 LIBFLAGS = -llapack -lblas -lfftw3

 F90 = mpif90
 FC = gfortran -ffixed-form -std-legacy
 F90FLAGS = $(GLOBALgfortranF90FLAGS) -fcray-pointer -fexternal-blas
 F77FLAGS = $(GLOBALgfortranF90FLAGS)

endif

# NOTE if linking against $(NAG_SHARED) you must set the environment variable
# LD_LIBRARY_PATH=/warwick/mathlib/nag/ifort-x86_64/lib:$LD_LIBRARY_PATH
# in order to run the code

REFINEOBJECTFILES = gmodules.o smodules.o message.o inputmodules.o errorcodes.o felixrefine.o util.o scatteringfactors.o in.o out.o readcif.o \
experimentalsetup.o imagesetup.o montagesetup.o structurefactorsetup.o readinput.o\
eigen.o invert.o microscopy.o diffractionpatterndefinitions.o diffractionpatterninitialisation.o \
crystallographyinitialisation.o crystallography.o bloch.o Ug.o image.o ciftbx.o simplex.o\
hash_funcs.o RefineWriteOut.o felixfunction.o refineutils.o writeoutput.o errorchecks.o\
symmetry.o quadpack_double.f90
#dqng.o xerror.o d1mach.o 

%.o: %.f90
	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -c $< $(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F90
	$(F90) $(F90FLAGS) -DF90 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.f
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

%.o: %.F
	$(FC) $(F77FLAGS) -DF77 -c $< -I$(INCLUDEDIR) $(LIBFLAGS)

all: felixrefine

felixrefine: $(REFINEOBJECTFILES)
	$(F90) $(F90FLAGS) -o $@ $(REFINEOBJECTFILES) -I$(INCLUDEDIR) $(LIBFLAGS)
#	$(F90) -o $@ $(REFINEOBJECTFILES) $(LIBFLAGS)


print:	
	a2ps -E -o print.ps gmodules.f90 smodules.f90 message.f90 felixsim.f90 util.f90scatterfactors.f90 \
in.f90 out.f90 readinput.f90 readcif.f90 experimentalsetup.f90 imagesetup.f90 \
montagesetup.f90 structurefactorsetup.f90 eigen.f90 invert.f90 microscopy.f90 \
diffraction.f90 simplex.f diffractionpatterninitialisation.f90 crystallographyinitialisation.f90 \
crystallography.f90 bloch.f90 Ug.f90 image.f90 ciftbx.f hash_funcs.f RefineWriteOut.f90 felixfunction.f90 \
writeoutput.f90 dqng.f xerror.f d1mach.f felixdraw.f90 felixrefine.f90 refineutils.f90 symmetry.f90 quadpack_double.f90\
makefile-CSC.GF; convert -density 150 print.ps print.pdf
	     	   	  
clean:	
	rm -f core *.mod *.o *.bin

#clean all:	
#	rm -f core *.mod *.o *.exe
