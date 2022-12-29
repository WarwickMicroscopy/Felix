#!/bin/csh

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# For Warwick SCRTP:
#
# intel compiler: module load intel imkl imp
# gfortran compiler: module load GCC/11.2.0 OpenMPI/4.1.1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#
# makefile for FELIX
# 
MYPRECISION="d"
STARTDIR=$(PWD)
MYSTARTDIR=$(STARTDIR)

include user.mk
MYPLATFORM=$(PLATFORM)
include makefiles/platforms.mk
include makefiles/default.mk

DIRCIFTBX=ciftbx/
DIRQUADPACK=quadpack/
DIRFELIX=felix/

.PHONY:  felix felixdir felixrefine ciftbx quadpack \
         clean clear exe main

all: gitversion ciftbx quadpack felixdir ourlibs felix

# setting version
gitversion: 
	@ echo --- GETTING GIT VERSION INFORMATION
	@ cd $(DIRFELIX);\
	sed "s/GITVERSION/`git describe --tags --long`/g" version_template.txt | \
	sed "s/GITBRANCH/`git branch --show-current`/g" | \
	sed "s/COMPILED/`date`/g" >version.txt

# compiling libraries
quadpack:
	@ echo --- COMPILING QUADPACK integration routines
	@ cd $(DIRQUADPACK); \
        $(MAKE) "ARITHMETIC=-D_DOUBLE_REAL_" "PRECISION=$(MYPRECISION)" "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)"

ciftbx:
	@ echo --- COMPILING CIFTBX routines
	@ cd $(DIRCIFTBX); \
        $(MAKE) "ARITHMETIC=-D_DOUBLE_REAL_" "PRECISION=$(MYPRECISION)" "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)"

# compiling executable
felixdir:
	@ echo --- COMPILING felix code
	@ cd $(DIRFELIX);\
        $(MAKE) "ARITHMETIC=-D_DOUBLE_REAL_" "PRECISION=$(MYPRECISION)" "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)"

# binding libraries	
ourlibs: quadpack ciftbx
	@ echo --- BINDING libraries
	@ cd ../lib;\
        $(MAKE) ourlibs \
        "MYSTARTDIR=$(STARTDIR)" "PRECISION=$(MYPRECISION)" "MYPLATFORM=$(PLATFORM)"\
        "BLAS=" \
        "LAPACK="

# linking executable

FINALEXE=$(MAIN).$(PLATFORM).$(MYPRECISION)
	
felix: $(FINALEXE)

$(FINALEXE):
	@ echo --- LINKING felix code
	@ cd $(DIRFELIX);\
	$(F90) $(F90FLAGS) -o ../$(MAIN).$(PLATFORM).$(MYPRECISION) *.o \
	-L$(LIBDIR) -l$(MYPRECISION)quadpack -l$(MYPRECISION)ciftbx $(LIBS)
#        $(MAKE) exe "MYSTARTDIR=$(STARTDIR)" "PRECISION=$(MYPRECISION)" \
#	"MYPLATFORM=$(PLATFORM)" \
#	"MYLIBCIFTBX=$(LIBCIFTBX)" "MYLIBQUADPACK=$(LIBQUADPACK)"
	@ echo executable $(FINALEXE) can be found in the current directory

print:
	a2ps -E --pro=color --chars-per-line=100 -o - felix/*.f90 | ps2pdf - felix.pdf

clean:
	rm -rf   */*.o */*.mod

clear:
	rm -rf */*.o */*.mod */*/*.o */D/* */Z/* */S/* */C/*;\
        rm -rf ../lib/$(PLATFORM)/lib*.a
