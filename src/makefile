#!/bin/csh
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
DIRFELIX=felixdir/

.PHONY:  felix felixdir felixrefine ciftbx quadpack \
         clean clear exe main

all: ciftbx quadpack felixdir ourlibs felix

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
felix: $(MAIN).$(PLATFORM).$(MYPRECISION)

$(MAIN).$(PLATFORM).$(MYPRECISION):
	@ echo --- LINKING felixrefine code
	@ cd $(DIRFELIX);\
	$(F90) -o ../$(MAIN).$(PLATFORM).$(MYPRECISION) *.o \
	-L$(LIBDIR) -l$(MYPRECISION)quadpack -l$(MYPRECISION)ciftbx $(LIBS)
#        $(MAKE) exe "MYSTARTDIR=$(STARTDIR)" "PRECISION=$(MYPRECISION)" \
#	"MYPLATFORM=$(PLATFORM)" \
#	"MYLIBCIFTBX=$(LIBCIFTBX)" "MYLIBQUADPACK=$(LIBQUADPACK)"
#	@ echo sample executable program can be found in the current directory


clean:
	rm -rf   */*.o

clear:
	rm -rf */*.o */*/*.o */D/* */Z/* */S/* */C/*;\
        rm -rf ../lib/$(PLATFORM)/lib*.a
