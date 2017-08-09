#!/bin/csh
#
# makefile for FELIX
# 
STARTDIR=$(PWD)
MYSTARTDIR=$(STARTDIR)

include user.mk
MYPLATFORM=$(PLATFORM)
include makefiles/platforms.mk
include makefiles/default.mk

DIRCIFTBX=ciftbx/
DIRQUADPACK=quadpack/
DIRFELIX=felixdir/

.PHONY:  felix felixdir ciftbx quadpack \
         clean clear 

all: felix felix ciftbx quadpack

# compiling libraries
quadpack:
	@ echo --- build QUADPACK integration routines
	@ cd $(DIRQUADPACK); \
        $(MAKE) "ARITHMETIC=-D_DOUBLE_REAL_" "PRECISION=d" "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)"

ciftbx:
	@ echo --- build CIFTBX routines
	@ cd $(DIRCIFTBX); \
        $(MAKE) "ARITHMETIC=-D_DOUBLE_REAL_" "PRECISION=d" "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)"

# compiling executable
felix:
	@ echo --- build felixrefine code
	@ cd $(DIRFELIX);\
        $(MAKE) "ARITHMETIC=-D_DOUBLE_REAL_" "PRECISION=d" "MYPLATFORM=$(PLATFORM)" \
        "MYSTARTDIR=$(STARTDIR)"
	@ echo sample executable program can be found in the current directory

# binding libraries	
ourlibs: quadpack ciftbx
	@ echo bind libraries
	@ cd ../lib;\
        $(MAKE) ourmain \
        "MYSTARTDIR=$(STARTDIR)" \
        "BLAS=" \
        "LAPACK="

clean:
	rm -rf   */*.o


clear:
	rm -rf */*.o */*/*.o */D/* */Z/* */S/* */C/*;\
        rm -rf ../lib/$(PLATFORM)/lib*.a
