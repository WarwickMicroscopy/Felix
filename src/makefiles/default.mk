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

# default EXE name
MAIN=felix

# -----------------------------------------------------

FELIX=\
$(DIRFELIX)$(PRECISION)read_dm3_mod.o \
$(DIRFELIX)$(PRECISION)gmodules.o \
$(DIRFELIX)$(PRECISION)smodules.o \
$(DIRFELIX)$(PRECISION)l_alert_mod.o \
$(DIRFELIX)$(PRECISION)message_mod.o \
$(DIRFELIX)$(PRECISION)utilities_mod.o \
$(DIRFELIX)$(PRECISION)read_cif_mod.o \
$(DIRFELIX)$(PRECISION)read_files_mod.o \
$(DIRFELIX)$(PRECISION)setup_reflections_mod.o \
$(DIRFELIX)$(PRECISION)set_scatter_factors_mod.o \
$(DIRFELIX)$(PRECISION)ug_matrix_mod.o \
$(DIRFELIX)$(PRECISION)crystallography_mod.o \
$(DIRFELIX)$(PRECISION)bloch_mod.o \
$(DIRFELIX)$(PRECISION)image_initialisation_mod.o \
$(DIRFELIX)$(PRECISION)setup_space_group_mod.o \
$(DIRFELIX)$(PRECISION)write_output_mod.o \
$(DIRFELIX)$(PRECISION)refinementcontrol_mod.o \
$(DIRFELIX)$(PRECISION)simplex_mod.o \
$(DIRFELIX)$(PRECISION)felixrefine.o \

QUADPACK= \
  $(DIRQUADPACK)$(PRECISION)quadpack_double.o

CIFTBX=\
$(DIRCIFTBX)$(PRECISION)ciftbx.o $(DIRCIFTBX)$(PRECISION)hash_funcs.o

SAMPLES=$(MAIN).o

#STARTDIR=$(PWD)
#MYSTARTDIR=$(STARTDIR)

# Linux
#LIBS= -lblas -llapack -lfftw3

# HP alpha
#LIBS=-ljadamilu -llapack -lcxml -lblas -lm -lc -lfor

# IBM AIX
#LIBS=-ljadamilu -llapack -lblas  -lm -lc

# SGI Altix
#LIBS=-ljadamilu  -lmkl_lapack -lmkl  -lm -lc  -lifcore -lguide



# where are the headers
INCDIR=$(MYSTARTDIR)/include

# where are the libraries
LIBDIR=$(STARTDIR)/../lib/$(PLATFORM)

LIBFELIX=/lib$(PRECISION)felix.a
LIBQUADPACK=lib$(PRECISION)quadpack.a
LIBCIFTBX=lib$(PRECISION)ciftbx.a

.SUFFIXES: .c .f .f90 .F .o .a
.DEFAULT: main

#%.o: %.c
#	$(CC)  $(CCFLAGS)  -I$(INCDIR)  $(FORTRANNAMES) $(ARITHMETIC) $(LONGINTEGER) -c -o $(PRECISION)$@ $<
#
#
#%.o: %.F
#	$(FCOMPILE)

$(PRECISION)%.o: %.f90
#	@ echo --- --- compiling F90+ code
	$(F90) $(F90FLAGS) -c -o $@ $<

$(PRECISION)%.o: %.f
#	@ echo --- --- compiling F77 code
	$(FF) $(FFFLAGS) -c -o $@ $<

$(PRECISION)%.o: %.c
	$(CC)  $(CCFLAGS)  -I$(INCDIR)  $(FORTRANNAMES) $(ARITHMETIC) $(LONGINTEGER) -c -o $@ $<

$(PRECISION)%.o: %.F
	$(FCOMPILE)

#exe: $(STARTDIR)/$(MAIN).$(PLATFORM).$(PRECISION)

#$(STARTDIR)/$(MAIN).$(PLATFORM).$(PRECISION): $(FELIX)
#	@ echo --- linking $(MAIN) executable
#	$(LD) $(LDFLAGS) -o $@ $(FELIX) -L$(LIBDIR) -l$(LIBQUADPACK)  -l$(LIBCIFTBX)


