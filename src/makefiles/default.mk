
# default EXE name
MAIN=felix

# -----------------------------------------------------

FELIX=\
$(DIRFELIX)$(PRECISION)gmodules.f90 \
$(DIRFELIX)$(PRECISION)smodules.f90 \
$(DIRFELIX)$(PRECISION)alert_mod.f90 \
$(DIRFELIX)$(PRECISION)message_mod.f90 \
$(DIRFELIX)$(PRECISION)utilities_mod.f90 \
$(DIRFELIX)$(PRECISION)read_cif_mod.f90 \
$(DIRFELIX)$(PRECISION)read_mod.f90 \
$(DIRFELIX)$(PRECISION)setup_reflections_mod.f90 \
$(DIRFELIX)$(PRECISION)setup_scattering_factors_mod.f90 \
$(DIRFELIX)$(PRECISION)crystallography_mod.f90 \
$(DIRFELIX)$(PRECISION)Ug_mod.f90 \
$(DIRFELIX)$(PRECISION)bloch_mod.f90 \
$(DIRFELIX)$(PRECISION)felixfunction_mod.f90 \
$(DIRFELIX)$(PRECISION)felixrefine.f90 \
$(DIRFELIX)$(PRECISION)image.f90 \
$(DIRFELIX)$(PRECISION)simplex.f90 \
$(DIRFELIX)$(PRECISION)RefineWriteOut.f90 \
$(DIRFELIX)$(PRECISION)symmetry.f90 \


QUADPACK= \
  $(DIRQUADPACK)$(PRECISION)quadpack_double.o

CIFTBX=\
$(DIRCIFTBX)$(PRECISION)ciftbx.o $(DIRCIFTBX)$(PRECISION)hash_funcs.o

SAMPLES=$(MAIN).o

#STARTDIR=$(PWD)
#MYSTARTDIR=$(STARTDIR)

# Linux
LIBS= -lblas -llapack -lfftw3

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


