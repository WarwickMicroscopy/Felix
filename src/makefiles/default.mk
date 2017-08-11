
# default EXE name
MAIN=felixrefine

# -----------------------------------------------------

FELIX=\
$(DIRFELIX)$(PRECISION)gmodules.o \
$(DIRFELIX)$(PRECISION)smodules.o \
$(DIRFELIX)$(PRECISION)message_mod.o \
$(DIRFELIX)$(PRECISION)felixrefine.o \
$(DIRFELIX)$(PRECISION)felixfunction.o \
$(DIRFELIX)$(PRECISION)bloch.o \
$(DIRFELIX)$(PRECISION)refineutils.o \
$(DIRFELIX)$(PRECISION)Ug.o \
$(DIRFELIX)$(PRECISION)crystallography.o \
$(DIRFELIX)$(PRECISION)image.o \
$(DIRFELIX)$(PRECISION)RefineWriteOut.o \
$(DIRFELIX)$(PRECISION)util.o \
$(DIRFELIX)$(PRECISION)diffractionpatterndefinitions.o \
$(DIRFELIX)$(PRECISION)in.o \
$(DIRFELIX)$(PRECISION)scatteringfactors.o \
$(DIRFELIX)$(PRECISION)writeoutput.o \
$(DIRFELIX)$(PRECISION)errorchecks.o \
$(DIRFELIX)$(PRECISION)simplex.o \
$(DIRFELIX)$(PRECISION)readcif.o \
$(DIRFELIX)$(PRECISION)symmetry.o
#$(DIRFELIX)$(PRECISION)obselete.o

QUADPACK= \
  $(DIRQUADPACK)$(PRECISION)quadpack_double.o

CIFTBX=\
$(DIRCIFTBX)$(PRECISION)ciftbx.o $(DIRCIFTBX)$(PRECISION)hash_funcs.o

SAMPLES=$(MAIN).o

#STARTDIR=$(PWD)
#MYSTARTDIR=$(STARTDIR)

# Linux
#LIBS=-ljadamilu  -llapack -lblas -lm -lc 

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


