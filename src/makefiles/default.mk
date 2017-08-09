
# default EXE name
MAIN=felixrefine

# -----------------------------------------------------

FELIX=\
$(DIRFELIX)$(PRECISION)gmodules.o \
$(DIRFELIX)$(PRECISION)smodules.o \
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
$(DIRFELIX)$(PRECISION)message_mod.o \
$(DIRFELIX)$(PRECISION)simplex.o \
$(DIRFELIX)$(PRECISION)readcif.o \
$(DIRFELIX)$(PRECISION)symmetry.o
#$(DIRFELIX)$(PRECISION)obselete.o

QUADPACK= \
  $(DIRQUADPACK)$(PRECISION)quadpack_double.o

CIFTBX=\
$(DIRCIFTBX)$(PRECISION)ciftbx.o $(DIRCIFTBX)hash_funcs.o

SAMPLES=$(MAIN).o

STARTDIR=$(PWD)
MYSTARTDIR=$(STARTDIR)

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
LIBDIR=$(MYSTARTDIR)/lib/$(MYPLATFORM)


.SUFFIXES: .c .f .f90 .F .o .a
.DEFAULT: main

#%.o: %.c
#	$(CC)  $(CCFLAGS)  -I$(INCDIR)  $(FORTRANNAMES) $(ARITHMETIC) $(LONGINTEGER) -c -o $(PRECISION)$@ $<
#
#
#%.o: %.F
#	$(FCOMPILE)

$(PRECISION)%.o: %.f %.f90
	@ echo --- --- compiling fortran code
	$(FF)  $(FFFLAGS)  -c -o $(PRECISION)$@ $<


$(PRECISION)%.o: %.c
	$(CC)  $(CCFLAGS)  -I$(INCDIR)  $(FORTRANNAMES) $(ARITHMETIC) $(LONGINTEGER) -c -o $@ $<


$(PRECISION)%.o: %.F
	$(FCOMPILE)



