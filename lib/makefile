#/usr/bin/bash

include ../src/makefiles/platforms.mk
include ../src/user.mk
include ../src/makefiles/default.mk

DIRCIFTBX=../src/ciftbx/
DIRQUADPACK=../src/quadpack/
DIRFELIX=../src/felix/

LIBFELIX=$(MYPLATFORM)/lib$(PRECISION)felix.a
LIBQUADPACK=$(MYPLATFORM)/lib$(PRECISION)quadpack.a
LIBCIFTBX=$(MYPLATFORM)/lib$(PRECISION)ciftbx.a

ourlibs: $(LIBFELIX) $(LIBQUADPACK) $(LIBCIFTBX)


$(LIBFELIX): $(FELIX) 
	mkdir -p $(MYPLATFORM)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

$(LIBQUADPACK):  $(QUADPACK)
	mkdir -p $(MYPLATFORM)
	$(AR) $(ARFLAGS) $@ $<
	$(RL) $(RLFLAGS) $@

$(LIBCIFTBX):  $(CIFTBX)
	mkdir -p $(MYPLATFORM)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

