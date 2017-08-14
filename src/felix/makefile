 
include ../makefiles/platforms.mk
include ../user.mk
include ../makefiles/default.mk


main: $(FELIX)


#$(PRECISION)*.o: ../include/fftw03.f03
#$(PRECISION)*.o: ../include/ciftbx-f90.cmn


print:
	a2ps -E --pro=color --chars-per-line=100 -o print.ps *.f90
# JR convert to pdf does not work correctly
# so I convert manually using ubunutu document viewer
