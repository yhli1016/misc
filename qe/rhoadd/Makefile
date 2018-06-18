FC = gfortran
FCFLAGS = -ffree-form -fopenmp -O3 -march=native -mtune=native
FPP = /usr/bin/cpp -ansi
FPPFLAGS = -DOPENMP

SRC = libcube.f90 sumrho.f90 rhoadd.f90
OBJ = $(SRC:.f90=.o)

.PHONY: default clean
default: rhoadd.x
clean:
	rm *.o *.x *.f

rhoadd.x: $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $^
%.o: %.f90
	$(FPP) $(FPPFLAGS) -o $*.f $<
	$(FC) $(FCFLAGS) -c -o $@ $*.f
