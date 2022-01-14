# Compiler settings
FC = gfortran
FPPFLAGS = -cpp -DQE6
FCFLAGS = -O2

# Top level pseudo-targets
.PHONY: default all clean
default: all
all: analyldos.x analypdos.x
clean:
	rm *.o *.mod *.x

# Executables
%.x: %.o sumdos.o projection.o
	$(FC) $(FPPFLAGS) $(FCFLAGS) -o $@ $^

# Objects
analyldos.o analypdos.o: sumdos.o projection.o
sumdos.o: projection.o
%.o: %.f90
	$(FC) $(FPPFLAGS) $(FCFLAGS) -c -o $@ $<
