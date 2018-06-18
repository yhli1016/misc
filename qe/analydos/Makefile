# Compiler settings
FC = gfortran
FCFLAGS = -O2

# Top level pseudo-targets
.PHONY: default all clean
default: all
all: analyldos.x analypdos.x
clean:
	rm *.o *.mod *.x

# Executables
analyldos.x: sumldos.o shareddata.o
analypdos.x: sumpdos.o shareddata.o
%.x: %.o
	$(FC) $(FCFLAGS) -o $@ $^

# Objects
analyldos.o analypdos.o sumldos.o sumpdos.o: shareddata.o
%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
