FC = gfortran
FCFLAGS = -O2

.PHONY: default clean
default: plotband.x
clean:
	rm plotband.x

plotband.x: plotband.f90
	$(FC) $(FCFLAGS) -o $@ $<