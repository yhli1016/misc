# Compiler settings
CC = gcc
CFLAGS = -O3 -march=native -mtune=native -fopenmp

# Top level pseudo-targets
.PHONY: default all clean
default: all
all: psi.x overlap.x rad.x
clean:
	rm *.o *.x

# Executables
psi.x: rnl.o ylm.o checkinp.o
overlap.x: rnl.o ylm.o
rad.x: rnl.o
%.x: %.o
	$(CC) $(CFLAGS) -lm -o $@ $^

# Objects
# Each library depends on its own header
rnl.o: rnl.h
ylm.o: ylm.h
checkinp.o: checkinp.h

# Objects using libraries
psi.o: rnl.h ylm.h checkinp.h
overlap.o: rnl.h ylm.h
rad.o: rnl.h

# Common rule to produce objects
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
