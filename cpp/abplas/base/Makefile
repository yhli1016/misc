# Compiler settings
CXX = g++
CXXFLAGS = -O3 -march=native -mtune=native

# Top level pseudo-targets
.PHONY: default all clean
default: all 
all: test_input.x test_lattice.x test_struct.x
clean:
	rm *.o *.x

# Executables
test_input.x: input.o lattice.o
test_lattice.x: input.o lattice.o
test_struct.x: input.o lattice.o
%.x: %.o 
	$(CXX) $(CXXFLAGS) -lm -o $@ $^

# Objects
lattice.o: lattice.h consts.h
input.o: input.h lattice.h consts.h
test_input.o: input.h
test_lattice.o: lattice.h
test_struct.o: input.h lattice.h

# Common rule to produce objects
%.o: %.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<
