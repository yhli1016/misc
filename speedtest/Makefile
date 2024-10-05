CXX = g++
CXXFLAGS = -O3 -fopenmp -DWITH_FULL_CSR

TBPLAS_ROOT = ${HOME}/proj/tbplas-cpp
INCFLAGS = -I${TBPLAS_ROOT}/sources -I${TBPLAS_ROOT}/3rdparty
LIBFLAGS = -L${TBPLAS_ROOT}/build/out -ltbpm -ldiag -lbuilder -lbase -ladapter

# Top level pseudo-targets
.PHONY: default all clean
default: all
all: speedtest
clean:
	rm *.o speedtest

# Executables
speedtest: model.o diag.o tbpm.o speedtest.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

# Objects
# Each library depends on its own header
model.o: model.h
diag.o: diag.h
tbpm.o: tbpm.h

# Objects using libraries
speedtest.o: model.h diag.h tbpm.h

# Common rule to produce objects
%.o: %.cpp
	$(CXX) $(INCFLAGS) $(CXXFLAGS) -c -o $@ $<
