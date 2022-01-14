cmake_minimum_required(VERSION 3.4)

project (psi)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CFLAGS "-O3 -march=native -mtune=native -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "-lm")

add_executable(psi.x psi.c checkinp.c)
add_executable(rad.x rad.c)
add_executable(overlap.x overlap.c)

target_link_libraries(psi.x func)
target_link_libraries(rad.x func)
target_link_libraries(overlap.x func)

add_library(func rnl.c ylm.c)

# Cmake can look up the '#include' statements in C/C++ source files of a target
# for any necessary recompilation of the component and the target. So you do
# not have to care about the header files.
#
# Unfortunately that does not work for FORTRAN. Use SCons instead.
