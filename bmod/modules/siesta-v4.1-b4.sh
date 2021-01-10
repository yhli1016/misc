siesta=$HOME/soft/dft/siesta-v4.1-b4
set_mod $1 pkg $siesta/Docs/build/flook/0.8.1
set_mod $1 pkg $siesta/Docs/build/hdf5/1.8.21
set_mod $1 pkg $siesta/Docs/build/netcdf/4.7.4
set_mod $1 pkg $siesta/Docs/build/zlib/1.2.11
set_mod $1 bin $siesta/Obj
unset siesta
reset_env $1 OMP_NUM_THREADS 1
