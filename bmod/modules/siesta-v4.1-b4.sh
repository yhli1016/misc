# Shared commands for loading and unloading
siesta=$HOME/soft/siesta-v4.1-b4
set_mod $1 pkg $siesta/Docs/build/flook/0.8.1
set_mod $1 pkg $siesta/Docs/build/hdf5/1.8.21
set_mod $1 pkg $siesta/Docs/build/netcdf/4.7.4
set_mod $1 pkg $siesta/Docs/build/zlib/1.2.11
set_mod $1 pkg $HOME/soft/elpa-2019.11.001
set_mod $1 bin $siesta/Obj
unset siesta

# Addtional settings
if [[ "$1" == "add" ]]; then
    export OMP_NUM_THREADS=1
else
    unset OMP_NUM_THREADS
fi
