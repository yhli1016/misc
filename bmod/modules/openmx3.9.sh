set_mod $1 pkg $HOME/soft/fftw-3.3.8
set_mod $1 bin $HOME/soft/openmx3.9/source

if [[ "$1" == "add" ]]; then
    export OMP_NUM_THREADS=1
else
    unset OMP_NUM_THREADS
fi
