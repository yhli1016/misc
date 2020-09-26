set_mod $1 pkg $HOME/soft/libxc-4.3.4
set_mod $1 pkg $HOME/soft/libvdwxc-0.4.0
set_mod $1 pkg $HOME/soft/elpa-2019.11.001
set_mod $1 pkg $HOME/soft/fftw-3.3.8

if [[ "$1" == "add" ]]; then
    export OMP_NUM_THREADS=1
else
    unset OMP_NUM_THREADS
fi
