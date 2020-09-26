set_mod $1 pkg $HOME/soft/libxc-4.3.4
set_mod $1 pkg $HOME/soft/libbeef
set_mod $1 pkg $HOME/soft/elpa-2019.11.001
set_mod $1 bin $HOME/qe-6.6/bin

if [[ "$1" == "add" ]]; then
    export OMP_NUM_THREADS=1
else
    unset OMP_NUM_THREADS
fi
