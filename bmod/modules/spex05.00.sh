set_mod $1 bin $HOME/soft/fleurMaXR3.1-nowan/build.GNU
set_mod $1 bin+lib $HOME/soft/wannier90-2.1.0
set_mod $1 bin $HOME/soft/spex05.00/bin

if [[ "$1" == "add" ]]; then
    export OMP_NUM_THREADS=1
else
    unset OMP_NUM_THREADS
fi
