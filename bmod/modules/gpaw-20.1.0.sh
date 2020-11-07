set_mod $1 pkg $HOME/soft/libxc-4.3.4
set_mod $1 pkg $HOME/soft/libvdwxc-0.4.0
set_mod rm pkg $HOME/soft/elpa-2018.11.001
set_mod $1 pkg $HOME/soft/elpa-2019.11.001
set_mod $1 pkg $HOME/soft/fftw-3.3.8
reset_env $1 OMP_NUM_THREADS 1
