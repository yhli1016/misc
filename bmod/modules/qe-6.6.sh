set_mod $1 pkg $HOME/soft/libxc-4.3.4
set_mod $1 pkg $HOME/soft/libbeef
set_mod rm pkg $HOME/soft/elpa-2018.11.001
set_mod $1 pkg $HOME/soft/elpa-2019.11.001
set_mod $1 bin $HOME/soft/qe-6.6/bin
reset_env $1 OMP_NUM_THREADS 1
