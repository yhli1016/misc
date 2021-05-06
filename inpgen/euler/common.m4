include(defs.m4)dnl
#!/bin/bash
#BSUB -n NCPU 
#BSUB -W TIME:00
#BSUB -J NAME 
#BSUB -o log 

# 1) set the obname (BSUB -J)
# 2) set scratchroot
# 3) set number of run 
# 4) set ncores / time
# 5) adjust in INCAR  < 8 CPU --> NPAR=1, 16 CPU --> NPAR=2, >= 32 CPU --> NPAR=4  

name=${LSB_JOBNAME}
run=RUN

scratchroot=/cluster/scratch/zhangwenj
scratch=${scratchroot}/${name}
