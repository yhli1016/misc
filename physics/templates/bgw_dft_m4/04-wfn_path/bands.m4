include(defs.m4)
#
# Input Template for Calculation 'bands'
#
# Version: 1.0.0
#
# Notes:
# 
# (1) This file applies for periodic systems only.
#
# (2) The built-in Bravis lattice types are rather unconvenient.
#     So we set ibrav to 0, celldm(1) to 1.8897, and specify the
#     coordinates of lattice vectors in ANGSTROMS under card
#     'CELL_PARAMETERS {alat}'.
#
# (3) Change K_POINTS to 'crystal_b' if you want to calculate 
#     band structure.
#
&CONTROL
    prefix           = 'PREFIX'
    calculation      = 'bands'
    restart_mode     = 'from_scratch'
    pseudo_dir       = './'
    outdir           = './'
    wfcdir           = './'
    verbosity        = 'high'
    wf_collect       = .TRUE.
/
&SYSTEM
    ibrav            = 0
    celldm(1)        = 1.889726125
    ntyp             = NTYP
    nat              = NAT
    ecutwfc          = ECUTWFC
    nbnd             = NBND_PATH
/
&ELECTRONS
    electron_maxstep = 250
    conv_thr         = 1.0d-10
    mixing_mode      = 'plain'
    mixing_beta      = 0.7
    mixing_ndim      = 8
    diagonalization  = 'david'
    diago_full_acc   = .TRUE.
    startingwfc      = 'random'
/
POS
KPT_PATH
