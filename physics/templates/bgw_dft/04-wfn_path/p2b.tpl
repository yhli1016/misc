#
# Input Template for pw2bgw.x
#
# Version: 1.0.0
#
# Notes:
# 
# (1) Complex version is always used (real_or_complex = 2).
#
# (2) This file is used to generate WFN_fi for inteqp.$flavor.x.
#
# (3) For band structure interpolation, wave function along high-symmetric
#     lines are required. For that purpose, we set wfng_kgrid to false.
#
&INPUT_PW2BGW
    prefix          = '<PREFIX>'
    real_or_complex = 2
    wfng_flag       = .true.
    wfng_file       = 'WFN_fi'
    wfng_kgrid      = .false.
/
