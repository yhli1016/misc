#
# Input Template for pw2bgw.x
#
# Version: 1.0.0
#
# Notes:
# 
# (1) Complex version is always used (real_or_complex = 2).
#
# (2) This file is used to generate WFN_fi.
#
&INPUT_PW2BGW
    prefix          = <PREFIX>
    real_or_complex = 2
    wfng_flag       = .true.
    wfng_file       = 'WFN_fi'
    wfng_kgrid      = .true.
    wfng_nk1        = <NK1_FI>
    wfng_nk2        = <NK2_FI>
    wfng_nk3        = <NK3_FI>
    wfng_dk1        = <DK1>
    wfng_dk2        = <DK2>
    wfng_dk3        = <DK3>
/
