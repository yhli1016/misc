#
# Input Template for pw2bgw.x
#
# Version: 1.0.0
#
# Notes:
# 
# (1) Complex version is always used (real_or_complex = 2).
#
# (2) This file is used to generate WFNq.
#
# (3) There are no flags related to q-shift, so add it to wfng_dki.
#     That is to say, wfng_dki = dki + nki * qki
#
&INPUT_PW2BGW
    prefix          = '<PREFIX>'
    real_or_complex = 2
    wfng_flag       = .true.
    wfng_file       = 'WFNq'
    wfng_kgrid      = .true.
    wfng_nk1        = <NK1Q>
    wfng_nk2        = <NK2Q>
    wfng_nk3        = <NK3Q>
    wfng_dk1        = <DK1Q>
    wfng_dk2        = <DK2Q>
    wfng_dk3        = <DK3Q>
/
