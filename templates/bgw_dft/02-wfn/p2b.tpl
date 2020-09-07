#
# Input Template for pw2bgw.x
#
# Version: 1.0.0
#
# Notes:
# 
# (1) Complex version is always used (real_or_complex = 2).
#
# (2) This file is used to generate WFN, RHO and vxc.dat.
#
&INPUT_PW2BGW
    prefix           = '<PREFIX>'
    real_or_complex  = 2
    wfng_flag        = .true.
    wfng_file        = 'WFN'
    wfng_kgrid       = .true.
    wfng_nk1         = <NK1>
    wfng_nk2         = <NK2>
    wfng_nk3         = <NK3>
    wfng_dk1         = <DK1>
    wfng_dk2         = <DK2>
    wfng_dk3         = <DK3>
    rhog_flag        = .true.
    rhog_file        = 'RHO'
    vxcg_flag        = .false.
    vxcg_file        = 'VXC'
    vxc_flag         = .true.
    vxc_file         = 'vxc.dat'
    vxc_diag_nmin    = <XC_MIN>
    vxc_diag_nmax    = <XC_MAX>
    vxc_offdiag_nmin = 0
    vxc_offdiag_nmax = 0
/
