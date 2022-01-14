#! /usr/bin/env python
"""Script for generating input files for BerkeleyGW."""
import os
import m4tools as m4
import pw2kgrid as p2k


def run_kgrid():
    p2k.main(qe_struct, "wfn.inp",    nk,    dk,    dq,    ng)
    p2k.main(qe_struct, "wfnq.inp",   nkq,   dkq,   dqq,   ng)
    p2k.main(qe_struct, "wfn_fi.inp", nk_fi, dk_fi, dq_fi, ng)
    for prefix in ["wfn", "wfnq", "wfn_fi"]:
        os.system("kgrid.x %s.inp %s.out %s.log" % (prefix, prefix, prefix))


def link_upf():
    prefixes = ["01-scf", "02-wfn", "03-wfnq", "04-wfn_path", "05-wfn_fi"]
    for prefix in prefixes:
        os.system("ln -sf ../share/*.UPF ../%s" % prefix)


def run_m4():
    m4.write_defs(m, outfile="defs.m4")
    args = "-I../share"
    prefixes = ["01-scf", "02-wfn", "03-wfnq", "04-wfn_path", "05-wfn_fi"]
    for prefix in prefixes:
        if prefix == "01-scf":
            m4.run_m4("../%s/scf.m4" % prefix, "../%s/scf.in" % prefix, args)
        else:
            m4.run_m4("../%s/bands.m4" % prefix, "../%s/bands.in" % prefix, args)
            m4.run_m4("../%s/p2b.m4" % prefix, "../%s/p2b.in" % prefix, args)


## kgrid.x

# 02-wfn
nk = 
dk = [0., 0., 0.]
dq = [0., 0., 0.]

# 03-wfnq
nkq = nk
dkq = dk
dqq = 

# 05-wfn_fi
nk_fi = 
dk_fi = dk
dq_fi = dq

# FFT grid
ng = 

# QE structure file
qe_struct = 

## pw.x
m = dict()

# CONTROL
m["PREFIX"] = 

# SYSTEM
m["NTYP"] = 
m["NAT"] = 
m["ECUTWFC"] = 

m["NBND"] = 
m["NBNDQ"] = 
m["NBND_FI"] = 
m["NBND_PATH"] = 

# ATOMIC_POSITIONS
m["POS"] = m4.include(qe_struct, )

# K_POINTS
m["KPT_SCF"] = 
m["KPT_PATH"] = 

## pw2bgw.x
m["XC_MIN"] = 
m["XC_MAX"] = 

for i in range(3):
    j = i + 1
    m["NK%d" % j] = nk[i]
    m["DK%d" % j] = dk[i] + nk[i] * dq[i]
    m["NK%dQ" % j] = nkq[i]
    m["DK%dQ" % j] = dkq[i] + nkq[i] * dqq[i]
    m["NK%d_FI" % j] = nk_fi[i]
    m["DK%d_FI" % j] = dk_fi[i] + nk_fi[i] * dq_fi[i]


## Generate input files
run_kgrid()
link_upf()
run_m4()
