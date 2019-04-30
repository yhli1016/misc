#! /usr/bin/env python
"""Script for generating input files for BerkeleyGW."""
import os
import shutil
import glob
import mace
import pw2kgrid


def run_kgrid():
    pw2kgrid.main(qe_struct, "wfn.inp",    nk,    dk,    dq,    ng)
    pw2kgrid.main(qe_struct, "wfnq.inp",   nkq,   dkq,   dqq,   ng)
    pw2kgrid.main(qe_struct, "wfn_fi.inp", nk_fi, dk_fi, dq_fi, ng)
    for prefix in ["wfn", "wfnq", "wfn_fi"]:
        os.system("kgrid.x %s.inp %s.out %s.log" % (prefix, prefix, prefix))


def copy_upf():
    prefixes = ["01-scf", "02-wfn", "03-wfnq", "04-wfn_path", "05-wfn_fi"]
    upfs = glob.glob("*.UPF")
    for prefix in prefixes:
        for upf in upfs:
            shutil.copyfile(upf, "../%s/%s" % (prefix, upf))


def run_mace():
    prefixes = ["01-scf", "02-wfn", "03-wfnq", "04-wfn_path", "05-wfn_fi"]
    for prefix in prefixes:
        if prefix == "01-scf":
            mace.main(m, "../%s/scf.tpl" % prefix, "../%s/scf.in" % prefix)
        else:
            mace.main(m, "../%s/bands.tpl" % prefix, "../%s/bands.in" % prefix)
            mace.main(m, "../%s/p2b.tpl" % prefix, "../%s/p2b.in" % prefix)


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

# Run kgrid.x
run_kgrid()


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
m["POS"] = mace.include(qe_struct, )

# K_POINTS
m["KPT_SCF"] = 
m["KPT"] = mace.include("wfn.out")
m["KPTQ"] = mace.include("wfnq.out")
m["KPT_FI"] = mace.include("wfn_fi.out")
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
copy_upf()
run_mace()
