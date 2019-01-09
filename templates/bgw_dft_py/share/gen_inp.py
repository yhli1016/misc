#! /usr/bin/env python
"""Script for generating input files for BerkeleyGW."""

import sys
import os
import mace


def include(filename, nl0=None, nl1=None):
    try:
        with open(filename, "r") as infile:
            content = infile.readlines()
    except IOError:
        print("ERROR: cannot include '%s'" % filename)
        sys.exit(-1)
    else:
        nl_start = nl0 if nl0 is not None else 1
        nl_end   = nl1 if nl1 is not None else len(content)
        longline = "".join(content[(nl_start-1):nl_end])
        return longline


def run_kgrid():
    command = "pw2kgrid.py %s wfn.inp" % qe_struct
    command = command + "".join([" %s" % i for i in nk])
    command = command + "".join([" %s" % i for i in dk])
    command = command + "".join([" %s" % i for i in dq])
    command = command + "".join([" %s" % i for i in ng])
    os.system(command)

    command = "pw2kgrid.py %s wfnq.inp" % qe_struct
    command = command + "".join([" %s" % i for i in nkq])
    command = command + "".join([" %s" % i for i in dkq])
    command = command + "".join([" %s" % i for i in dqq])
    command = command + "".join([" %s" % i for i in ng])
    os.system(command)

    command = "pw2kgrid.py %s wfn_fi.inp" % qe_struct
    command = command + "".join([" %s" % i for i in nk_fi])
    command = command + "".join([" %s" % i for i in dk_fi])
    command = command + "".join([" %s" % i for i in dq_fi])
    command = command + "".join([" %s" % i for i in ng])
    os.system(command)

    os.system("kgrid.x wfn.inp    wfn.out    wfn.log")
    os.system("kgrid.x wfnq.inp   wfnq.out   wfnq.log")
    os.system("kgrid.x wfn_fi.inp wfn_fi.out wfn_fi.log")


def run_mace():
    os.system("cp *.UPF ../01-scf")
    mace.main(macro, "../01-scf/scf.tpl", "../01-scf/scf.in")
    for prefix in ["02-wfn", "03-wfnq", "04-wfn_path", "05-wfn_fi"]:
        os.system("cp *.UPF ../%s" % prefix)
        mace.main(macro, "../%s/bands.tpl" % prefix, "../%s/bands.in" % prefix)
        mace.main(macro, "../%s/p2b.tpl" % prefix, "../%s/p2b.in" % prefix)


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
m = macro = dict()

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
m["POS"] = include(qe_struct, )

# K_POINTS
m["KPT_SCF"] = include(qe_struct, )
m["KPT"] = include("wfn.out")
m["KPTQ"] = include("wfnq.out")
m["KPT_FI"] = include("wfn_fi.out")
m["KPT_PATH"] = """
add your settings
"""

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
run_mace()
