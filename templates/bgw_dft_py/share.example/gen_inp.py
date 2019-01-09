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
nk = [9, 9, 1]
dk = [0., 0., 0.]
dq = [0., 0., 0.]

# 03-wfnq
nkq = nk
dkq = dk
dqq = [0.001, 0., 0.]

# 05-wfn_fi
nk_fi = [18, 18, 1]
dk_fi = dk
dq_fi = dq

# FFT grid
ng = [45, 45, 216]

# QE structure file
qe_struct = "scf.in"

# Run kgrid.x
run_kgrid()


## pw.x
m = macro = dict()

# CONTROL
m["PREFIX"] = "0"

# SYSTEM
m["NTYP"] = 2
m["NAT"] = 4
m["ECUTWFC"] = 70

m["NBND"] = 1301
m["NBNDQ"] = 9
m["NBND_FI"] = 20
m["NBND_PATH"] = 20

# ATOMIC_POSITIONS
m["POS"] = include(qe_struct, 47, 58)

# K_POINTS
m["KPT_SCF"] = include(qe_struct, 60, 61)
m["KPT"] = include("wfn.out")
m["KPTQ"] = include("wfnq.out")
m["KPT_FI"] = include("wfn_fi.out")
m["KPT_PATH"] = """
K_POINTS {crystal_b}
   4
   0.000000000   0.000000000   0.000000000  40
   0.500000000   0.000000000   0.000000000  40
   0.666666667   0.333333333   0.000000000  40
   0.000000000   0.000000000   0.000000000  40
"""

## pw2bgw.x
m["XC_MIN"] = 1
m["XC_MAX"] = 40

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
