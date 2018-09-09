#! /usr/bin/env python
"""Macro definitions."""

import sys


################################################################################
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
################################################################################


m = macro = dict()

# CONTROL
m["PREFIX"] = "''"

# SYSTEM
m["NTYP"] = 
m["NAT"] = 
m["ECUTWFC"] = 

m["NBND"] = 
m["NBNDQ"] = 
m["NBND_FI"] = 
m["NBND_PATH"] = 

# ATOMIC_POSITIONS
m["POS"] = include("../share/scf.in", [,])

# K_POINTS
m["KPT_SCF"] = include("../share/scf.in", [,])

m["KPT"] = include("../share/wfn.out")
m["KPTQ"] = include("../share/wfnq.out")
m["KPT_FI"] = include("../share/wfn_fi.out")

m["KPT_PATH"] = """
# add your data
"""

# 02-wfn
# DK = dk + dq * nk
# dk, dq, nk are from run_kgrid.sh
m["NK1"] = 
m["NK2"] = 
m["NK3"] = 
m["DK1"] = 0
m["DK2"] = 0
m["DK3"] = 0
m["XC_MIN"] = 
m["XC_MAX"] = 

# 03-wfnq
# sharing NK with 02-wfn
m["DK1Q"] = 
m["DK2Q"] = 
m["DK3Q"] = 

# 05-wfn_fi
# sharing DK with 02-wfn
m["NK1_FI"] = 
m["NK2_FI"] = 
m["NK3_FI"] = 
