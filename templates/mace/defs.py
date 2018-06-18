#! /usr/bin/env python
"""Macro definitions."""

import sys


################################################################################
def include(filename, nl=[0,0]):
    try:
        with open(filename, "r") as infile:
            content = infile.readlines()
    except IOError:
        print("ERROR: cannot include '%s'" % filename)
        sys.exit(-1)
    else:
        nl_start = nl[0] if nl[0] != 0 else 1
        nl_end   = nl[1] if nl[1] != 0 else len(content)
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
m["NK1"] = 
m["NK2"] = 
m["NK3"] = 
m["DK1"] = 
m["DK2"] = 
m["DK3"] = 
m["XC_MIN"] = 
m["XC_MAX"] = 

# 03-wfnq
m["DK1Q"] = 
m["DK2Q"] = 
m["DK3Q"] = 

# 05-wfn_fi
m["NK1_FI"] = 
m["NK2_FI"] = 
m["NK3_FI"] = 
