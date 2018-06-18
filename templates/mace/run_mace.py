#! /usr/bin/env python
"""Script generating files using mace as a module."""

import sys
sys.path.append("/home/yhli/soft/bin")

import mace


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


m = macro = {}

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
m["POS"] = include("scf.in", [,])

# K_POINTS
m["KPT_SCF"] = include("scf.in", [,])

m["KPT"] = include("wfn.out")
m["KPTQ"] = include("wfnq.out")
m["KPT_FI"] = include("wfn_fi.out")

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

################################################################################
for pref in ["01-scf", "02-wfn", "03-wfnq", "04-wfn_path", "05-wfn_fi"]:
    workdir = "../" + pref + "/"
    if pref == "01-scf":
        mace.main(macro, workdir+"scf.tpl", workdir+"scf.in")
    else:
        mace.main(macro, workdir+"bands.tpl", workdir+"bands.in")
        mace.main(macro, workdir+"p2b.tpl", workdir+"p2b.in")
