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


m = macro = {}

# Macros can be strings, integers and floats.
m["PREFIX"] = "'test'"
m["NTYP"] = 1
m["NAT"] = 4
m["ECUTWFC"] = 60.5

# A macro may be dependent on other macros.
m["ECUTRHO"] = m["ECUTWFC"] * 4
m["TEST"] = "# TEST"
m["TEST2"] = m["TEST"] + "2"

# String macros may consist of multiple lines.
m["DP1"] = """
    tefield          = .TRUE.
    dipfield         = .TRUE.
"""

m["DP2"] = """
    edir             = 3
    emaxpos          = 0.90
    eopreg           = 0.05
    eamp             = 0.0
"""

m["POS"] = """
ATOMIC_SPECIES
   P  30.97   P.pbe-mt_fhi.UPF
CELL_PARAMETERS {alat}
   3.343472076   0.000000000   0.000000000
   0.000000000   4.675443939   0.000000000
   0.000000000   0.000000000  20.000000000
ATOMIC_POSITIONS {crystal}
   P   0.800035018   0.092073099   0.567952767
   P   0.800156942   0.258150821   0.460485652
   P   0.300083492   0.575667132   0.474444252
   P   0.299987230   0.775400398   0.576765798
"""

# Including a file as macro is also supported.
m["KPT"] = include("kpt.inc")
