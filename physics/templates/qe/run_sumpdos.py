#! /usr/bin/env python

#prefix = "se"
#atomid = [1,2]
#atomsym = "In"
#wfcid = 3 
#lsym = "d"
#outfile = "Ind.dat"
#efermi = 0.0

import sys
import os
ENV_PATH_SYS = os.environ["PATH"].split(":")
ENV_PATH_USR = [path for path in ENV_PATH_SYS if path.find("home") is not -1]
sys.path.extend(ENV_PATH_USR)
from sumpdos import *

prefix = "in"
# atomid is generated automatically from atomlist and atomsym
atomlist = extract_atomlist("scf.in")
atomsymlist = ["In", "Se"]
wfcidlist = [1, 2, 3, 4]
lsymlist = ["s", "p", "d", "f"]
# outfile is generated automatically from wfcid and lsym
efermi = -2.1510

for atomsym in atomsymlist:
    atomid = extract_atomid(atomlist, atomsym)
    for i in range(len(wfcidlist)):
        wfcid = wfcidlist[i]
        lsym = lsymlist[i]
        outfile = atomsym + str(wfcid) + lsym + ".dat"
        main(prefix, atomid, atomsym, wfcid, lsym, outfile, efermi)
