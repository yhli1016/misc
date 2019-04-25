#! /usr/bin/env python
from sumpdos import main, extract_atomlist, extract_atomid


prefix = "0"
# atomid is generated automatically from atomlist and atomsym
atomlist = extract_atomlist("scf.in")
atomsymlist = ["P"]
wfcidlist = [1, 2, 3, 4]
lsymlist = ["s", "p", "d", "f"]
# outfile is generated automatically from wfcid and lsym
efermi = -2.0177

for atomsym in atomsymlist:
    atomid = extract_atomid(atomlist, atomsym)
    for i in range(len(wfcidlist)):
        wfcid = wfcidlist[i]
        lsym = lsymlist[i]
        outfile = "%s%s%s.dat" % (atomsym, str(wfcid), lsym)
        main(prefix, atomid, atomsym, wfcid, lsym, outfile, efermi)
