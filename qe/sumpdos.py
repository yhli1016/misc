#! /usr/bin/env python

"""
This program sums up the PDOS on given atomic states.

This program has two working modes.

In 'program' mode it requires an input file named inp.py,
which must contain the following variables:
  prefix      # prefix in pwscf calculations
  atomid      # list holding the indices of given atoms
  atomsym     # atomic symbol
  wfcid       # index of atomic state in pseudo-potentials
  lsym        # s/p/d/f
  outfile     # file name of output file
  efermi      # fermi energy in eV

In 'library' mode, these variables are passed to 'main' as parameters
directly. There are also two auxiliary functions 'extract_atomlist'
and 'extract_atomid' for this mode.
"""

import sys

def readpdos(pdosname):
    with open(pdosname, "r") as ifile:
        pdoscontent = ifile.readlines()
    # the first line (comment) is omitted
    pdosmat = [[float(f) for f in line.split()] for line in pdoscontent[1:]]
    return pdosmat

def writedat(tpdos, outfile):
    with open(outfile, "w") as ofile:
        for row in tpdos:
            ofile.write("%7.3f" % row[0])
            for f in row[1:]:
                ofile.write("%11.3e" % f)
            ofile.write("\n")

def main(prefix, atomid, atomsym, wfcid, lsym, outfile, efermi):
    # begin with the first pdos file
    iname = prefix + ".pdos_atm#" + str(atomid[0]) + "(" + atomsym + ")_wfc#" + str(wfcid) + "(" + lsym + ")"
    tpdos = readpdos(iname)
    nrow = len(tpdos)
    ncol = len(tpdos[0])
    
    # then deal with remaining files
    for i in atomid[1:]:
        iname = prefix + ".pdos_atm#" + str(i) + "(" + atomsym + ")_wfc#" + str(wfcid) + "(" + lsym + ")"
        tpdosi = readpdos(iname)
        for j in range(nrow):
            # the first column (energies) is omitted
            for k in range(1, ncol):
                tpdos[j][k] += tpdosi[j][k]
    
    # shift efermi to 0 eV
    for row in tpdos:
        row[0] -= efermi

    # output
    writedat(tpdos, outfile)
    return 0

# auxiliary functions for library mode
# Extract the list of atoms from the input of Quantum-ESPRESSO
def extract_atomlist(scfname):
    with open(scfname, "r") as scffile:
        scfinp = scffile.readlines()
    for line in scfinp:
        if line.find("ATOMIC_POSITIONS") != -1:
            nl0 = scfinp.index(line)
        if line.find("K_POINTS") != -1:
            nl1 = scfinp.index(line)
    atomlist = [line.split()[0] for line in scfinp[nl0+1:nl1]]
    return atomlist

# Extract the id of atoms of given atomic symbol
def extract_atomid(atomlist, atomsym):
    atomid = [i+1 for i in range(len(atomlist)) if atomlist[i] == atomsym]
    return atomid

# program mode
if __name__ == "__main__":
    sys.path.append(".")
    from inp import *
    sys.exit(main(prefix, atomid, atomsym, wfcid, lsym, outfile, efermi))
