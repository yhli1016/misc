"""
Library for summing up the PDOS on given atomic states.

Most of the time you only need to call the main() function, which accepts the
following arguments:
  prefix      # string, prefix in pwscf calculations
  atomid      # list of strings,  holding the indices of given atoms
  atomsym     # string, atomic symbol
  wfcid       # integer, index of atomic state in pseudo-potentials
  lsym        # string, s/p/d/f
  outfile     # string, file name of output file
  efermi      # float, fermi energy in eV
 
Two auxiliary functions extract_atomlist() and extract_atomid() are provided
for setting up atomid.
"""


def readpdos(pdosname):
    with open(pdosname, "r") as ifile:
        pdoscontent = ifile.readlines()
    pdosmat = [[float(f) for f in line.split()] for line in pdoscontent[1:]]
    return pdosmat


def writedat(tpdos, outfile):
    with open(outfile, "w") as ofile:
        for row in tpdos:
            ofile.write("%7.3f" % row[0])
            for f in row[1:]:
                ofile.write("%11.3e" % f)
            ofile.write("\n")


def extract_atomlist(scfname):
    """Extract the list of atoms from the input of Quantum-ESPRESSO"""
    with open(scfname, "r") as scffile:
        scfinp = scffile.readlines()
    for line in scfinp:
        if line.find("ATOMIC_POSITIONS") != -1:
            nl0 = scfinp.index(line)
        if line.find("K_POINTS") != -1:
            nl1 = scfinp.index(line)
    atomlist = [line.split()[0] for line in scfinp[nl0+1:nl1]]
    return atomlist


def extract_atomid(atomlist, atomsym):
    """Extract the id of atoms of given atomic symbol"""
    atomid = [i+1 for i in range(len(atomlist)) if atomlist[i] == atomsym]
    return atomid


def main(prefix, atomid, atomsym, wfcid, lsym, outfile, efermi):
    # Begin with the first pdos file
    iname = "%s.pdos_atm#%s(%s)_wfc#%s(%s)" % (prefix, str(atomid[0]), atomsym,
                                               str(wfcid),  lsym)
    tpdos = readpdos(iname)
    nrow = len(tpdos)
    ncol = len(tpdos[0])
    
    # Then deal with remaining files
    for i in atomid[1:]:
        iname = "%s.pdos_atm#%s(%s)_wfc#%s(%s)" % (prefix, str(i), atomsym,
                                                   str(wfcid),  lsym)
        tpdosi = readpdos(iname)
        for j in range(nrow):
            # The first column (energies) is omitted
            for k in range(1, ncol):
                tpdos[j][k] += tpdosi[j][k]
    
    # Shift efermi to 0 eV
    for row in tpdos:
        row[0] -= efermi

    # Output
    writedat(tpdos, outfile)
