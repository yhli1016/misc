#! /usr/bin/env python

"""
This program should be run in the top directory.
"""

def extract_eng(filename, bandlist):
    # read EIGENVAL
    with open(filename, "r") as f:
        eigenval = f.readlines()
    
    # determine the entry for Gamma point
    nl0 = 0
    for line in eigenval:
        s = line.split()
        if len(s) == 4 and float(s[0]) <= 1 and float(s[1]) <= 1 \
        and float(s[2]) <= 1:
            nl0 = eigenval.index(line) + 1
            break

    #  energies for i in [bmin, bmax]
    eng = []
    for ib in bandlist:
        for line in eigenval[nl0:]:
            s = line.split()
            if int(s[0]) == ib:
                eng.append(float(s[1]))
                break
    return eng

def extract_proj(filename, bandlist, atomlist):
    # open PROCAR
    with open(filename, "r") as f:
        procar = f.readlines()
    
    # determine number of ions
    nion = int(procar[1].split()[-1])
    
    # determine the entries for all bands
    band_entries = [i for i in range(len(procar)) \
                    if procar[i].find("ion") != -1 \
                    and procar[i].find("tot") != -1]
                    
    # extract projection
    proj = []
    for ib in bandlist:
        proj_band = 0
        
        # loop over all lines of this band to collect contribution of given
        # atoms
        nl0 = band_entries[ib-1] + 1
        nl1 = nl0 + nion
        for line in procar[nl0:nl1]:
            s = line.split()
            if int(s[0]) in atomlist:
                proj_band += float(s[-1])
        # normalize proj_band
        proj_band /= float(procar[nl1].split()[-1])
        proj.append(proj_band)
    
    return proj
    
def main():
    postmin = 1001
    postmax = 1001
    bandlist = [1,2,3,4]
    atomlist = [1,2,3,4]
    
    data = []
    for post in range(postmin, postmax+1):
        dirname = str("%04d" % post)
        eng = extract_eng(dirname + "/EIGENVAL", bandlist)
        proj = extract_proj(dirname + "/PROCAR", bandlist, atomlist)
        data.append([post, eng, proj])
    
    with open("eng.dat", "w") as ofile:
        for row in data:
            ofile.write("%6d" % row[0])
            for f in row[1]:
                ofile.write("%16.6f" % f)
            ofile.write("\n")
    
    with open("proj.dat", "w") as ofile:
        for row in data:
            ofile.write("%6d" % row[0])
            for f in row[2]:
                ofile.write("%7.3f" % f)
            ofile.write("\n")
    
    with open("max.dat", "w") as ofile:
        for row in data:
            maxband = row[2].index(max(row[2]))
            ofile.write("%6d" % row[0])
            ofile.write("%6d" % bandlist[maxband])
            ofile.write("%16.6f" % row[1][maxband])
            ofile.write("%7.3f\n" % row[2][maxband])
    
if __name__ == "__main__":
    main()
