#! /usr/bin/env python
"""
This program is similar to filthr.py with -f option. It extracts hij(R) for
given i,j and R and outputs as a matrix for physical inspiration.

Usage: dumphr.py seedname_hr.dat Rn(1:3) imin imax jmin jmax
       where imin/imax and jmin/jmax of the lower and upper bounds of local
       orbital indices of <bra| and |ket> respectively.
"""

import sys


def readhr(hrname):
    hrfile = open(hrname)
    rawhr = hrfile.readlines()
    hrfile.close()
    for line in rawhr[1:]:
        s = line.split()
        if len(s) == 7 and int(s[0]) <= 0:
            nl0 = rawhr.index(line)
            break
    hr = []
    for line in rawhr[nl0:]:
        s = line.split()
        hr.append([int(s[3]), int(s[4]), int(s[0]), int(s[1]), int(s[2]), float(s[5]), float(s[6])])
    return hr


def extrhop_hopterm(hr, hopterm):
    for hri in hr:
        if hri[0:5] == hopterm:
            return hri


def main(hr, r1, r2, r3, imin, imax, jmin, jmax):
    # echo
    atomilist = range(imin, imax+1)
    atomjlist = range(jmin, jmax+1)
    #print(atomilist)
    #print(atomjlist)

    # extract hij 
    hij_real = []
    hij_imag = []
    for atomi in atomilist:
        row_real = []
        row_imag = []
        for atomj in atomjlist:
            hopterm = [atomi, atomj, r1, r2, r3]
            hop = extrhop_hopterm(hr, hopterm)
            row_real.append(hop[5])
            row_imag.append(hop[6])
        hij_real.append(row_real)
        hij_imag.append(row_imag)

    return hij_real, hij_imag


if __name__ == "__main__":
    # parse command-line parameters
    try:
        hrname = sys.argv[1]
        r1, r2, r3 = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        imin, imax = int(sys.argv[5]), int(sys.argv[6])
        jmin, jmax = int(sys.argv[7]), int(sys.argv[8])
    except:
        print("Usage: dumphr.dat seedname_hr.dat Rn(1:3) imin imax jmin jmax")
        sys.exit(-1)

    # load seedname_hr.dat
    try:
        hr = readhr(hrname)
    except:
        print("ERROR: Cannot open " + hrname)
        sys.exit(-1)

    # extract and print matrix element
    hij_real, hij_imag = main(hr, r1, r2, r3, imin, imax, jmin, jmax)
    for row in hij_real:
        for hij in row:
            print("%9.3f" % hij, end="")
        print("")
