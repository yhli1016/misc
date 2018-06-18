#! /usr/bin/env python

# Wannier90 requires the outer and frozen energy windows in the input if
# disentanglement is enabled. This program extracts the energy window for
# given band from the output of pw.x.
#
# Usage: engwin.py bands.out iband, where iband is counting from 1

import sys

pwoutfn = sys.argv[1]
iband = int(sys.argv[2])

pwoutfile = open(pwoutfn, "r")
pwout = pwoutfile.readlines()
pwoutfile.close()

# extract line numbers of kpoints
nlindex = []
for line in pwout:
    if line.find("k =") != -1 and line.find("bands") != -1:
        nlindex.append(pwout.index(line))


# extract full band structure
band = []
# loop over kpoints
for nl in nlindex:
    band_k = []
    # loop over bands
    for line in pwout[nl+2:]:
        sline = line.split()
        # check if this line contains data
        if len(sline) > 0:
            for eng in sline:
                band_k.append(float(eng))
        else:
            break
    # add band_k to band
    band.append(band_k)

# extract band structure for given band
band_i = []
for band_k in band:
    band_i.append(band_k[iband-1])

# echo
print min(band_i), max(band_i)
