#! /usr/bin/env python

import sys


# Parse command-line parameters
if len(sys.argv) == 4:
    dos_filename = sys.argv[1]
    band_filename = sys.argv[2]
    out_filename = sys.argv[3]
else:
    print("Usage: analydos_pp.py dosfn bndfn outfn")
    sys.exit()

# Load and parse the output of analy*dos.x
with open(dos_filename, "r") as dos_file:
    dos_raw = dos_file.readlines()[1:]
if len(dos_raw[0].split()) != 3:
    raise IOError("The format of dos file is wrong.")
dos = [[int(line.split()[0]), int(line.split()[1]), float(line.split()[2])]
        for line in dos_raw]
dos.sort(key=lambda x: x[1])

# Load and parse band structure
with open(band_filename, "r") as band_file:
    band_raw = band_file.readlines()
kpath = [float(line.split()[0]) for line in band_raw]
energy = [[float(f) for f in line.split()[1:]] for line in band_raw]

# Check consistency between dos and band structure
ikmin_dos = min([row[0] for row in dos])
ikmax_dos = max([row[0] for row in dos])
ibmin_dos = min([row[1] for row in dos])
ibmax_dos = max([row[1] for row in dos])
nk_band = len(kpath)
nb_band = len(energy[0])
if ikmin_dos > 1 or ikmax_dos < nk_band:
    print("Warning: There are less kpoints in dos than in band structure.")
if ibmin_dos > 1 or ibmax_dos < nb_band:
    print("Warning: There are less bands in dos than in band structure.")

# Ouput
with open(out_filename, "w") as out_file:
    for ib in range(nb_band):
        for ik in range(nk_band):
            x = kpath[ik]
            y = energy[ik][ib]
            z = 0.0
            for row in dos:
                if (ik == row[0] - 1) and (ib == row[1] - 1):
                    z = row[2]
                    dos.remove(row)
                    break
            out_file.write("%12.6f%12.6f%12.6f\n" % (x, y, z))
        out_file.write("\n")
