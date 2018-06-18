#! /usr/bin/env python

def shift_eng(filename, dE):
    with open(filename, "r") as f:
        dat_raw = f.readlines()
    dat = [[float(s) for s in line.split()] for line in dat_raw]
    #
    for row in dat:
        row[0] += dE
    #
    with open("s" + filename, "w") as f:
        for row in dat:
            f.write("%7.3f" % row[0])
            for s in row[1:]:
                f.write("%11.3e" % s)
            f.write("\n")

dE = 0.4196
for species in ["In", "Se"]:
    for wfc in ["1s", "2p", "3d", "4f"]:
        filename = species + wfc + ".dat"
        shift_eng(filename, dE)
