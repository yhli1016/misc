#! /usr/bin/env python

from sys import *

ifnm = argv[1]
ofnm = argv[2]

ifile = open(ifnm, "r")
content = ifile.readlines()
ifile.close()

ofile = open(ofnm, "w")
for line in content:
    sline = line.split()
    kpath = float(sline[0])
    energy = []
    for data in sline[1:]:
        energy.append(float(data))
    energy.sort()

    ofile.write("%12.5f" % kpath)
    for data in energy:
        ofile.write("%15.9f" % data)
    ofile.write("\n")
ofile.close()
