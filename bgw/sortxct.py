#! /usr/bin/env python

# This program reads the output of summarize_eigenvectors.$flavor.x, parses
# and sorts the transition weights for each exciton, then writes to file.
#
# usage: sortxct.py xctname sxctname

import sys

def extrStateIndex(xctweight):
    stateindex = []
    for line in xctweight:
        if line.find("Special analysis for state") != -1:
            stateindex.append(xctweight.index(line))
    return stateindex

def extrNVC(xctweight):
    for line in xctweight:
        if line.find("ns, nv, nc, nk") != -1:
            nv = int(line.split()[6])
            nc = int(line.split()[7])
            return nv * nc

def procState(xctweight, sxctname):
    stateindex = extrStateIndex(xctweight)
    nvc = extrNVC(xctweight)

    # open file for output
    sxctname = sys.argv[2]
    sxctfile = open(sxctname, "w")

    # loop over all states
    for nl0 in stateindex:
        # extract and sort weights
        unsorted_weights = []
        for line in xctweight[(nl0+2):(nl0+2+nvc)]:
            s = line.split()
            unsorted_weights.append([int(s[0]), int(s[1]), float(s[2]), float(s[3]), int(s[4])])
        sorted_weights = sorted(unsorted_weights, key=lambda x: x[3], reverse=True)

        # write to file
        sxctfile.write(xctweight[nl0])
        sxctfile.write(xctweight[nl0+1])
        for w in sorted_weights:
            sxctfile.write("%5d%5d%10.5f%10.5f%10d\n" % (w[0], w[1], w[2], w[3], w[4]))
        sxctfile.write("\n\n")

    # close file
    sxctfile.close()

def main():
    xctname = sys.argv[1]
    sxctname = sys.argv[2]
    xctfile = open(xctname)
    xctweight = xctfile.readlines()
    xctfile.close()
    procState(xctweight, sxctname)

if __name__ == "__main__":
    main()
