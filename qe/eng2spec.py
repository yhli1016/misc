#! /usr/bin/env python

# This program converts energy levels in matrix form to the format suitable for
# plotting energy spectrum with gnuplot.
#
# Usage: eng2spec.py infilename outfilename width_data width_space
#        where width_data and width_space apply for one repeating unit in the x
#        direction only.
#
# Example input:
#-17.2861  -8.2078  -8.2078 -17.2861
# -6.5110  -2.5149  -2.5149  -6.5110
# -6.5110  -2.5149  -2.5149  -6.5110
# -6.5110  -2.5149  -2.5149  -6.5110
#
# Corresponding output (with wdata = 2 and wspace = 1)
#
#   0   -17.28610   3    -8.20780   6    -8.20780   9   -17.28610
#   2   -17.28610   5    -8.20780   8    -8.20780  11   -17.28610
#
#   0    -6.51100   3    -2.51490   6    -2.51490   9    -6.51100
#   2    -6.51100   5    -2.51490   8    -2.51490  11    -6.51100
#
#   0    -6.51100   3    -2.51490   6    -2.51490   9    -6.51100
#   2    -6.51100   5    -2.51490   8    -2.51490  11    -6.51100
#
#   0    -6.51100   3    -2.51490   6    -2.51490   9    -6.51100
#   2    -6.51100   5    -2.51490   8    -2.51490  11    -6.51100

import sys

# parse cli-parameters
try:
    ifname = sys.argv[1]
    ofname = sys.argv[2]
    wdata = float(sys.argv[3])
    wspace = float(sys.argv[4])
except:
    print("Usage: eng2spec.py ifname ofname wdata wspace")
    sys.exit(-1)

# read and parse infile
ifile = open(ifname, "r")
rawcontent = ifile.readlines()
ifile.close()
data = []
for line in rawcontent:
    s = line.split()
    row = []
    for f in s:
        row.append(float(f))
    data.append(row)

# write to outfile
wtot = wdata + wspace
ofile = open(ofname, "w")
for row in data:
    # write the coordinates of starting points
    for icol in range(len(row)):
        x = wtot * icol
        ofile.write("%4d%12.5f" % (x, row[icol]))
    ofile.write('\n')
    # write the coordinates of ending points
    for icol in range(len(row)):
        x = wtot * icol + wdata
        ofile.write("%4d%12.5f" % (x, row[icol]))
    ofile.write('\n')
    # a blank line is required for gnuplot
    ofile.write('\n')
ofile.close()
