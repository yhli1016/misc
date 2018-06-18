#! /usr/bin/env python

from sys import *
from math import *

ifilenm = argv[1]
ofilenm = argv[2]

ifile = open( ifilenm, "r" )
kpoints = ifile.readlines()
ifile.close()

kpath = []
kpath.append( 0.0 )
for i in range( 1, len( kpoints ) ):
    dkx = float(kpoints[i].split()[0]) - float(kpoints[i-1].split()[0])
    dky = float(kpoints[i].split()[1]) - float(kpoints[i-1].split()[1])
    dkz = float(kpoints[i].split()[2]) - float(kpoints[i-1].split()[2])
    dki = sqrt( dkx**2 + dky**2 + dkz**2 )
    kpath.append( kpath[i-1] + dki )

ofile = open( ofilenm, "w" )
for kpathi in kpath:
    ofile.write( "%12.5f\n" % kpathi )
ofile.close()
