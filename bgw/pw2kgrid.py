#! /usr/bin/env python

# Last modified on 2015-12-20
#
# This program generates input for kgrid.x from the input of pw.x.
#
# Usage: pw2kgrid.py in out nx ny nz dx dy dz qx qy qz gx gy gz
#   in is the input of pwscf,
#   qx/ny/nz is the size of MP kgrid,
#   dx/dy/dz is the shift of MO kgrid
#   qx/qy/qz is the shift for epsilon and absorption calculation
#   gx/gy/gz is the size of FFT grid determined from pwscf calculation or gsphere.py
#
# Notes:
# (1) Only scf.in, nscf.in and bands.in are acceptable.
# (2) Lattice vectors must be specified manually, i.e. ibrav = 0.
# (3) Atomic positions must be fractional coordinates, i.e. ATOMIC_POSITIONS {crystal}.
# (4) Specifying FFT grid size is mandatory, otherwise the following GW/BSE calculations
#     are likely to CRASH due to symmetry reasons.

import sys

def msg( message ):
    print "\a\n" + message + "\n"

def readinp( inpname ):
    inpcont = []
    inpfile = open( inpname )
    rawcont = inpfile.readlines()
    inpfile.close()
    for line in rawcont:
        if line.find("#") == -1 and line.find("!") == -1:
            inpcont.append(line)
    return inpcont

def checkfrac( inpcont ):
    lfrac = False
    for line in inpcont:
        if line.find("ATOMIC_POSITIONS") != -1 and line.find("crystal") != -1:
            lfrac = True
    return lfrac

def checkibrav( inpcont ):
    lbrav = False
    for line in inpcont:
        if line.find("ibrav") != -1:
            sline = line.split()
            if int(sline[2]) == 0:
                lbrav = True
    return lbrav

def extnat( inpcont ):
    nat = 0
    for line in inpcont:
        if line.find("nat") != -1:
            nat = int( line.split()[2] )
    return nat
    
def extspecies( inpcont ):
    species = []
    for line in inpcont:
        if line.find("UPF") != -1:
            species.append( line.split()[0] )
    return species

def extvector( inpcont ):
    vector = []
    for line in inpcont:
        if line.find("CELL_PARAMETERS") != -1:
            indstart = inpcont.index(line) + 1
            for i in range( indstart, indstart+3 ):
                sline = inpcont[i].split()
                vector.append( [ float(sline[0]), float(sline[1]), float(sline[2]) ] )
    return vector
    
def extpos( inpcont ):
    atomlist = []
    atompos  = []
    indstart = 0
    indend   = 0
    for line in inpcont:
        if line.find("ATOMIC_POSITIONS") != -1:
            indstart = inpcont.index(line) + 1
    nat = extnat( inpcont )
    indend = indstart + nat
    for i in range( indstart, indend ):
        sline = inpcont[i].split()
        atomlist.append(sline[0])
        atompos.append([ float(sline[1]), float(sline[2]), float(sline[3]) ])
    return atomlist, atompos

def frac2cart( FracCoord, BaseVec ):
    CartCoord = []
    for i in FracCoord:
        x0 = i[0]
        y0 = i[1]
        z0 = i[2]
        x1 = x0 * BaseVec[0][0] + y0 * BaseVec[1][0] + z0 * BaseVec[2][0]
        y1 = x0 * BaseVec[0][1] + y0 * BaseVec[1][1] + z0 * BaseVec[2][1]
        z1 = x0 * BaseVec[0][2] + y0 * BaseVec[1][2] + z0 * BaseVec[2][2]
        CartCoord.append( [ x1, y1, z1 ] )
    return CartCoord

def main( argv = None ):
    
    try:
        # parse command-line parameters
        argv = sys.argv
        inpname = argv[1]
        outname = argv[2]
        nx = int( argv[3] )
        ny = int( argv[4] )
        nz = int( argv[5] )
        dx = float( argv[6] )
        dy = float( argv[7] )
        dz = float( argv[8] )
        qx = float( argv[9] )
        qy = float( argv[10] )
        qz = float( argv[11] )
        gx = int( argv[12] )
        gy = int( argv[13] )
        gz = int( argv[14] )
    except:
        msg("Usage: pw2kgrid.py in out nx ny nz dx dy dz qx qy qz gx gy gz")
        return 1
        
    try:
        # read PWSCF input
        inpcont = readinp( inpname )
    except:
        msg( "ERROR: cannot open " + inpname )
        return 2
    
    try:
        # check PWSCF input
        if checkfrac( inpcont ) == False:
            msg( "ERROR: Atomic positions are not fractional!" )
            return 3
        elif checkibrav( inpcont ) == False:
            msg( "ERROR: ibrav != 0" )
            return 4
        else:
            # extract other information from PWSCF input
            nat = extnat( inpcont )
            species = extspecies( inpcont )
            basevec = extvector( inpcont )
            atomlist, atompos = extpos( inpcont )
            cartpos = frac2cart( atompos, basevec )
            try:
                # output
                ofile = open( outname, "w" )
            except:
                msg("ERROR: cannot open " + outname)
                return 5
            # write kgrid
            ofile.write( "%8d%8d%8d\n" % ( nx, ny, nz ) )
            ofile.write( "%8.3f%8.3f%8.3f\n" % ( dx, dy, dz ) )
            ofile.write( "%8.3f%8.3f%8.3f\n" % ( qx, qy, qz ) )
            ofile.write( "\n" )
            # write basevec
            for i in range(3):
                ofile.write( "%16.9f%16.9f%16.9f\n" % ( basevec[i][0], basevec[i][1], basevec[i][2] ) )
            ofile.write( "\n" )
            # write atoms
            ofile.write( "%4d\n" % ( nat ) )
            for i in range( nat ):
                ofile.write( "%4d%16.9f%16.9f%16.9f\n" % ( species.index(atomlist[i])+1, \
                            cartpos[i][0], cartpos[i][1], cartpos[i][2] ) )
            ofile.write( "\n" )
            # write FFT grid and symmetry options
            ofile.write( "%4d%4d%4d\n" % ( gx, gy, gz ) )
            ofile.write( "%8s\n%8s\n" % ( ".false.", ".false." ) )
            ofile.close()
    except:
        msg("ERROR: illegal content of " + inpname )
        return 6
    
    return 0

if __name__ == "__main__":
   sys.exit(main())
