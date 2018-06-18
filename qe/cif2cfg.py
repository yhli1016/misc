#! /usr/bin/env python

# Last modified on 2015-10-13
#
# Max atoms supported: 999
# Max lattice length supported: 9.99 nm, 99.9 angstrom
# Otherwise output would be wrong
#
# This program reads in cif file, shifts coordinates of atoms according to
# specified parameters, then outputs to cfg file together with information
# of crystal cell.
#
# Usage: cif2cfg.py cif cfg 0/1/2/3 0.0(anyfloat)
#   where 0/1/2/3 stands for dimension of the system which has to be consist
#   with truncation mode in MBPT calculations, and anyfloat stands for the
#   geometery center of the system: ( f, f, f ) for 0D, ( f, f, Z ) for 1D,
#   ( X, Y, f ) for 2D

import sys
from math import sin, cos, sqrt, fabs

def msg( message ):
    print "\a\n" + message + "\n"

def GetLattConst( CifContent ):
    LattConst = []
    for line in CifContent:
        if line.find("_cell_length_a") != -1:
            LattConst.append(float(line.split()[1]))
            break
    for line in CifContent:
        if line.find("_cell_length_b") != -1:
            LattConst.append(float(line.split()[1]))
            break
    for line in CifContent:
        if line.find("_cell_length_c") != -1:
            LattConst.append(float(line.split()[1]))
            break
    for line in CifContent:
        if line.find("_cell_angle_alpha") != -1:
            LattConst.append(float(line.split()[1]))
            break
    for line in CifContent:
        if line.find("_cell_angle_beta") != -1:
            LattConst.append(float(line.split()[1]))
            break
    for line in CifContent:
        if line.find("_cell_angle_gamma") != -1:
            LattConst.append(float(line.split()[1]))
            break
    return LattConst

def GetAtomCoord( CifContent ):
    AtomList  = []
    AtomCoord = []
    for line in CifContent:
        if line.find("Uiso") != -1:
            s = line.split()
            AtomList.append(s[1])
            AtomCoord.append([float(s[2]),float(s[3]),float(s[4])])
    return AtomList, AtomCoord

def CalcGeoCenter( AtomCoord ):
    CoordSum = [ 0.0, 0.0, 0.0 ]
    for i in AtomCoord:
        CoordSum[0] = CoordSum[0] + i[0]
        CoordSum[1] = CoordSum[1] + i[1]
        CoordSum[2] = CoordSum[2] + i[2]
    N = len(AtomCoord)
    GeoCenter = [ CoordSum[0]/N, CoordSum[1]/N, CoordSum[2]/N ]
    return GeoCenter

def ShiftCoord( AtomCoord, Dimension, ShiftStd ):
    GeoCenter = CalcGeoCenter(AtomCoord)
    ShiftVec  = [ 0.0, 0.0, 0.0 ]
    if Dimension == 0:
        ShiftVec = [ShiftStd-GeoCenter[0],ShiftStd-GeoCenter[1],ShiftStd-GeoCenter[2]]
    else:
        if Dimension == 1:
            ShiftVec = [ShiftStd-GeoCenter[0],ShiftStd-GeoCenter[1],0.0]
        else:
            if Dimension == 2:
                ShiftVec = [0.0,0.0,ShiftStd-GeoCenter[2]]
            else:
                ShiftVec = [0.0,0.0,0.0]
    ShiftedCoord = []
    for i in AtomCoord:
        ShiftedCoord.append([i[0]+ShiftVec[0],i[1]+ShiftVec[1],i[2]+ShiftVec[2]])
    return ShiftedCoord

def CalcLattVec( LattConst ):
    deg2rad = 1.74532925199433e-2
    a = LattConst[0]
    b = LattConst[1]
    c = LattConst[2]
    alpha = LattConst[3] * deg2rad
    beta  = LattConst[4] * deg2rad
    gamma = LattConst[5] * deg2rad
    LattVec = []
    LattVec.append([ a, 0.0, 0.0 ])
    LattVec.append([ b*cos(gamma), b*sin(gamma), 0.0 ])
    cx = cos(beta)
    cy = (cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
    #cz = sqrt(1.0 - cx**2 - cy**2)
    cz = sqrt( 1 + 2*cos(alpha)*cos(beta)*cos(gamma) - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 ) / sin(gamma)
    LattVec.append([ c*cx, c*cy, c*cz ])
    # in case where componets of vectors are close to zero
    for i in range(0,3):
        for j in range(0,3):
            if fabs(LattVec[i][j]) < 1.0e-13:
                LattVec[i][j] = fabs(LattVec[i][j])
    return LattVec
    
def GetTypeInfo( AtomList ):
    TypeList = []
    NumList  = []
    TypeList.append(AtomList[0])
    for i in AtomList:
        if TypeList.count(i) == 0:
            TypeList.append(i)
    for i in TypeList:
        NumList.append(AtomList.count(i))
    return TypeList, NumList

def Frac2Cart( FracCoord, BaseVec ):
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

def main(argv = None):

    try:
        # parse cli parameters
        argv      = sys.argv
        CifName   = argv[1]
        CfgName   = argv[2]
        Dimension = int(argv[3])
        ShiftStd  = float(argv[4])
    except:
        msg("Usage: cif2cfg.py cif cfg 0/1/2/3 0.0(anyfloat)")
        return 1
    
    # convertion coefficient
    ang2bohr = 1.8897261246
    deg2rad  = 1.74532925199433e-2
    
    try:
        # read cif file
        CifFile = open( CifName, "r" )
        CifContent = CifFile.readlines()
        CifFile.close()
    except:
        msg( "ERROR: cannot open " + CifName )
        return 2
    
    try:
        # get lattice parameters from cifcontent. the lattice parameters contains 6
        # float numbers which stand for a, b, c, alpha, beta, gamma
        LattConst = GetLattConst( CifContent )
        
        # extract atom type and atom coordinates from cifcontent
        AtomList, AtomCoord = GetAtomCoord( CifContent )
        
        # determine number of type, number of atoms and number of each type
        TypeList, NumList = GetTypeInfo( AtomList )
        NAtom = len(AtomList)
        NType = len(TypeList)
        
        # calculate lattice vectors of the system
        LattVec = CalcLattVec(LattConst)
        
        # shift to center according to dimension
        ShiftedCoord = ShiftCoord(AtomCoord,Dimension,ShiftStd)
        
        # convert to cartesian coordinates
        CartCoord = Frac2Cart( ShiftedCoord, LattVec )
    except:
        msg( "ERROR: illegal content of " + CifName )
        return 3
        
    try:
        # open cfg file for output
        CfgFile = open( CfgName, "w" )
        
        # write lattice constants
        CfgFile.write("lattice constants:\n")
        CfgFile.write("%7s%9.4f%9.4f\n" % ("a",     LattConst[0], LattConst[0]*ang2bohr))
        CfgFile.write("%7s%9.4f%9.4f\n" % ("b",     LattConst[1], LattConst[1]*ang2bohr))
        CfgFile.write("%7s%9.4f%9.4f\n" % ("c",     LattConst[2], LattConst[2]*ang2bohr))
        CfgFile.write("%7s%9.4f%9.4f\n" % ("alpha", LattConst[3], LattConst[3]*deg2rad))
        CfgFile.write("%7s%9.4f%9.4f\n" % ("beta",  LattConst[4], LattConst[4]*deg2rad))
        CfgFile.write("%7s%9.4f%9.4f\n" % ("gamma", LattConst[5], LattConst[5]*deg2rad))
        CfgFile.write("\n")
        
        # write number of type, number of atoms and number of each type
        CfgFile.write("%-16s%4d\n" % ("number of types:",NType))
        CfgFile.write("%-16s%4d\n" % ("number of atoms:",NAtom))
        CfgFile.write("number of atoms for each type:\n")
        for i in range(0,NType):
            CfgFile.write("%4s%4d\n" % (TypeList[i],NumList[i]))
        CfgFile.write("\n")
        
        # write lattice vectors in ang
        CfgFile.write("begin_vectors_ang\n")
        for v in LattVec:
            CfgFile.write("%14.9f%14.9f%14.9f\n" % (v[0],v[1],v[2]))
        CfgFile.write("end_vectors_ang\n")
        CfgFile.write("\n")
        
        # write lattice vectors in bohr
        CfgFile.write("begin_vectors_bohr\n")
        for v in LattVec:
            CfgFile.write("%16.9f%16.9f%16.9f\n" % (v[0]*ang2bohr,v[1]*ang2bohr,v[2]*ang2bohr))
        CfgFile.write("end_vectors_bohr\n")
        CfgFile.write("\n")
        
        # write frac coordinates
        CfgFile.write("begin_coordinates_frac\n")
        for i in range(0,NAtom):
            CfgFile.write("%4s%10.5f%10.5f%10.5f\n" % (AtomList[i],ShiftedCoord[i][0],ShiftedCoord[i][1],ShiftedCoord[i][2]))
        CfgFile.write("end_coordinates_frac\n")
        CfgFile.write("\n")
        
        # write cartesian coordinates in ang
        CfgFile.write("begin_coordinates_ang\n")
        for i in range(0,NAtom):
            CfgFile.write("%4s%10.5f%10.5f%10.5f\n" % (AtomList[i],CartCoord[i][0],CartCoord[i][1],CartCoord[i][2]))
        CfgFile.write("end_coordinates_ang\n")
        CfgFile.write("\n")
        
        # write cartesian coordinates in bohr
        CfgFile.write("begin_coordinates_bohr\n")
        for i in range(0,NAtom):
            CfgFile.write("%4s%12.5f%12.5f%12.5f\n" % (AtomList[i],CartCoord[i][0]*ang2bohr,CartCoord[i][1]*ang2bohr,CartCoord[i][2]*ang2bohr))
        CfgFile.write("end_coordinates_bohr\n")
        CfgFile.write("\n")
        
        # write command-line parameters
        CfgFile.write(CifName)
        CfgFile.write("\n")
        CfgFile.write(CfgName)
        CfgFile.write("\n")
        CfgFile.write("%-4d\n" % Dimension)
        CfgFile.write("%-9.4f\n" % ShiftStd)
        
        # close and quit
        CfgFile.close()
    except:
        msg( "ERROR: cannot open " + CfgName )
        return 4

    return 0

if __name__ == "__main__":
   sys.exit(main())
