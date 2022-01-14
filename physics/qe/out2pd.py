#! /usr/bin/env python

# Last modified on 2016-08-18
#
# This program extracts namelists and relaxed lattice constants / atomic 
# coordinates from relax and vc-relax runs.
#
# Usage: out2pd.py in out pd

import sys

def ReadInp(inpname):
    infile = open(inpname, "r")
    infile_content = infile.readlines()
    infile.close()
    inp = []
    for line in infile_content:
        if line.find("!") == -1 and line.find("#") == -1:
            inp.append(line)
    return inp

def ReadOut(outname):
    infile = open(outname, "r")
    infile_content = infile.readlines()
    infile.close()
    nl0 = infile_content.index("Begin final coordinates\n")
    nl1 = infile_content.index("End final coordinates\n")
    return infile_content[nl0:nl1]

def ExtractNamelist(inp, flaglist):
    result = []
    for flag in flaglist:
        for line in inp:
            if line.find(flag) != -1:
                result.append(line)
    return result

def ExtractUPF(inp):
    result = []
    for line in inp:
        if line.find("UPF") != -1:
            result.append(line)
    return result

def ExtractCellInp(inp):
    for line in inp:
        if line.find("CELL_PARAMETERS") != -1:
            nl = inp.index(line)
    return inp[nl+1:nl+4]

def ExtractCellOut(out):
    for line in out:
        if line.find("CELL_PARAMETERS") != -1:
            nl = out.index(line)
    return out[nl+1:nl+4]

def ExtractCoord(out):
    for line in out:
        if line.find("ATOMIC_POSITIONS") != -1:
            nl = out.index(line)
    return out[nl+1:]

def Norm(x):
    while x < 0.0:
        x += 1.0
    while x > 1.0:
        x -= 1.0
    return x

def main():
    # parse command-line parameters
    inpname = sys.argv[1]
    outname = sys.argv[2]
    pdname  = sys.argv[3]
    
    # read inp and out
    inp = ReadInp(inpname)
    out = ReadOut(outname)
    
    # extract namelists
    foo = ["prefix", "tefield", "dipfield"]
    nml_control = ExtractNamelist(inp, foo)
    foo = ["ibrav", "celldm", "ntyp", "nat", "edir",
          "emaxpos", "eopreg", "eamp", "input_dft"]
    nml_system = ExtractNamelist(inp, foo)
    
    # extract UPF information
    upf_info = ExtractUPF(inp)
    
    # check if this is a vc-relax run
    isvc = False
    for line in inp:
        if line.find("vc-relax") != -1:
            isvc = True
    
    # extract CELL_PARAMETERS according to isvc
    if isvc == True:
        cell_param = ExtractCellOut(out)
    else:
        cell_param = ExtractCellInp(inp)
    
    # extract coordinates
    coord = ExtractCoord(out)
    
    # output
    ofile = open(pdname, "w")
    
    # CONTROL
    ofile.write("&CONTROL\n")
    for line in nml_control:
        ofile.write(line)
    ofile.write("\\\n")
    
    # SYSTEM
    ofile.write("&SYSTEM\n")
    for line in nml_system:
        ofile.write(line)
    ofile.write("\\\n")
    
    # UPF
    for line in inp:
        if line.find("ATOMIC_SPECIES") != -1:
            ofile.write(line)
    for line in upf_info:
        ofile.write(line)
        
    # CELL_PARAMETERS
    for line in inp:
        if line.find("CELL_PARAMETERS") != -1:
            ofile.write(line)
    for line in cell_param:
        ofile.write(line)
    
    # ATOMIC_POSITIONS
    for line in inp:
        if line.find("ATOMIC_POSITIONS") != -1:
            ofile.write(line)
    for line in coord:
        s = line.split()
        ofile.write("%4s%14.9f%14.9f%14.9f\n" % (s[0], Norm(float(s[1])), Norm(float(s[2])), Norm(float(s[3]))))
    #
    ofile.close()
    return 0

if __name__ == "__main__":
   sys.exit(main())
