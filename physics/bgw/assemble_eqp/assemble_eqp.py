#! /usr/bin/env python

# The sigma code of BerkeleyGW runs parallelly on bands, but serially on kpoints.
# In some cases it causes huge RAM usage, and we have to break a large task into
# several smaller subtasks with respect to bands. The aim of this program is to
# assemble the eqp*.dat files produced by such sub-tasks and write to one file.
#
# Usage: assemble_eqp.py output input1 input2 input3 ......
#
# NOTE: THE ORDER OF INPUT FILES MUST BE CONSISTENT WITH BAND INDICES!!!

import sys
import os

def ExtractKpoint(eqp):
    kptlist = []
    for line in eqp:
        s = line.split()
        x = float(s[0])
        y = float(s[1])
        z = float(s[2])
        if abs(x) < 1 and abs(y) < 1 and abs(z) < 1:
            kptlist.append([x, y, z])
    return kptlist

def main():
    # display notice
    print("")
    print("NOTE: The order of eqp*.dat must be consistent with the band indices.")
    print("Also, the k-points must be consistent among files. Please check them")
    print("carefully.")
    print("")
    
    # parse cli-parameters
    outfnm = sys.argv[1]
    eqpfnm = sys.argv[2:]

    # check if outfile exists
    if os.path.exists(outfnm):
        print("ERROR: " + outfnm + " already exists!")
        print("Program aborted to avoid accidental overwritting.")
        sys.exit(-1)
    
    # read eqp*.dat
    eqp_raw = []
    for fnm in eqpfnm:
        try:
            infile = open(fnm, "r")
            content = infile.readlines()
            infile.close()
        except:
            print("ERROR: cannot open " + fnm)
            sys.exit(-1)
        else:
            eqp_raw.append(content)
    
    # extract k-points from eqp_raw[0]
    kptlist = ExtractKpoint(eqp_raw[0])
    nk = len(kptlist)

    # open file for output
    ofile = open(outfnm, "w")
    
    # loop over kpoints
    for ik in range(nk):
        kpt = kptlist[ik]
        
        # loop over eqp to assemble eqp_k
        eqp_k = []
        for eqp in eqp_raw:
            nbnd = int(eqp[0].split()[3])
            nl0 = ik * (nbnd + 1) + 1
            nl1 = nl0 + nbnd
            for line in eqp[nl0:nl1]:
                eqp_k.append(line)
        
        # write to file
        ofile.write("%13.9f%13.9f%13.9f%8d\n" % (kpt[0], kpt[1], kpt[2], len(eqp_k)))
        for line in eqp_k:
            ofile.write(line)
    
    # close outfile
    ofile.close()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
