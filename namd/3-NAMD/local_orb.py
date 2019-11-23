#! /usr/bin/env python

from local import *

def read_max(maxname):
    with open(maxname, "r") as f:
        max_raw = f.readlines()
    return [int(line.split()[1]) for line in max_raw]

def get_bmin(shprop_name):
    with open(shprop_name) as shprop_file:
        shprop_raw = shprop_file.readlines()
    for line in shprop_raw:
        if line.find("BMIN") != -1:
            bmin = int(line.split()[3])
    return bmin

def calc_localorb(maxdat, bmin, postfix):
    shprop = read_shprop("SHPROP." + postfix)
    nstep = len(shprop)
    
    nl0 = int(postfix) - 1
    nl1 = nl0 + nstep
    max_post = maxdat[nl0:nl1]
    
    localorb = [ [i, shprop[i][max_post[i]-bmin+2]] for i in range(nstep) ]
    return localorb

def main():
    # read max.dat
    maxname = "max.dat"
    maxdat = read_max(maxname)
    
    # get the list of SHPROP.n files
    postfixlist = get_shprop_postfix()
    
    # average data
    bmin = get_bmin("SHPROP." + postfixlist[0])
    mat0 = calc_localorb(maxdat, bmin, postfixlist[0])
    for postfix in postfixlist[1:]:
        mat1 = calc_localorb(maxdat, bmin, postfix)
        sum_mat(mat0, mat1)
    average_mat(mat0, len(postfixlist))
    
    # output
    with open("local_orb.dat", "w") as f:
        for datarow in mat0:
            for data in datarow:
                f.write("%20.10e" % data)
            f.write("\n")

if __name__ == "__main__":
    main()