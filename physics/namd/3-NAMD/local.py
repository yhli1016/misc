#! /usr/bin/env python

def read_proj(projname):
    with open(projname, "r") as projfile:
        proj_raw = projfile.readlines()
    proj = []
    for line in proj_raw:
        proj.append([float(f) for f in line.split()])
    return proj

def get_shprop_postfix():
    with open("INICON") as inicon_file:
        inicon = inicon_file.readlines()
    postfix = [line.split()[0] for line in inicon]
    return postfix

def read_shprop(shprop_name):
    with open(shprop_name) as shprop_file:
        shprop_raw = shprop_file.readlines()
    
    # determine tnumber of KS states involved in the propagation
    comment  = [line for line in shprop_raw if line.find("#") != -1]
    for line in comment:
        if line.find("BMIN") != -1:
            bmin = int(line.split()[3])
        if line.find("BMAX") != -1:
            bmax = int(line.split()[3])
    nband = bmax - bmin + 1

    """
    CAUTION: there may be possible memory issues in the current implementation.
    But I cannot figure out any better algorithm.
    """
    
    # join all data lines into a very very long line and split it
    shprop_str = "".join(shprop_raw[len(comment):]).split()
    shprop_list = [float(f) for f in shprop_str]
    
    # convert shprop_list to matrix
    ndata = len(shprop_list)
    ndata_per_step = nband + 2
    shprop = [shprop_list[i:i+ndata_per_step] for i in range(ndata) \
              if i % ndata_per_step == 0]
    return shprop

def calc_localization(proj, postfix):
    shprop = read_shprop("SHPROP." + postfix)
    nstep = len(shprop)
    
    # extract projection data for SHPROP.postfix
    nl0 = int(postfix) - 1
    nl1 = nl0 + nstep
    proj_post = proj[nl0:nl1]
        
    # loop over all steps of SHPROP.post
    localization = []
    for i in range(nstep):
        rowp = proj_post[i][1:]
        rows = shprop[i][2:]
        localization.append([i, dot_product(rowp, rows)])
    return localization

def dot_product(a, b):
    result = 0.0
    for i in range(len(a)):
        result += a[i] * b[i]
    return result

def sum_mat(mat1, mat2):
    nrow, ncol = len(mat1), len(mat1[0])
    for i in range(nrow):
        for j in range(ncol):
            mat1[i][j] += mat2[i][j]

def average_mat(mat, dividend):
    nrow, ncol = len(mat), len(mat[0])
    for i in range(nrow):
        for j in range(ncol):
            mat[i][j] /= dividend

def main():
    # read proj.dat
    projname = "proj.dat"
    proj = read_proj(projname)
    
    # get the list of SHPROP.n files
    postfixlist = get_shprop_postfix()
    
    # average data
    mat0 = calc_localization(proj, postfixlist[0])
    for postfix in postfixlist[1:]:
        mat1 = calc_localization(proj, postfix)
        sum_mat(mat0, mat1)
    average_mat(mat0, len(postfixlist))
    
    # output
    with open("local.dat", "w") as f:
        for datarow in mat0:
            for data in datarow:
                f.write("%20.10e" % data)
            f.write("\n")

if __name__ == "__main__":
    main()