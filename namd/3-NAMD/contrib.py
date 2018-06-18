#! /usr/bin/env python

from local import *

def calc_contrib(proj, postfix):
    shprop = read_shprop("SHPROP." + postfix)
    nstep = len(shprop)
    
    # extract projection data for SHPROP.postfix
    nl0 = int(postfix) - 1
    nl1 = nl0 + nstep
    proj_post = proj[nl0:nl1]
    
    # calculate difference of proj_post shprop
    dproj = diff(proj_post)
    dsh = diff(shprop)
    
    """
    NOTE: the last row of proj_post and shprop is dropped when evaluating the
    difference, so dproj and dsh are one row shorter than their counter-parts.
    That's why i is in range(nstep-1), not range(nstep) in the code blocks below.
    """
    # the adiabatic contribution is sum(dproj * shprop)
    contrib_ad = []
    for i in range(nstep-1):
        rowp = dproj[i][1:]
        rows = shprop[i][2:]
        contrib_ad.append([i, dot_product(rowp, rows)])
    
    # the non-adiabatic contribution is sum(proj_post * dsh)
    contrib_nonad = []
    for i in range(nstep-1):
        rowp = proj_post[i][1:]
        rows = dsh[i][2:]
        contrib_nonad.append([i, dot_product(rowp, rows)])
    
    # integrate contrib_ad and contrib_nonad
    contrib_ad_int = integrate(contrib_ad)
    contrib_nonad_int = integrate(contrib_nonad)
    
    return contrib_ad_int, contrib_nonad_int
    
def diff(mat):
    """
    Note: For all columns df = f(i+1) - f(i). The last row is dropped.
    
    The 1st column (time) for proj and the 1-2 columns (time and energy)
    for shprop won't be a problem, as they are not involved in the summation.
    """
    nrow, ncol = len(mat), len(mat[0])
    df = [ [mat[i+1][j] - mat[i][j] for j in range(ncol)] \
            for i in range(nrow-1) ]
    return df

def integrate(mat):
    """
    NOTE: It is assumed that mat has two columns. The first column is kept unchaged
    during the integration.
    """
    nrow = len(mat)
    x = [mat[i][0] for i in range(nrow)]
    y = [mat[i][1] for i in range(nrow)]
    integ = [[x[i], sum(y[:i+1])] for i in range(nrow)]
    return integ
    
def main():
    # read proj.dat
    projname = "proj.dat"
    proj = read_proj(projname)
    
    # get the list of SHPROP.n files
    postfixlist = get_shprop_postfix()
    
    # average data
    contrib_ad, contrib_nonad = calc_contrib(proj, postfixlist[0])
    for postfix in postfixlist[1:]:
        temp_ad, temp_nonad = calc_contrib(proj, postfix)
        sum_mat(contrib_ad, temp_ad)
        sum_mat(contrib_nonad, temp_nonad)
    average_mat(contrib_ad, len(postfixlist))
    average_mat(contrib_nonad, len(postfixlist))
    
    # output
    with open("contrib_ad.dat", "w") as f:
        for datarow in contrib_ad:
            for data in datarow:
                f.write("%20.10e" % data)
            f.write("\n")
    #
    with open("contrib_nonad.dat", "w") as f:
        for datarow in contrib_nonad:
            for data in datarow:
                f.write("%20.10e" % data)
            f.write("\n")

if __name__ == "__main__":
    main()
