#! /usr/bin/env python

from local import *

def main():
    # extact the postfixes of SHPROP.* files
    postfixlist = get_shprop_postfix()
    
    # average data
    mat0 = read_shprop("SHPROP." + postfixlist[0])
    for postfix in postfixlist[1:]:
        mat1 = read_shprop("SHPROP." + postfix)
        sum_mat(mat0, mat1)
    average_mat(mat0, len(postfixlist))
    
    # output
    with open("average.dat", "w") as f:
        for datarow in mat0:
            for data in datarow:
                f.write("%20.10e" % data)
            f.write("\n")

if __name__ == "__main__":
    main()
