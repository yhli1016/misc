#! /usr/bin/env python


# This program generates a list of hopping terms up to a given distance
# by looping over all orbitals in [0,0,0] and in neighbouring unit cells.
#
# Onsite terms are defined as i == j and R == [0, 0, 0]. If 'reduce_onsite'
# is set to True in the inp, they will removed from the list.
#
# Equivalent conjugate terms are defined as <i|H|j+R> = <j|H|i-R>*.
# If 'reduce_conjugate' is set to True in the inp, they will be removed
# from the list.


import math
import operator
import inp


def writedat(filename, dat):
    ofile = open(filename, "w")
    for line in dat:
        ofile.write("%14.9f%14.9f%14.9f\n" % (line[0], line[1], line[2]))
    ofile.close()


def frac2cart(basevec, fraccoord):
    x0 = fraccoord[0]
    y0 = fraccoord[1]
    z0 = fraccoord[2]
    x1 = x0 * basevec[0][0] + y0 * basevec[1][0] + z0 * basevec[2][0]
    y1 = x0 * basevec[0][1] + y0 * basevec[1][1] + z0 * basevec[2][1]
    z1 = x0 * basevec[0][2] + y0 * basevec[1][2] + z0 * basevec[2][2]
    cartcoord = [x1, y1, z1]
    return cartcoord


def check_onsite(hopi):
    flag = True
    i1 = hopi[0]
    j1 = hopi[1]
    Rn = [hopi[2], hopi[3], hopi[4]]
    if (i1==j1) and (Rn==[0,0,0]):
        flag = False
    return flag


def check_conj(hop, hopi):
    flag = True
    i1 = hopi[0]
    j1 = hopi[1]
    #
    # check hopi against exsiting items in hop, <i|H|j+R> is equivalent to <j|H|i-R>
    # as <i|H|j+R> = <j|H|i-R>*.
    #
    # Here we use a tricky algorithm. As equivalent hopping terms have equal 
    # distances and all the hopping terms have been sorted in hop, it will avoid
    # a lot of unnecessary comparisons if we reverse the loop and start from the
    # tail of hop.
    #
    for hopk in reversed(hop):
        i2 = hopk[0]
        j2 = hopk[1]
        dR = [hopi[2]+hopk[2], hopi[3]+hopk[3], hopi[4]+hopk[4]]
        if (i1==j2) and (j1==i2) and (dR==[0,0,0]):
            flag = False
            break
    return flag


def calc_dr(orbi, orbj, i, j, k):
    rfrac = [inp.orb[orbj][0] - inp.orb[orbi][0] + i,\
             inp.orb[orbj][1] - inp.orb[orbi][1] + j,\
             inp.orb[orbj][2] - inp.orb[orbi][2] + k]
    rcart = frac2cart(inp.lat, rfrac)
    dr = math.sqrt(rcart[0]**2 + rcart[1]**2 + rcart[2]**2)
    return dr

def main():
    # create the unfiltered hopping list
    norb = len(inp.orb)
    rawhop = [ [orbi, orbj, i, j, k, calc_dr(orbi, orbj, i, j, k)] \
               for orbi in range(norb) \
               for orbj in range(norb) \
               for i in range(inp.namin, inp.namax+1) \
               for j in range(inp.nbmin, inp.nbmax+1) \
               for k in range(inp.ncmin, inp.ncmax+1) \
               if calc_dr(orbi, orbj, i, j, k) <= inp.dmax ]

    # sort rawhop
    sorthop = sorted(rawhop, key=operator.itemgetter(5))

    # filt sorthop
    if inp.reduce_onsite is True:
        filthop1 = []
        for hopi in sorthop:
            if check_onsite(hopi):
                filthop1.append(hopi)
    else:
        filthop1 = sorthop
    #
    if inp.reduce_conjugate is True:
        filthop2 = []
        for hopi in filthop1:
            # CAUTION: Here we check whether hopi is in filthop2, not in filthop1!
            if check_conj(filthop2, hopi):
                filthop2.append(hopi)
    else:
        filthop2 = filthop1

    # write filthop2 to hopterm.dat
    hopfile = open("hopterm.dat", "w")
    if inp.output_flavor == "matlab":
        idcount = 1
    elif inp.output_flavor == "python":
        idcount = 0
    else:
        print("ERROR: Unknown output_flavor!")
        return -1
    for hopi in filthop2:
        hopfile.write("%4d%4d%4d%4d%4d%14.9f\n" % (hopi[0]+idcount, hopi[1]+idcount, hopi[2], hopi[3], hopi[4], hopi[5]))
    hopfile.close()
    
    # if specified, write lat and orb to lat.dat and orb.dat
    if inp.write_lat is True:
        writedat("lat.dat", inp.lat)
    if inp.write_orb is True:
        writedat("orb.dat", inp.orb)

    return 0

if __name__ == "__main__":
    main()
