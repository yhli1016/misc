#! /usr/bin/env python

# Last modified on 2015-10-13
#
# This program searches for optimal numbers of bands to be used in epsilon
# and sigma calculations for BerkeleyGW.
#
# Usage: optimal.py -c ncpu -v nval -a nbandmin -b nbandmax
# where ncpu is the number of cpu, nval is the number of valence bands,
# nbandmin and nbandmax are lower and upper bounds for search.

import sys
import getopt

def usage():
    print
    print "optimal.py usage:"
    print
    print "-c, --ncpu: number of cpu"
    print "-v, --nval: number of valence bands"
    print "-a, --nmin: lower bound of nband"
    print "-b, --nmax: upper bound of nband"
    print


def parse(argv):
    ncpu = 1
    nval = 1
    nmin = 1
    nmax = 1
    try:
        opts, args = getopt.getopt(argv[1:], "c:v:a:b:", ["ncpu=", "nval=", "nmin=", "nmax="])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(1)
    for opt, val in opts:
        if opt in ("-c", "--ncpu"):
            ncpu = int(val)
        elif opt in ("-v", "--nval"):
            nval = int(val)
        elif opt in ("-a", "--nmin"):
            nmin = int(val)
        elif opt in ("-b", "--nmax"):
            nmax = int(val)
        else:
            print "unhandled option", opt
            usage()
            sys.exit(2)
    return ncpu, nval, nmin, nmax


def main(argv):
    
    ncpu, nval, nmin, nmax = parse(argv)
    
    # determine optimal nbnd for epsilon
    print "optimal number of bands for epsilon:"
    for i in range(nmin, nmax+1):
        if (i-nval)%ncpu == 0:
            print "%8d%8d%8d" % (i, nval, i-nval)
    
    # determine optimal nbnd for sigma
    print "optimal number of bands for sigma:"
    for i in range(nmin, nmax+1):
        if i%ncpu == 0:
            print "%8d%8d%8d" %(i, nval, i-nval)


if __name__ == '__main__':
    main(sys.argv)
