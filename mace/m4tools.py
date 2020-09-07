import os
from os import system


def write_defs(macro, outfile="defs.m4"):
    with open(outfile, "w") as m4_file:
        m4_file.write("changequote([,])dnl\n")
        for key in sorted(macro.keys()):
            m4_file.write("define([%s], [%s])dnl\n" % (key, macro[key]))


def run_m4(infile, outfile, args=""):
    os.system("m4 %s %s | awk 'NF>0' > %s" % (args, infile, outfile))
