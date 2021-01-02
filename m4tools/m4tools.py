import os


def include(filename, nl0=None, nl1=None):
    with open(filename, "r") as infile:
        content = infile.readlines()
    nl_start = nl0 if nl0 is not None else 1
    nl_end   = nl1 if nl1 is not None else len(content)
    longline = "".join(content[(nl_start-1):nl_end])
    return longline


def write_defs(macro, outfile="defs.m4"):
    with open(outfile, "w") as m4_file:
        m4_file.write("changequote([,])dnl\n")
        for key in sorted(macro.keys()):
            m4_file.write("define([%s], [%s])dnl\n" % (key, macro[key]))


def dict2args(macro):
    args = "".join(["-D%s=%s " % (key, macro[key]) for key in macro.keys()])
    return args


def run_m4(infile, outfile, args=""):
    os.system("m4 %s %s | awk 'NF>0' > %s" % (args, infile, outfile))
