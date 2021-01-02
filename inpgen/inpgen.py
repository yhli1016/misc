#! /usr/bin/env python
import os


M4_DIR = "%s/inpgen/euler" % os.environ["HOME"]


def write_defs(macro, outfile="defs.m4"):
    with open(outfile, "w") as m4_file:
        m4_file.write("changequote([,])dnl\n")
        for key in sorted(macro.keys()):
            m4_file.write("define([%s], [%s])dnl\n" % (key, macro[key]))


def run_m4(infile, outfile, args=""):
    os.system("m4 %s %s > %s" % (args, infile, outfile))
    os.system("rm defs.m4")


def get_param():
    param = dict()

    # Get the type of job
    while True:
        param["job"] = input("Please input type of job (opt/neb): ")
        if param["job"] in ("opt", "neb"):
            break

    # Get general job information
    param["name"] = input("Please input name of the job: ")
    param["ncpu"] = input("Please input the number of cpu to use: ")
    param["time"] = input("Please input the time limit for the job (in hours): ")
    while True:
        param["restart"] = input("Please input whether to restart (yes/no): ")
        if param["restart"] in ("yes", "no"):
            break
    param["run"] = input("Please input the number of run: ")
    
    # Get specific job information
    if param["job"] == "neb":
        param["num_image"] = int(input("Please input the number of transition states: "))
    else:
        param["num_image"] = 0
    return param


def main():
    # Get arguments from stdin
    param = get_param()

    # Create the macro definition dictionary
    macro = dict()

    # General job information
    macro["NAME"] = param["name"]
    macro["NCPU"] = param["ncpu"]
    macro["TIME"] = param["time"]
    if param["restart"] == "yes":
        macro["RESTART"] = "YES"
    macro["RUN"] = param["run"]
    
    # Job specific information
    if param["job"] == "neb":
        macro["DIR_TS"] = "".join(["%02d " % _ for _ in range(1, param["num_image"] + 1)])
        macro["DIR_TOT"] =  "".join(["%02d " % _ for _ in range(param["num_image"] + 2)])
        macro["NMAX"] = "%02d" % (param["num_image"] + 1)

    # Run m4 to generate input file
    write_defs(macro)
    run_m4("%s/run_%s.m4" % (M4_DIR, param["job"]), "run_%s.sh" % param["job"])
    print("Script written to '%s'" % "run_%s.sh" % param["job"])


if __name__ == "__main__":
    main()
