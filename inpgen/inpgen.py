#! /usr/bin/env python
import os


def write_defs(macro, outfile="defs.m4"):
    with open(outfile, "w") as m4_file:
        m4_file.write("changequote([,])dnl\n")
        m4_file.write("changecom(\"/*\", \"*/\")dnl\n")
        for key in macro.keys():
            m4_file.write("define([%s], [%s])dnl\n" % (key, macro[key]))


def run_m4(infile, outfile, args=""):
    os.system("m4 %s %s > %s" % (args, infile, outfile))


def get_input():
    macro = dict()

    # Get the type of job
    while True:
        job = input("\nInput type of job (opt/neb/bader): ")
        if job in ("opt", "neb", "bader"):
            macro[".job"] = job
            break

    # Get specific job information
    if job == "neb":
        num_image = int(input("\nInput number of transition states: "))
        macro["NIMAGE"] = num_image
        macro["DIR_TS"] = "".join(["%02d " % _ for _ in range(1, num_image + 1)])
        macro["DIR_TOT"] =  "".join(["%02d " % _ for _ in range(num_image + 2)])
        macro["NMAX"] = "%02d" % (num_image + 1)

    # Get general job information
    macro["NAME"] = input("\nInput job name: ")
    macro["NCPU"] = input("\nInput number of cpu to use: ")
    macro["TIME"] = input("\nInput time limit (in hours): ")
    while True:
        restart = input("\nInput whether to restart (yes/no): ")
        if restart in ("yes", "no"):
            break
    if restart == "yes":
        macro["RESTART"] = True
    macro["RUN"] = input("\nInput number of run: ")
    while True:
        write_incar = input("\nInput whether to write incar (yes/no): ")
        if write_incar in ("yes", "no"):
            macro[".write_incar"] = write_incar
            break
    return macro


def main():
    # Location of template files
    inpgen_dir = "%s/soft/inpgen" % os.environ["HOME"]
    script_dir = "%s/euler" % inpgen_dir
    incar_dir = "%s/incar" % inpgen_dir

    # Get arguments from stdin
    macro = get_input()

    # Run m4 to generate input file
    write_defs(macro)
    run_m4("%s/run_%s.m4" % (script_dir, macro[".job"]), "run_%s.sh" % macro[".job"])
    if macro[".write_incar"] == "yes":
        run_m4("%s/INCAR_%s.m4" % (incar_dir, macro[".job"]), "INCAR")
    os.system("rm defs.m4")

    # Prompt the user
    print("\nScript written to '%s'" % "run_%s.sh" % macro[".job"])
    if macro[".write_incar"] == "yes":
        print("\nInput written to 'INCAR'")


if __name__ == "__main__":
    main()
