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
    os.system("rm defs.m4")


def main():
    print("*****************************************************************")
    print("\n    NEVER FORGET TO UPDATE THE TEMPLATES AFTER INSTALLATION!\n")
    print("*****************************************************************")
    # Get arguments from stdin
    macro = dict()

    # Get the type of job
    while True:
        job = input("Please input type of job (opt/neb): ")
        if job in ("opt", "neb"):
            break

    # Get specific job information
    if job == "neb":
        num_image = int(input("Please input the number of transition states: "))
        macro["DIR_TS"] = "".join(["%02d " % _ for _ in range(1, num_image + 1)])
        macro["DIR_TOT"] =  "".join(["%02d " % _ for _ in range(num_image + 2)])
        macro["NMAX"] = "%02d" % (num_image + 1)

    # Get general job information
    macro["NAME"] = input("Please input name of the job: ")
    macro["NCPU"] = input("Please input the number of cpu to use: ")
    macro["TIME"] = input("Please input the time limit for the job (in hours): ")
    while True:
        restart = input("Please input whether to restart (yes/no): ")
        if restart in ("yes", "no"):
            break
    if restart == "yes":
        macro["RESTART"] = True
    macro["RUN"] = input("Please input the number of run: ")

    # Run m4 to generate input file
    write_defs(macro)
    m4_dir = "%s/proj/misc/inpgen/euler" % os.environ["HOME"]
    run_m4("%s/run_%s.m4" % (m4_dir, job), "run_%s.sh" % job, args="-I./ -I%s" % m4_dir)
    print("Script written to '%s'" % "run_%s.sh" % job)


if __name__ == "__main__":
    main()
