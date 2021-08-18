#! /usr/bin/env python
import re
import os


def include(filename, nl0=None, nl1=None):
    """
    Read selected lines [nl0, nl1] of given file and merge into a long line.

    :param filename: string
        name of the file to read
    :param nl0: integer
        line number of the starting line
    :param nl1: integer
        line number of the ending line
    :return: one_line, string
        selected content of file in a long line
    """
    with open(filename, "r") as infile:
        content = infile.readlines()
    nl_start = nl0 if nl0 is not None else 1
    nl_end = nl1 if nl1 is not None else len(content)
    one_line = "".join(content[(nl_start-1):nl_end])
    return one_line


def expand_line(macro, line):
    """
    Expand macros in given line and reports whether modification has been made.

    :param macro: dictionary
        definition of macros, with keys being macro name and value being macro
        content
    :param line: string
        line to be searched for macros and expanded
    :return line: string
        modified line if it contains macros, otherwise the original line
    :return flag: boolean
        whether the line has been modified
    """
    flag = False
    match_result = re.search(r"<[a-zA-Z0-9_]+>", line)
    if match_result is not None:
        mac_name = match_result.group()[1:-1]
        if mac_name in macro.keys():
            mac_text = str(macro[mac_name]).lstrip("\n").rstrip("\n")
            line = re.sub(r"<%s>" % mac_name, mac_text, line)
            flag = True
    return line, flag


def expand_file(macro, template, output):
    """
    Expand macros in given template file and write to output, for generation
    of input files for ab initio codes.

    :param macro: dictionary
        definition of macros, with keys being macro name and value being macro
        content
    :param template: string
        file name of template
    :param output: string
        file name of output
    :return: None
    """
    with open(template, "r") as in_file:
        content_raw = in_file.readlines()
    content_expanded = []
    for line in content_raw:
        while True:
            line, flag = expand_line(macro, line)
            if not flag:
                break
        content_expanded.append(line)
    with open(output, "w") as out_file:
        for line in content_expanded:
            out_file.write(line)


def get_option(prompt, options):
    """
    Get one option from a list of options.

    :param prompt: string
        prompt to the user
    :param options: tuple or list of strings
        available options
    :return: option, string
        chosen option
    """
    while True:
        option = input(prompt)
        if option in options:
            break
    return option


def get_input():
    """
    Get input parameters as a macro.

    NOTE: keys in lowercase and starting with '.' are reserved for
    special purposes. DO NOT use them as macro names in templates.

    :return: macro, dictionary
        all the input parameters
    """
    macro = dict()

    # Get the type of job
    macro[".job"] = get_option("\nInput type of job (opt/neb/bader): ",
                               ("opt", "neb", "bader"))

    # Get specific job information
    if macro[".job"] == "neb":
        num_image = int(input("\nInput number of transition states: "))
        macro["NIMAGE"] = num_image
        macro["DIR_TS"] = "$(seq -f '%%02g' %d %d)" % (1, num_image)
        macro["DIR_TOT"] = "$(seq -f '%%02g' %d %d)" % (0, num_image + 1)
        macro["NMAX"] = "%02d" % (num_image + 1)

    # Get general job information
    macro["NAME"] = input("\nInput job name: ")
    macro["NCPU"] = input("\nInput number of cpu to use: ")
    macro["TIME"] = input("\nInput time limit (in hours): ")
    restart = get_option("\nInput whether to restart (yes/no): ",
                         ("yes", "no"))
    if restart == "yes":
        macro["RESTART"] = 1
    else:
        macro["RESTART"] = 0
    macro["RUN"] = input("\nInput number of run: ")
    macro[".incar"] = get_option("\nInput whether to write incar (yes/no): ",
                                 ("yes", "no"))
    return macro


def main():
    # Location of template files
    root_dir = "%s/soft/inpgen" % os.environ["HOME"]
    script_dir = "%s/euler" % root_dir
    incar_dir = "%s/incar" % root_dir

    # Get arguments
    macro = get_input()

    # Generate input file
    expand_file(macro, "%s/run_%s.m4" % (script_dir, macro[".job"]),
                "run_%s.sh" % macro[".job"])
    if macro[".incar"] == "yes":
        expand_file(macro, "%s/INCAR_%s.m4" % (incar_dir, macro[".job"]),
                    "INCAR")

    # Notify the user
    print("\nScript written to 'run_%s.sh'" % macro[".job"])
    if macro[".incar"] == "yes":
        print("\nInput written to 'INCAR'")


if __name__ == "__main__":
    main()
