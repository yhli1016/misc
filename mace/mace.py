"""
Library for generating files using MACro Expansion.

Usage: main(macro, template, ouput)
    macro: dictionary, macro definitions
    template: string, filename of the template
    output: string, filename of the output file

An auxiliary function include() is provided for including a file and defining
a macro from it.
"""


def expand_line(macro, line):
    flag = False
    for mac_name, mac_text in macro.items():
        word = "<%s>" % mac_name
        if line.find(word) is not -1:
            line = line.replace(word, str(mac_text).lstrip("\n").rstrip("\n"))
            flag = True
    return line, flag


def include(filename, nl0=None, nl1=None):
    with open(filename, "r") as infile:
        content = infile.readlines()
    nl_start = nl0 if nl0 is not None else 1
    nl_end   = nl1 if nl1 is not None else len(content)
    longline = "".join(content[(nl_start-1):nl_end])
    return longline


def main(macro, template, output):
    with open(template, "r") as in_file:
        content_raw = in_file.readlines()

    content_expanded = []
    for line in content_raw:
        flag = True
        while flag is True:
            line, flag = expand_line(macro, line)
        content_expanded.append(line)

    with open(output, "w") as out_file:
        for line in content_expanded:
            out_file.write(line)
