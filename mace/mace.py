#! /usr/bin/env python
"""MACro Expansion program."""

import sys
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pythonpath", type=str, action="store",
                        help="search PYTHONPATH for macro definition")
    parser.add_argument("template", type=str, action="store",
                        help="template for generating file")
    parser.add_argument("output", type=str, action="store",
                        help="file generated from macro and template")
    args = parser.parse_args()
    return args


def expand_content(macro, content_raw):
    content_expanded = []
    for line in content_raw:
        flag = True
        while flag is True:
            line, flag = expand_line(macro, line)
        content_expanded.append(line)
    return content_expanded


def expand_line(macro, line):
    flag = False
    for mac_name, mac_text in macro.items():
        word = "<%s>" % mac_name
        if line.find(word) is not -1:
            line = line.replace(word, str(mac_text).lstrip("\n").rstrip("\n"))
            flag = True
    return line, flag


def main(macro, template, output):
    try:
        with open(template, "r") as in_file:
            content_raw = in_file.readlines()
    except IOError:
        print("ERROR: cannot read '%s'" % template)
        sys.exit(-1)
    content_expanded = expand_content(macro, content_raw)
    try:
        with open(output, "w") as out_file:
            for line in content_expanded:
                out_file.write(line)
    except IOError:
        print("ERROR: cannot write to '%s'" % output)
        sys.exit(-1)


if __name__ == "__main__":
    args = parse_args()
    try:
        sys.path.append(args.pythonpath)
        import defs
    except ImportError:
        print("ERROR: Cannot import defs.py")
        sys.exit(-1)
    else:
        main(defs.macro, args.template, args.output)
