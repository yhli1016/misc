import re
import os
from typing import Dict, Tuple


def include(filename: str, nl0: int = None, nl1: int = None) -> str:
    """
    Read given file return the contents as a single line.

    :param filename: name of the file
    :param nl0: starting line number, counting from 1
    :param nl1: ending line number, counting from 1
    :return:
    """
    with open(filename, "r") as infile:
        content = infile.readlines()
    nl_start = nl0 if nl0 is not None else 1
    nl_end = nl1 if nl1 is not None else len(content)
    line = "".join(content[(nl_start-1):nl_end])
    return line


def _expand_line(macro: Dict[str, str], line: str) -> Tuple[str, bool]:
    """
    Expand macros in given line.

    :param macro: definition of macros, with keys being names and values being
        contents
    :param line: line to be expanded
    :return: (line, status) where line is the expanded line and flag is whether
        expansion has been performed
    """
    result = re.search(r"<[a-zA-Z0-9_]+>", line)
    status = False
    if result is not None:
        name = result.group()[1:-1]
        try:
            text = str(macro[name]).lstrip("\n").rstrip("\n")
        except KeyError:
            pass
        else:
            line = re.sub(f"<{name}>", text, line)
            status = True
    return line, status


def expand_file(macro: Dict[str, str], template: str, output: str) -> None:
    """
    Expand macros in template and write to output.

    :param macro: definition of macros, with keys being names and values being
        contents
    :param template: filename of template
    :param output: filename of output
    :return: None
    """
    with open(template, "r") as in_file:
        content_raw = in_file.readlines()
    content_expanded = []
    for line in content_raw:
        status = True
        while status:
            line, status = _expand_line(macro, line)
        content_expanded.append(line)
    with open(output, "w") as out_file:
        out_file.writelines(content_expanded)


def write_defs(macro: Dict[str, str], defs: str = "defs.m4") -> None:
    """
    Write macro definitions to a m4 file.

    :param macro: definition of macros, with keys being names and values being
        contents
    :param defs: filename of the definition file
    :return: None
    """
    with open(defs, "w") as m4_file:
        m4_file.write("changequote([,])dnl\n")
        for key, value in macro.items():
            m4_file.write(f"define([{key}], [{value}])dnl\n")


def run_m4(template: str, output: str, args: str = "") -> None:
    """
    Run m4 to expand macros in template and write to output.

    :param template: filename of template
    :param output: filename of output
    :param args: additional arguments for m4
    :return:
    """
    os.system(f"m4 {args} {template} | awk 'NF>0' > {output}")


def main() -> None:
    """
    Main function when running mace as a script.

    :return: None
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,  action="store", required=True)
    parser.add_argument("-o", "--output", type=str, action="store", required=True)
    parser.add_argument("-m", "--macro", type=str, action="store", nargs="*")
    args = parser.parse_args()

    macro = dict()
    for argv in args.macro:
        s = argv.split("=")
        if len(s) == 1:
            macro[s[0]] = ""
        else:
            macro[s[0]] = s[1]
    expand_file(macro, args.input, args.output)


if __name__ == "__main__":
    main()
