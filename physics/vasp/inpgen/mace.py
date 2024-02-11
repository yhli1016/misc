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


class Mace:
    """
    Class for macro expansion.

    Attributes
    ----------
    _macro: Dict[str, str]
        macro definitions
    _pattern: re.Pattern
        pre-compiled regular expressions
    """
    def __init__(self, macro: Dict[str, str]) -> None:
        """
        :param macro: dictionary containing macro definitions
        """
        self._macro = macro
        self._pattern = re.compile(r"<\w+>")

    def _expand_line(self, line: str) -> Tuple[str, bool]:
        """
        Expand macros in given line using regular expression.

        :param line: line to be expanded
        :return: (line, status) where line is the expanded line and flag is
            whether expansion has been performed
        """
        result = re.search(self._pattern, line)
        status = False
        if result is not None:
            name = result.group()[1:-1]
            try:
                text = str(self._macro[name]).lstrip("\n").rstrip("\n")
            except KeyError:
                pass
            else:
                line = re.sub(f"<{name}>", text, line)
                status = True
        return line, status

    def expand_re(self, template: str, output: str) -> None:
        """
        Expand macros in template and write to output using regular expression.

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
                line, status = self._expand_line(line)
            content_expanded.append(line)
        with open(output, "w") as out_file:
            out_file.writelines(content_expanded)

    def expand_m4(self, template: str, output: str, args: str = "") -> None:
        """
        Expand macros in template and write to output using m4.

        :param template: filename of template
        :param output: filename of output
        :param args: additional arguments for m4
        :return: None
        """
        with open("defs.m4", "w") as m4_file:
            m4_file.write("changequote([,])dnl\n")
            for key, value in self._macro.items():
                m4_file.write(f"define([{key}], [{value}])dnl\n")
        os.system(f"m4 {args} {template} | awk 'NF>0' > {output}")


def main() -> None:
    """
    Main function when running mace as a script.

    :return: None
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,  action="store",
                        required=True)
    parser.add_argument("-o", "--output", type=str, action="store",
                        required=True)
    parser.add_argument("-m", "--macro", type=str, action="store",
                        nargs="*")
    args = parser.parse_args()

    macro = dict()
    for argv in args.macro:
        s = argv.split("=")
        if len(s) == 1:
            macro[s[0]] = ""
        else:
            macro[s[0]] = s[1]
    Mace(macro).expand_re(args.input, args.output)


if __name__ == "__main__":
    main()
