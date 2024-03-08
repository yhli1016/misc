import re
import os
from typing import Dict, Tuple, List, Any


def include(filename: str, nl0: int = None, nl1: int = None) -> str:
    """
    Read file and return the content in a single line.

    :param filename: name of the file
    :param nl0: starting line number, counted from 1
    :param nl1: ending line number, counted from 1
    :return: file content
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
    _macro: Dict[str, Any]
        macro definitions
    _pattern: re.Pattern
        pre-compiled regular expressions
    """
    def __init__(self, macro: Dict[str, Any] = None) -> None:
        """
        :param macro: dictionary containing macro definitions
        """
        if macro is not None:
            self._macro = macro
        else:
            self._macro = dict()
        self._pattern = re.compile(r"<\w+>")

    def __getitem__(self, key: str) -> str:
        """
        Get macro content.

        :param key: name of the macro
        :return: content of the macro
        """
        return self._macro[key]

    def __setitem__(self, key: str, value: Any) -> None:
        """
        Add or update a macro.

        :param key: name of the macro
        :param value: content of the macro
        :return: None
        """
        self._macro[key] = value

    def _expand(self, line: str) -> Tuple[str, bool]:
        """
        Expand macros in given line once.

        :param line: line to be expanded
        :return: (line, status) where line is the expanded line and flag is
            whether expansion has been performed
        """
        result = re.findall(self._pattern, line)
        status = False
        for item in result:
            name = item[1:-1]
            try:
                text = str(self._macro[name]).lstrip("\n").rstrip("\n")
            except KeyError:
                pass
            else:
                line = re.sub(f"<{name}>", text, line)
                status = True
        return line, status

    def expand_lines(self, lines: List[str]) -> List[str]:
        """
        Expand macros in incoming lines.

        :param lines: incoming lines to expand
        :return: expanded lines
        """
        lines_expanded = []
        for line in lines:
            status = True
            while status:
                line, status = self._expand(line)
            lines_expanded.append(line)
        return lines_expanded

    def expand_file(self, template: str, output: str) -> None:
        """
        Expand macros in template file and save to output.

        :param template: name of template file
        :param output: name of output file
        :return: None
        """
        with open(template, "r") as in_file:
            content_raw = in_file.readlines()
        content_expanded = self.expand_lines(content_raw)
        with open(output, "w") as out_file:
            out_file.writelines(content_expanded)

    def expand_file_m4(self, template: str, output: str, args: str = "") -> None:
        """
        Expand macros in template file and save to output.

        :param template: name of template file
        :param output: name of output file
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

    macro = Mace()
    for argv in args.macro:
        s = argv.split("=")
        if len(s) == 1:
            macro[s[0]] = ""
        else:
            macro[s[0]] = s[1]
    macro.expand_file(args.input, args.output)


if __name__ == "__main__":
    main()
