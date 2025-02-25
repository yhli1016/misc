import re
import os
from collections import defaultdict
from typing import Dict, Tuple, List, Any


def read_lines(filename: str, nl0: int = None, nl1: int = None) -> List[str]:
    """
    Read and select file content according to line numbers.

    :param filename: name of the file
    :param nl0: starting line number, counted from 1
    :param nl1: ending line number, counted from 1
    :return: list of selected lines
    """
    with open(filename, "r") as infile:
        content = infile.readlines()
    nl_start = nl0 if nl0 is not None else 1
    nl_end = nl1 if nl1 is not None else len(content)
    return content[(nl_start-1):nl_end]


def include(filename: str, nl0: int = None, nl1: int = None) -> str:
    """
    Read and select file content according to line numbers.

    :param filename: name of the file
    :param nl0: starting line number, counted from 1
    :param nl1: ending line number, counted from 1
    :return: single line merged from selected lines
    """
    line = "".join(read_lines(filename, nl0, nl1))
    return line


def parse_block(filename: str) -> Dict[str, List[str]]:
    """
    Read and split file content into blocks according to tags.

    Beginning tag: /* begin xxx */
    Ending tag: /* end xxx */

    :param filename: name of the file
    :return: dict with keys being block names and values begin block contents
    """
    # Read raw content
    with open(filename, "r") as infile:
        content = infile.readlines()

    # Determine starting and ending line numbers
    begin_pattern = re.compile(r"^\s*/\s*\*\s*begin\s+(\w+)\s*\*\s*/$", re.I)
    end_pattern = re.compile(r"^\s*/\s*\*\s*end\s+(\w+)\s*\*\s*/$", re.I)
    begin_nl = dict()
    end_nl = dict()
    for i, line in enumerate(content):
        result = re.search(begin_pattern, line)
        if result is not None:
            begin_nl[result.group(1)] = i
        result = re.search(end_pattern, line)
        if result is not None:
            end_nl[result.group(1)] = i

    # Check if any block has missing line numbers
    diff = set(begin_nl.keys()).difference(end_nl.keys())
    if len(diff) > 0:
        for key in diff:
            raise RuntimeError(f"Ending tag of {key} not found")
    diff = set(end_nl.keys()).difference(begin_nl.keys())
    if len(diff) > 0:
        for key in diff:
            raise RuntimeError(f"Beginning tag of {key} not found")

    # Split content into blocks
    blocks = dict()
    for key, nl0 in begin_nl.items():
        nl1 = end_nl[key]
        blocks[key] = content[nl0:nl1+1]
    return blocks


class Mace:
    """
    Class for macro expansion.

    Attributes
    ----------
    _macro: Dict[str, Any]
        macro definitions
    _pattern_var: re.Pattern
        regular expression for detecting variable references
    _pattern_func: re.Pattern
        regular expression for detecting function calls
    """
    def __init__(self, macro: Dict[str, Any] = None) -> None:
        """
        :param macro: dictionary containing macro definitions
        """
        if macro is None:
            self._macro = dict()
        elif isinstance(macro, dict) or isinstance(macro, defaultdict):
            self._macro = macro
        else:
            raise TypeError("macro should be dict or defaultdict")
        self._pattern_var = re.compile(r"<\s*(\w+)\s*>")
        self._pattern_func = re.compile(r"<\s*(\w+\s*:.+)\s*>")

    def __getitem__(self, key: str) -> Any:
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
        Expand variable references in given line.

        :param line: line to expand
        :return: (line, status) expanded line and whether it differs from the
            original line
        """
        status = False
        result = re.findall(self._pattern_var, line)
        for item in result:
            try:
                var_text = self._macro[item]
            except KeyError:
                pass
            else:
                var_text = str(var_text).lstrip("\n").rstrip("\n")
                line = re.sub(f"<{item}>", var_text, line)
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
            for key, value in self._macro.items():
                m4_file.write(f"define(`{key}', `{value}')dnl\n")
        os.system(f"m4 {args} {template} | awk 'NF>0' > {output}")

    def line_replace(self, template: str, output: str) -> None:
        """
        Perform simple line replace on template file and save to output.

        :param template: name of template file
        :param output: name of output file
        :return: None
        """
        with open(template, "r") as in_file:
            raw_content = in_file.readlines()
        patterns = {_: re.compile(rf"^\s*<\s*{_}\s*>\s*$")
                    for _ in self._macro.keys()}
        with open(output, "w") as out_file:
            for line in raw_content:
                expanded = False
                for key, value in self._macro.items():
                    if re.search(patterns[key], line) is not None:
                        out_file.writelines(value)
                        expanded = True
                if not expanded:
                    out_file.write(line)

    def clear(self) -> None:
        """
        Clear all macro definitions.

        :return: None
        """
        if isinstance(self._macro, dict):
            self._macro = dict()
        else:
            self._macro = defaultdict(list)


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
