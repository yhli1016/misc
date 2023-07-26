"""Classes regular expressions."""
import re
from typing import Union, List, Set


__all__ = ["Rules", "DependRules"]


class Rules:
    """
    Regex rules for extracting definitions and references.

    Attributes
    ----------
    _def_start: List[Tuple[str, re.Pattern, int]]
        rules for detecting starting line of definition
    _def_end: List[Tuple[str, re.Pattern, int]]
        rules for detecting ending line of definition
    _ref: List[Tuple[str, re.Pattern, int]]
        rules for detecting references
    _exclude: Set[str]
        reserved keywords to exclude from matched results
    """
    def __init__(self):
        # CAUTION: interface to operators must be placed prior to interface to
        # subroutines and functions.
        self._def_start = [
            ("sub", re.compile(r"^\s*(pure)?\s*subroutine\s+(\w+)", re.I), 2),
            ("func", re.compile(r"^\s*(pure)?\s*function\s+(\w+)", re.I), 2),
            ("int_op", re.compile(r"^\s*interface\s+operator\s*\((\s*\.\w+\.\s*)\)", re.I), 1),
            ("int_sub", re.compile(r"^\s*interface\s+(\w+)", re.I), 1),
        ]
        self._def_end = [
            ("sub", re.compile(r"^\s*end\s+subroutine", re.I), 1),
            ("func", re.compile(r"^\s*end\s+function", re.I), 1),
            ("int_op", re.compile(r"^\s*end\s+interface", re.I), 1),
            ("int_sub", re.compile(r"^\s*end\s+interface", re.I), 1),
        ]
        self._ref = [
            ("sub", re.compile(r"^\s*call\s+(\w+)", re.I), 1),
            ("func", re.compile(r"^\s*\w+\s*=\s*(\w+)\(.*\)", re.I), 1),
            ("op", re.compile(r"^\s*[^!].+(\.\w+\.)", re.I), 1),
            ("int", re.compile(r"^\s*module\s+procedure\s+([\w\s,]+)", re.I), 1),
        ]
        self._exclude = {".true.", ".false.", ".and.", ".or.", ".g.", ".gt.",
                         ".ge.", ".eq.", "operator", "abs", "sin", "cos", "dsin",
                         "dcos", "acos", "sqrt", "dcmplx", "cmplx", "size",
                         "aimag", "exp"}

    def match_def_start(self, line: str) -> Union[str, None]:
        """
        Check if the line is a starting line of symbol definition.

        :param line: line to check
        :return: matched symbol, otherwise none
        """
        sym_def = None
        for rule in self._def_start:
            kind, pattern, idx = rule
            result = re.search(pattern, line)
            if result is not None:
                s = result.group(idx).lstrip().rstrip().lower()
                if s not in self._exclude:
                    sym_def = s
                    break
        return sym_def

    def match_def_end(self, line: str) -> Union[str, None]:
        """
        Check if the line is an ending line of symbol definition.

        :param line: line to check
        :return: matched symbol, otherwise none
        """
        sym_def = None
        for rule in self._def_end:
            kind, pattern, idx = rule
            result = re.search(pattern, line)
            if result is not None:
                sym_def = "matched"
                break
        return sym_def

    def match_ref(self, content: List[str]) -> Set[str]:
        """
        Extract references from a list of lines.

        :param content: lines to check
        :return: extracted references
        """
        sym_refs = []
        for line in content:
            for rule in self._ref:
                # A single line may contain multiple references.
                # DO NOT add any break for this loop.
                kind, pattern, idx = rule
                result = re.search(pattern, line)
                if result is not None:
                    if kind != "int":
                        ref = result.group(idx).lstrip().rstrip().lower()
                        sym_refs.append(ref)
                    else:
                        ref = result.group(idx).split(",")
                        for s in ref:
                            sym_refs.append(s.lstrip().rstrip().lower())
        sym_refs = set(sym_refs).difference(self._exclude)
        return sym_refs


class DependRules(Rules):
    """Extended rules with module support."""
    def __init__(self):
        super().__init__()
        self._def_start.insert(0, ("mod", re.compile(r"^\s*module\s+(\w+)", re.I), 1))
        self._def_end.insert(0, ("mod", re.compile(r"^\s*end\s+module\s+(\w+)", re.I), 1))
        self._ref.insert(0, ("mod", re.compile(r"^\s*use\s+(\w+)", re.I), 1))
