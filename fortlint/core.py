"""Core module of the package."""

import re
import glob
import pickle
from typing import Tuple, Union, Set, Callable


class SymbolRules:
    """
    Class for extracting definitions and references of symbols from lines
    of FORTRAN source code.

    Attributes
    ----------
    rules: List[Tuple[re.Pattern, int, bool]], rules for analyzing the
        source code
    """
    def __init__(self) -> None:
        self.rules = []

    def add_rule(self, pattern: str, idx: int = 1,
                 is_def: bool = False) -> None:
        """
        Add a rule to the rules.

        :param pattern: regular expression of the rule
        :param idx: index to extract the symbol from the matched result
        :param is_def: whether the matched result contains a definition
        :return: None. The 'rules' attribute is updated.
        """
        pattern = re.compile(pattern, re.IGNORECASE)
        self.rules.append((pattern, idx, is_def))

    def match(self, line: str) -> Tuple[Union[str, None], bool]:
        """
        Loop over the rules to match given line of FORTRAN source code.

        :param line: line of FORTRAN source code
        :return: the extracted symbol (None if not matched) and if the
            matched result contains a definition
        """
        symbol, is_def = None, False
        for rule in self.rules:
            result = re.search(rule[0], line)
            if result is not None:
                symbol = result.group(rule[1]).lstrip().rstrip().lower()
                is_def = rule[2]
                break
        return symbol, is_def


# Predefined matching rules
RULES = SymbolRules()
# module definition
RULES.add_rule(r"^\s*module\s+(\w+)", is_def=True)
# subroutine definition
RULES.add_rule(r"^\s*(pure)?\s*subroutine\s+(\w+)", idx=2, is_def=True)
# function definition
RULES.add_rule(r"^\s*(pure)?\s*function\s+(\w+)", idx=2, is_def=True)
# interface to subroutines and functions
RULES.add_rule(r"^\s*interface\s+(\w+)", is_def=True)
# interface to operators
RULES.add_rule(r"^\s*interface operator\s*\((\s*\.\w+\.\s*)\)", is_def=True)
# module usage
RULES.add_rule(r"^\s*use\s+(\w+)")
# subroutine call
RULES.add_rule(r"^\s*call\s+(\w+)")
# function call
RULES.add_rule(r"^\s*\w+\s*=\s*(\w+)\(.*\)")
# operator call
RULES.add_rule(r"^\s*[^!].+(\.\w+\.)")


class Source:
    """
    Class representing a FORTRAN source file.

    Attributes
    ----------
    definitions: set of strings, containing the definitions in this file
    references: set of strings, containing the references to the definitions
    dependencies: set of strings, containing the file names where the
        definitions referenced within this file can be found
    """
    def __init__(self) -> None:
        self.definitions = set()
        self.references = set()
        self.dependencies = set()

    def add_definition(self, definition: str) -> None:
        """
        Add a definition to the definitions.

        :param definition: the definition
        :return: None.
        """
        self.definitions.add(definition)

    def add_reference(self, reference: str) -> None:
        """
        Add a reference to the references.

        :param reference: the reference
        :return: None.
        """
        self.references.add(reference)

    def add_dependency(self, dependency: str) -> None:
        """
        Add a dependent source file to the dependencies.

        :param dependency: name of the source file
        :return: None.
        """
        self.dependencies.add(dependency)


class SourceTree:
    """
    Class representing a tree of sources files.

    Attributes
    ----------
    sources: Dict[str, Source], collection instances of 'Source' class
    """
    def __init__(self):
        self.sources = dict()

    @staticmethod
    def parse_source(src_name: str) -> Source:
        """
        Parse a FORTRAN source file.

        NOTE: all the definitions and references in the source file will be
        converted to lower case. See also the 'find_symbol' method. This is
        because identifiers in FORTRAN are case-insensitive.

        :param src_name: name of the source file
        :return: the source object created from the file
        """
        source = Source()
        with open(src_name, "r") as in_file:
            content = in_file.readlines()
        for line in content:
            symbol, is_def = RULES.match(line)
            if symbol is not None:
                if is_def:
                    source.add_definition(symbol)
                else:
                    source.add_reference(symbol)
        return source

    def parse_source_tree(self, dir_name: str = ".") -> None:
        """
        Parse all the source files under given directory.

        NOTE: the source name WILL NOT be converted to lower case, unlike
        the 'parse_source' method. This is because file names in UNIX are
        case-sensitive.

        :param dir_name: name of the directory
        :return: None. The 'sources' attribute is updated.
        """
        pattern = re.compile(r"^(\S+)\.([fF]+\d*)$", re.IGNORECASE)
        # Omit the directory name for pwd. Keep it otherwise.
        if dir_name in (".", "./"):
            all_files = sorted(glob.glob("*"))
        else:
            all_files = sorted(glob.glob(f"{dir_name}/*"))
        for file in all_files:
            result = re.search(pattern, file)
            if result is not None:
                src_name = result.group(1)
                self.sources[src_name] = self.parse_source(file)

    def resolve_dependencies(self) -> None:
        """
        Resolve the dependencies between given source files.

        :return: None. The 'Source' instances in 'sources' attribute
            are updated.
        """
        # Build the definition tables
        # key: the definitions
        # value: list of sources where the definitions can be found
        def_table = dict()
        for src_name, src_obj in self.sources.items():
            for item in src_obj.definitions:
                try:
                    def_table[item].add(src_name)
                except KeyError:
                    def_table[item] = {src_name}

        # Build dependencies
        for src_name, src_obj in self.sources.items():
            for item in src_obj.references:
                try:
                    candidates = def_table[item]
                    if src_name not in candidates:  # Skip self-dependence.
                        if len(candidates) == 1:
                            src_obj.add_dependency(list(candidates)[0])
                        else:  # In the face of ambiguity, refuse the temptation to guess.
                            print(f"WARNING: skipping multiple definition for"
                                  f" {item} in {src_name}:\n\t{candidates}")
                except KeyError:
                    pass

        # Check for cyclic dependencies
        for src_name, src_obj in self.sources.items():
            dep_tree = list(src_obj.dependencies)
            for name2 in dep_tree:
                for dep in self.sources[name2].dependencies:
                    if dep not in dep_tree:
                        dep_tree.append(dep)
            if src_name in dep_tree:
                print(f"WARNING: cyclic dependency detected for {src_name}")

    def save_cache(self, file_name: str = "sources.pkl") -> None:
        """
        Saves the sources to file.

        :param file_name: name of the file in which the sources are saved.
        :return: None.
        """
        with open(file_name, "wb") as f:
            pickle.dump(self.sources, f, pickle.HIGHEST_PROTOCOL)

    def load_cache(self, file_name: str = "sources.pkl") -> None:
        """
        Load saved sources from file.

        :param file_name: name of the file in which the sources are saved.
        :return: None. The 'sources' attribute is updated.
        """
        with open(file_name, 'rb') as f:
            self.sources = pickle.load(f)

    def write_dot(self, dot_name: str = "dep.dot",
                  color_func: Callable[[str], str] = None) -> None:
        """
        Write the dependencies between given source files to dot file for
        visualization.

        :param dot_name: name of the dot file
        :param color_func: function receiving the src_name and returning
            the color
        :return: None.
        """
        with open(dot_name, "w") as dot:
            dot.write("digraph Dependency {\nrankdir=LR\n")
            for src_name, src_obj in self.sources.items():
                if color_func is not None:
                    color = color_func(src_name)
                    dot.write(f"\"{src_name}\""
                              f" [color={color},fontcolor={color}]\n")
                for dep in src_obj.dependencies:
                    dot.write(f"\"{dep}\" -> \"{src_name}\"\n")
            dot.write("}\n")

    def write_make(self, make_name: str = "make.depend") -> None:
        """
        Write the dependencies between given source files to makefile.

        :param make_name: name of the makefile
        :return: None.
        """
        with open(make_name, "w") as make:
            for src_name, src_obj in self.sources.items():
                for dep in src_obj.dependencies:
                    make.write(f"{src_name}.o: {dep}.o\n")

    def find_symbol(self, symbol: str, kind: str = "def") -> Set[str]:
        """
        Find the definitions or references of given symbol in sources.

        :param symbol: the symbol to search
        :param kind: definition or reference to search
        :return: set of the file names which contain definitions or reference
            of the symbol
        """
        candidates = set()
        symbol = symbol.lower()
        if kind == "def":
            for src_name, src_obj in self.sources.items():
                if symbol in src_obj.definitions:
                    candidates.add(src_name)
        else:
            for src_name, src_obj in self.sources.items():
                if symbol in src_obj.references:
                    candidates.add(src_name)
        return candidates

    def find_relevant_nodes(self, node: str, direction="in") -> Set[str]:
        """
        Find the relevant nodes of given node in the dependency digraph.

        :param node: name of the given node
        :param direction: 'in' for nodes on which the given node depends
            and 'out' for nodes which depend on the given node
        :return: set of names of relevant nodes
        """
        candidates = set()
        pattern = re.compile(f"^{node}$", re.IGNORECASE)
        if direction == "in":
            for src_name, src_obj in self.sources.items():
                if re.search(pattern, src_name) is not None:
                    candidates = src_obj.dependencies
        else:
            for src_name, src_obj in self.sources.items():
                for dep in src_obj.dependencies:
                    if re.search(pattern, dep) is not None:
                        candidates.add(src_name)
                        break
        return candidates
