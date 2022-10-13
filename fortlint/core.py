"""
Core module of the package

CONSTANTS
---------
    PATTERNS: dictionary for parsing the source files
    KEYS_DEF: dictionary for extracting definitions from source file

Classes
-------
    Source: developer class for representing a source file
    SourceTree: user class for representing a tree of sources files
"""

import re
import glob
import pickle
from typing import Set, Callable


# Patterns for detecting the symbols defined and referenced in the source file.
# The precise pattern for 'func' should be r"^\s*(pure\s+)?\s*function\s+(\w+)"
# but that makes the 'parse_source' method of 'SourceTree' class complicated
# when extracting the symbol. So we have to make a compromise.
PATTERNS = {
    'mod': re.compile(r"^\s*module\s+(\w+)", re.IGNORECASE),
    'sub': re.compile(r"^\s*subroutine\s+(\w+)", re.IGNORECASE),
    'func': re.compile(r"^\s*[pure]?\s*function\s+(\w+)", re.IGNORECASE),
    'interface': re.compile(r"^\s*interface\s+(\w+)", re.IGNORECASE),
    'interface_op': re.compile(r"^\s*interface operator\s*\((\s*\.\w+\.\s*)\)",
                               re.IGNORECASE),
    'use': re.compile(r"^\s*use\s+(\w+)", re.IGNORECASE),
    'call': re.compile(r"^\s*call\s+(\w+)", re.IGNORECASE),
    'call_func': re.compile(r"^\s*\w+\s*=\s*(\w+)\(.*\)", re.IGNORECASE),
    'call_op': re.compile(r"^\s*[^!].+(\.\w+\.)", re.IGNORECASE)
}
KEYS_DEF = ('mod', 'sub', 'func', 'interface', 'interface_op')


class Source:
    """
    Class representing a FORTRAN source file.

    Attributes
    ----------
    symbols: set of strings, containing the modules and subroutines defined
        in the source file
    references: set of strings, containing the modules in 'use' and subroutines
        in 'call' statements
    dependencies: set of strings, containing the names of files where
        the dependent modules and subroutines are defined
    """
    def __init__(self) -> None:
        self.symbols = set()
        self.references = set()
        self.dependencies = set()

    def add_symbol(self, symbol: str):
        """
        Add a symbol to the set of symbols.

        :param symbol: the symbol
        :return: None.
        """
        self.symbols.add(symbol)

    def add_reference(self, reference: str):
        """
        Add a reference to the set of references.

        :param reference: the reference
        :return: None.
        """
        self.references.add(reference)

    def add_dependency(self, dependency: str):
        """
        Add a dependent source file to the set of dependencies.

        :param dependency: name of the source file
        :return:
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
    def parse_source(source_name: str) -> Source:
        """
        Parse a FORTRAN source file.

        NOTE: all the symbols and references in the source file will be
        converted to lower case. See also the 'find_symbol' method. This
        is because identifiers in FORTRAN are case-insensitive.

        :param source_name: name of the source file
        :return: the source object created from the file
        """
        source = Source()
        with open(source_name, "r") as in_file:
            content = in_file.readlines()
        for line in content:
            for key, val in PATTERNS.items():
                result = re.search(val, line)
                if result is not None:
                    pattern = result.group(1).lstrip().rstrip().lower()
                    if key in KEYS_DEF:
                        source.add_symbol(pattern)
                    else:
                        source.add_reference(pattern)
        return source

    def parse_source_tree(self, dir_name: str = "*") -> None:
        """
        Parse all the source files under given directory.

        NOTE: the source name WILL NOT be converted to lower case, unlike
        the 'parse_source' method. This is because file names in UNIX are
        case-sensitive.

        :param dir_name: name of the directory
        :return: None. The 'sources' attribute is updated.
        """
        pattern = re.compile(r"^(\S+)\.([fF]+\d*)$", re.IGNORECASE)
        all_files = glob.glob(f"{dir_name}/*")
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
        # Build the symbol tables
        # key: the symbols
        # value: list of sources where the symbols are defined
        symbol_table = dict()
        for src_name, src_obj in self.sources.items():
            for symbol in src_obj.symbols:
                try:
                    symbol_table[symbol].add(src_name)
                except KeyError:
                    symbol_table[symbol] = {src_name}

        # Build dependencies
        for src_name, src_obj in self.sources.items():
            for symbol in src_obj.references:
                try:
                    candidates = symbol_table[symbol]
                    if src_name not in candidates:  # Skip self-dependence.
                        if len(candidates) == 1:
                            src_obj.add_dependency(list(candidates)[0])
                        else:  # In the face of ambiguity, refuse the temptation to guess.
                            print(f"WARNING: skipping multiple definition for"
                                  f" {symbol} in {src_name}:\n\t{candidates}")
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
        :param kind: definition of reference to search
        :return: set of the file names which contain definitions or reference
            of the symbol
        """
        candidates = set()
        symbol = symbol.lower()
        if kind == "def":
            for src_name, src_obj in self.sources.items():
                if symbol in src_obj.symbols:
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
