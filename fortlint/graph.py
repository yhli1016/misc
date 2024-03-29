"""Classes for dependency and call graphs."""
import re
import glob
import pickle
from collections import OrderedDict, defaultdict
from typing import Set, Callable

from .rules import Rules, DependRules


__all__ = ["DependGraph", "CallGraph"]


class Source:
    """
    Class representing a FORTRAN source file.

    Attributes
    ----------
    definitions: Set[str]
        symbols defined in this file
    references: Set[str]
        symbols referenced in this file
    dependencies: Set[str]
        file names where the referenced symbols are defined
    """
    def __init__(self) -> None:
        self.definitions = set()
        self.references = set()
        self.dependencies = set()


class DependGraph:
    """
    Class representing a tree of sources files.

    Attributes
    ----------
    _sources: OrderedDict[str, Source]
        collection of instances of the 'Source' class,
        keys: file names, values: instances
    _rules: ExtRules
        regex rules for detecting symbols
    """
    def __init__(self, rules: Rules = None) -> None:
        """
        :param rules: regex rules for detecting symbols
        """
        self._sources = OrderedDict()
        if rules is None:
            rules = DependRules()
        self._rules = rules

    def parse_source(self, src_name: str) -> Source:
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
            symbol = self._rules.match_def_start(line)
            if symbol is not None:
                source.definitions.add(symbol)
        source.references = self._rules.match_ref(content)
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
        # Get the file list under the directory
        # Omit the directory name for pwd. Keep it in other cases.
        if dir_name in (".", "./"):
            all_files = sorted(glob.glob("*"))
        else:
            all_files = sorted(glob.glob(f"{dir_name}/*"))

        # Parse FORTRAN source files
        pattern = re.compile(r"^(\S+)\.([fF]+\d*)$", re.IGNORECASE)
        for file in all_files:
            result = re.search(pattern, file)
            if result is not None:
                src_name = result.group(1)
                self._sources[src_name] = self.parse_source(file)

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
        for src_name, src_obj in self._sources.items():
            for item in src_obj.definitions:
                try:
                    def_table[item].add(src_name)
                except KeyError:
                    def_table[item] = {src_name}

        # Build dependencies
        for src_name, src_obj in self._sources.items():
            for item in src_obj.references:
                try:
                    candidates = def_table[item]
                    if src_name not in candidates:  # Skip self-dependence.
                        if len(candidates) == 1:
                            src_obj.dependencies.add(list(candidates)[0])
                        else:  # In the face of ambiguity, refuse to guess.
                            print(f"WARNING: skipping multiple definition for"
                                  f" {item} in {src_name}:\n\t{candidates}")
                except KeyError:
                    pass

        # Check for cyclic dependencies
        for src_name, src_obj in self._sources.items():
            dep_tree = list(src_obj.dependencies)
            for name2 in dep_tree:
                for dep in self._sources[name2].dependencies:
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
            pickle.dump(self._sources, f, pickle.HIGHEST_PROTOCOL)

    def load_cache(self, file_name: str = "sources.pkl") -> None:
        """
        Load saved sources from file.

        :param file_name: name of the file in which the sources are saved.
        :return: None. The 'sources' attribute is updated.
        """
        with open(file_name, 'rb') as f:
            self._sources = pickle.load(f)

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
            for src_name, src_obj in self._sources.items():
                if color_func is not None:
                    color = color_func(src_name)
                    dot.write(f"\"{src_name}\""
                              f" [color={color},fontcolor={color}]\n")
                for dep in sorted(src_obj.dependencies):
                    dot.write(f"\"{dep}\" -> \"{src_name}\"\n")
            dot.write("}\n")

    def write_make(self, make_name: str = "make.depend") -> None:
        """
        Write the dependencies between given source files to makefile.

        :param make_name: name of the makefile
        :return: None.
        """
        with open(make_name, "w") as make:
            for src_name, src_obj in self._sources.items():
                for dep in sorted(src_obj.dependencies):
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
            for src_name, src_obj in self._sources.items():
                if symbol in src_obj.definitions:
                    candidates.add(src_name)
        else:
            for src_name, src_obj in self._sources.items():
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
            for src_name, src_obj in self._sources.items():
                if re.search(pattern, src_name) is not None:
                    candidates = src_obj.dependencies
        else:
            for src_name, src_obj in self._sources.items():
                for dep in src_obj.dependencies:
                    if re.search(pattern, dep) is not None:
                        candidates.add(src_name)
                        break
        return candidates


class CallGraph:
    """
    Class for generating call graph.

    Attributes
    ----------
    _call_table: Dict[str, Set[str]
        symbols and references called by each symbol
    _rules: Rules
        regex rules for detecting symbols
    """
    def __init__(self, rules: Rules = None):
        """
        :param rules: regex rules for detecting symbols
        """
        self._call_table = defaultdict(set)
        if rules is None:
            rules = Rules()
        self._rules = rules

    def parse_source(self, src_name: str) -> None:
        """
        Parse a FORTRAN source file.

        NOTE: all the definitions and references in the source file will be
        converted to lower case. This is because identifiers in FORTRAN are
        case-insensitive.

        :param src_name: name of the source file
        :return: None
        """
        with open(src_name, "r") as in_file:
            content = in_file.readlines()

        # Get start and end line indices of definitions
        def_ln, num_lines = dict(), len(content)
        for i, line in enumerate(content):
            sym_def = self._rules.match_def_start(line)
            if sym_def is not None:
                for j in range(i+1, num_lines):
                    if self._rules.match_def_end(content[j]) is not None:
                        def_ln[sym_def] = (i, j)
                        break

        # Extract references within each definition
        for sym_def, ln in def_ln.items():
            refs = self._rules.match_ref(content[ln[0]+1:ln[1]])
            self._call_table[sym_def] = self._call_table[sym_def].union(refs)

    def parse_source_tree(self, dir_name: str = ".") -> None:
        """
        Parse all the source files under given directory.

        NOTE: the source name WILL NOT be converted to lower case, unlike
        the 'parse_source' method. This is because file names in UNIX are
        case-sensitive.

        :param dir_name: name of the directory
        :return: None. The 'sources' attribute is updated.
        """
        # Get the file list under the directory
        # Omit the directory name for pwd. Keep it in other cases.
        if dir_name in (".", "./"):
            all_files = sorted(glob.glob("*"))
        else:
            all_files = sorted(glob.glob(f"{dir_name}/*"))

        # Parse FORTRAN source files
        pattern = re.compile(r"^(\S+)\.([fF]+\d*)$", re.IGNORECASE)
        for file in all_files:
            result = re.search(pattern, file)
            if result is not None:
                self.parse_source(file)

    def write_dot(self, dot_name: str = "call.dot",
                  internal: bool = False) -> None:
        """
        Write the call graph to dot file for visualization.

        :param dot_name: name of the dot file
        :param internal: whether to output internal symbols only
        :return: None.
        """
        with open(dot_name, "w") as dot:
            symbols = self._call_table.keys()
            dot.write("digraph CallGraph {\nrankdir=LR\n")
            for symbol, ref in self._call_table.items():
                if internal:
                    ref = ref.intersection(symbols)
                for i in ref:
                    dot.write(f"\"{i}\" -> \"{symbol}\"\n")
            dot.write("}\n")
