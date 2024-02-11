from collections import OrderedDict
from typing import List, Any


class TagSet:
    """
    Class for holding "key = value" pairs, for generating input files or
    printing calculation parameters to stdout.

    The input file consists of several nodes, each being in the following
    format:

    [comment]
    [head]
    key1 [=] value1
    key2 [=] value2
    ...      ...
    [tail]

    Format for QE:
    ! Input for SCF calculation    <- comment
    &CONTROL                       <- head
    prefix = 'si'                  <- key1, value1
    calculation = 'scf'            <- key2, value2
    /                              <- tail

    Format for CP2K
    ! Input for SCF calculation    <- comment
    &CONTROL                       <- head
    prefix = 'si'                  <- key1, value1
    calculation = 'scf'            <- key2, value2
    /END CONTROL                   <- tail

    Format for VASP
    # Input for SCF calculation    <- comment
    prefix = 'si'                  <- key1, value1
    calculation = 'scf'            <- key2, value2
                                   <- tail

    SIESTA is similar to VASP, but with '=' omitted.

    Attributes
    ----------
    _style: str
        pre-defined tag style
    _comment: str
        comment line for the tag-set
    _head: str
        headline for the tag-set
    _tail: str
        tail line for the tag-set
    _assign_marker: str
        marker for assignment
    _tags: OrderedDict[str, Any]
        names and values of tags
    _subsets: List[TagSet]
        child tag-sets
    """
    def __init__(self, style: str = "qe",
                 comment: str = None,
                 name: str = None) -> None:
        """
        :param style: pre-defined style for programs
        :param comment: comment line
        :param name: name of the tag set, for generating head and tail
        """
        # Set default value for all attributes.
        self._style = style
        self._comment, self._head, self._tail = None, None, None
        self._assign_marker = "="
        self._tags = OrderedDict()
        self._subsets = []

        # Update attributes according to style.
        if style in ("qe", "cp2k"):
            # QE and CP2K use FORTRAN name lists
            if comment is not None:
                self._comment = f"! {comment}"
            if name is not None:
                self._head = f"&{name}"
                self._tail = "/" if style == "qe" else f"/END {name}\n"
            else:
                raise RuntimeError(f"Name required for {style} style tag sets")
            if style == "cp2k":
                self._assign_marker = ""

        # QE_SUB and CP2K_SUB for representing subsets in FORTRAN name lists
        elif style in ("qe_sub", "cp2k_sub"):
            if comment is not None:
                self._comment = f"! {comment}"
            if name is not None:
                self._head = f"! {name}"
                self._tail = ""
            if style == "cp2k_sub":
                self._assign_marker = ""

        # VASP and SIESTA have rather free input format.
        elif style in ("vasp", "siesta"):
            if comment is not None:
                self._comment = f"# {comment}"
            if name is not None:
                self._head = f"# {name}"
                self._tail = ""
            if style == "siesta":
                self._assign_marker = ""
        else:
            raise NotImplementedError(f"Style {style} not implemented" % style)

    def __setitem__(self, key: str, value: Any) -> None:
        self._tags[key] = value

    def __getitem__(self, key: str) -> Any:
        return self._tags[key]

    def add_subset(self, subset: Any) -> None:
        self._subsets.append(subset)

    def get_max_tag_length(self, with_subsets: bool = True) -> int:
        """
        Get the maximum length of tags.

        :param with_subsets: whether to consider the tags of subsets
        :return: maximum length of tags
        """
        tag_length = [len(_) for _ in self._tags.keys()]
        if with_subsets:
            for subset in self._subsets:
                tag_length.append(subset.get_max_tag_length(with_subsets))
        return max(tag_length)

    def write(self, out_file,
              indent_level: int = 0,
              align_tags: bool = True,
              align_subsets: bool = True,
              tag_length: int = -1) -> None:
        """
        Write tags to file.

        :param out_file: output file of stdout
        :param indent_level: indent level
        :param align_tags: whether to align tags for better appearance
        :param align_subsets: whether to consider subsets in the alignment
        :param tag_length: length of tags, 0 for automatic detection
        :return: None
        """
        # Determine indents
        # Comment and head use indent_head, while (key, tag) pairs uses
        # indent_tag.
        num_space = 4

        # QE and QE_SUB use fixed indent
        if self._style == "qe":
            indent_head = 0
            indent_tag = num_space
        elif self._style == "qe_sub":
            indent_head = num_space
            indent_tag = num_space
        # CP2K and CP2K_SUB use incremental indent
        elif self._style == "cp2k":
            indent_head = indent_level * num_space
            indent_tag = indent_head + num_space
        elif self._style == "cp2k_sub":
            indent_head = indent_level * num_space
            indent_tag = indent_head
        # VASP/SIESTA use no indent
        else:
            indent_head = 0
            indent_tag = 0

        # Write comment, head and tags to out file.
        if self._comment is not None:
            out_file.write(" " * indent_head)
            out_file.write(f"{self._comment}\n")
        if self._head is not None:
            out_file.write(" " * indent_head)
            out_file.write(f"{self._head}\n")
        if align_tags and tag_length == -1:
            tag_length = self.get_max_tag_length(align_subsets)
        for key, value in self._tags.items():
            out_file.write(" " * indent_tag)
            if align_tags:
                line = f"{key:<{tag_length}s} {self._assign_marker} {value}\n"
            else:
                line = f"{key} {self._assign_marker} {value}\n"
            out_file.write(line)

        # Write subsets
        for subset in self._subsets:
            subset.write(out_file,
                         indent_level=indent_level+1,
                         align_tags=align_tags,
                         align_subsets=align_subsets,
                         tag_length=tag_length)

        # Write tail
        if self._tail is not None:
            out_file.write(" " * indent_head)
            out_file.write(f"{self._tail}\n")


def main():
    import sys
    ts = TagSet(style="cp2k", name="CONTROL")
    ts["a"] = 10
    ts["b"] = 11
    ts2 = TagSet(style="cp2k", name="subset")
    ts2["xx"] = 1
    ts2["abc"] = 3
    ts.add_subset(ts2)
    ts.write(sys.stdout, align_subsets=True)


if __name__ == "__main__":
    main()
