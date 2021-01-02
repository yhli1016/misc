class TagSet(object):
    """
    Class for holding "key = value" pairs, for generating input files or printing calculation
    parameters to stdout.

    self.style is the pre-defined style of input file, currently supports qe, cp2k, vasp, siesta,
    and other programs that follow a similar input file syntax.

    The input file consists of several nodes, each being in the following format:
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

    Format for cp2k
    ! Input for SCF calculation    <- comment
    &CONTROL                       <- head
    prefix = 'si'                  <- key1, value1
    calculation = 'scf'            <- key2, value2
    /END CONTROL                   <- tail

    Format for vasp
    # Input for SCF calculation    <- comment
    prefix = 'si'                  <- key1, value1
    calculation = 'scf'            <- key2, value2
                                   <- tail

    Siesta is similar to vasp, but with '=' omitted.

    Tag set can have its own subsets (children), which can have different styles than their parent. Particularly,
    QE and CP2K may have subsets in VASP/SIESTA style. So we provided two pre-defined styles, namely qe_plain and
    p2k_plain, for QE and CP2K respectively. These styles are similar to VASP/SIESTA, but with comments starting
    with '!' instead of '#'.
    """
    def __init__(self, style="qe", comment=None, name=None):
        """
        :param style: string, pre-defined style for programs
        :param comment: string, comment
        :param name: string, name of the tag set, for generating self.head and self.tail
        """
        # Set default value for all attributes.
        self.style = style
        self.comment, self.head, self.tail = None, None, None
        self.assign_marker = "="
        self.tags = []
        self.subsets = []

        # Update attributes according to style.
        if style in ("qe", "cp2k"):
            # QE and CP2K uses name list, which has strict rules.
            if comment is not None:
                self.comment = "! %s" % comment
            if name is not None:
                self.head = "&%s" % name
                self.tail = "/" if style == "qe" else "/END %s\n" % name
            else:
                raise RuntimeError("Name is required for '%s' style tag sets." % style)
            if style == "cp2k":
                self.assign_marker = ""

        # QE_PLAIN and CP2K_PLAIN are aimed to present VASP style subsets in input of
        # QE and CP2K, so this piece of code is similar to that of VASP/SIESTA. The only
        # difference is that we use '!' instead of '#'.
        elif style in ("qe_plain", "cp2k_plain"):
            if comment is not None:
                self.comment = "! %s" % comment
            if name is not None:
                self.head, self.tail = "! %s" % name, ""
            if style == "cp2k_plain":
                self.assign_marker = ""

        # VASP and SIESTA have rather free input format. The different between the codes is
        # that SIESTA does not require an '='.
        elif style in ("vasp", "siesta"):
            if comment is not None:
                self.comment = "# %s" % comment
            if name is not None:
                self.head, self.tail = "# %s" % name, ""
            if style == "siesta":
                self.assign_marker = ""

        # Style for reporting calculation details, in which we use no '#' or '!'.
        elif style in ("report",):
            if comment is not None:
                self.comment = comment
            if name is not None:
                self.head, self.tail = name, ""
            self.assign_marker = ":"

        else:
            raise NotImplementedError("Style '%s' not implemented yet." % style)

    def add_tag(self, tag_name, tag_val=""):
        """
        Add an tag to self.tags.

        :param tag_name: string, tag name
        :param tag_val: any data_kind that can be converted to string, tag value
        """
        self.tags.append([tag_name, tag_val])

    def add_subset(self, subset):
        """
        Add existing tag set to self.subsets.

        :param subset: instance of TagSet class, sunset
        :return: None
        """
        self.subsets.append(subset)

    def get_max_tag_length(self, include_subsets=True):
        """
        Get the maximum length of tags.

        :param include_subsets: boolean, whether to take the tags of subsets in consideration
        :return: integer, maximum length of tags
        """
        tag_length = [len(tag[0]) for tag in self.tags]
        if include_subsets:
            for subset in self.subsets:
                tag_length.append(subset.get_max_tag_length(include_subsets))
        return max(tag_length)

    @staticmethod
    def write_indent(out_file, indent):
        """
        Write indent to file.

        :param out_file: file object or stdout, destination of writing
        :param indent: integer, number of spaces
        :return:
        """
        for i in range(indent):
            out_file.write(" ")

    def write_tags(self, out_file, indent_level=0, align_tags=True, align_to_subsets=True,
                   max_tag_length="auto"):
        """
        Write tags to file or stdout.

        :param out_file: file object or stdout, destination of writing
        :param indent_level: integer, indent level
        :param align_tags: boolean, whether to align tags for better appearance
        :param align_to_subsets: boolean, whether to consider subsets in the alignment
        :param max_tag_length: integer, maximum tag length
        :return: None
        """
        # Calculate indents
        # Comment and head use indent_head, while (key, tag) pairs uses indent_tag.
        space_indent = 4
        # QE and QE_PLAIN use fixed indent
        if self.style == "qe":
            indent_head, indent_tag = 0, space_indent
        elif self.style == "qe_plain":
            indent_head = indent_tag = space_indent

        # CP2K and CP2K_PLAIN use incremental indent
        elif self.style == "cp2k":
            indent_head, indent_tag = indent_level * space_indent, (indent_level + 1) * space_indent
        elif self.style == "cp2k_plain":
            indent_head = indent_tag = indent_level * space_indent

        # VASP/SIESTA use no indent
        elif self.style in ("vasp", "siesta"):
            indent_head = indent_tag = 0

        # Report use similar style as QE
        elif self.style == "report":
            indent_head, indent_tag = 0, space_indent

        else:
            raise NotImplementedError("Style '%s' not implemented yet." % self.style)

        # Determine tag format.
        if align_tags:
            # If no max_tag_length is given, compute the max tag length of self.
            # Otherwise, keep it as given in kwargs.
            if max_tag_length == "auto":
                max_tag_length = self.get_max_tag_length(align_to_subsets)
            format_tag = "%%-%ds %s %%s\n" % (max_tag_length, self.assign_marker)
        else:
            format_tag = "%%s %s %%s\n" % self.assign_marker

        # Write comment, head and tags to out file.
        if self.comment is not None:
            self.write_indent(out_file, indent_head)
            out_file.write("%s\n" % self.comment)
        if self.head is not None:
            self.write_indent(out_file, indent_head)
            out_file.write("%s\n" % self.head)
        for tag in self.tags:
            self.write_indent(out_file, indent_tag)
            out_file.write(format_tag % (tag[0], tag[1]))

        # Write subsets to file.
        for subset in self.subsets:
            indent_level_subset = indent_level + 1
            if not align_to_subsets:
                max_tag_length = "auto"
            #out_file.write("\n")
            subset.write_tags(out_file, indent_level=indent_level_subset, align_tags=align_tags,
                              align_to_subsets=align_to_subsets, max_tag_length=max_tag_length)

        # Write tail
        if self.tail is not None:
            self.write_indent(out_file, indent_head)
            out_file.write("%s\n" % self.tail)
