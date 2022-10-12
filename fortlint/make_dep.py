#! /usr/bin/env python
"""Script for generating the dependency digraph."""

import re
import argparse

from fortlint import SourceTree


def color_func(src_name):
    """Custom coloring function."""
    if re.match(r"tests/", src_name) is not None:
        color = "red"
    else:
        color = "black"
    return color


def main():
    # Parse command-line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file-name", type=str, action="store",
                        default="sources.pkl")
    parser.add_argument("-o", "--output", type=str, action="store",
                        default="dep.dot")
    args = parser.parse_args()

    sources = SourceTree()
    sources.load_cache(file_name=args.file_name)
    sources.write_dot(dot_name=args.output, color_func=color_func)


if __name__ == "__main__":
    main()