#! /usr/bin/env python
"""Script for generating the dependency digraph using Graphviz."""

import re
import argparse
import os

from fortlint import DependGraph


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
                        default="dep.svg")
    args = parser.parse_args()

    graph = DependGraph()
    graph.load_cache(file_name=args.file_name)
    dot_name = args.output.split(".")[0] + ".dot"
    graph.write_dot(dot_name=dot_name, color_func=color_func)
    os.system(f"dot -Tsvg -o{args.output} {dot_name}")


if __name__ == "__main__":
    main()
