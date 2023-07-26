#! /usr/bin/env python
"""Script for generating the call digraph using Graphviz."""

import argparse
import os

from fortlint import Rules, CallGraph


class MyRules(Rules):
    def __init__(self):
        super().__init__()
        exclude = {"hop_ind", "kq_map", "x",  "q_point", "omegas", "kmesh",
                   "eng"}
        self._exclude = self._exclude.union(exclude)


def main():
    # Parse command-line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file-name", type=str, action="store",
                        default="sources.pkl")
    parser.add_argument("-o", "--output", type=str, action="store",
                        default="call.svg")
    args = parser.parse_args()

    rules = MyRules()
    sources = CallGraph(rules)
    sources.parse_source_tree(".")
    dot_name = args.output.split(".")[0] + ".dot"
    sources.write_dot(dot_name=dot_name)
    os.system(f"dot -Tsvg -o{args.output} {dot_name}")


if __name__ == "__main__":
    main()
