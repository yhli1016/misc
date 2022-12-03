#! /usr/bin/env python
"""
Script for finding the relevant nodes of given node in the dependency digraph.
"""

import argparse
import re

from fortlint import SourceTree


def main():
    # Parse command-line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--out", action="store_true", default=False)
    parser.add_argument("-f", "--file-name", type=str, action="store",
                        default="sources.pkl")
    parser.add_argument("node", type=str, action="store")
    args = parser.parse_args()

    # Load cache
    sources = SourceTree()
    sources.load_cache(file_name=args.file_name)

    # Remove leading "./" and trailing '.' if any
    node = args.node
    if node[:2] == "./":
        node = node[2:]
    while node[-1] == ".":
        node = node[:-1]

    # Remove suffix if any
    pattern = re.compile(r"^(\S+)\.\w+$", re.IGNORECASE)
    result = re.search(pattern, node)
    if result is not None:
        node = result.group(1)

    # Search for the nodes and echo
    if args.out:
        candidates = sources.find_relevant_nodes(node, 'out')
        print(f"'{node}' -> {len(candidates)} nodes:")
        for item in candidates:
            print(f"{'':4s}{item}")
    else:
        candidates = sources.find_relevant_nodes(node, 'in')
        print(f"{len(candidates)} nodes -> '{node}':")
        for item in candidates:
            print(f"{'':4s}{item}")


if __name__ == "__main__":
    main()
