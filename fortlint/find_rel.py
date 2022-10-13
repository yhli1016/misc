#! /usr/bin/env python
"""
Script for finding the relevant nodes of given node in the dependency digraph.
"""

import argparse

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

    # Search for the symbol and echo
    node = args.node
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
