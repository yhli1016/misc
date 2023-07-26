#! /usr/bin/env python
"""Script for generating cache for other scripts."""

import argparse

from fortlint import DependGraph


def main():
    # Parse command-line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file-name", type=str, action="store",
                        default="sources.pkl")
    parser.add_argument("paths", type=str, action="store", nargs="*")
    args = parser.parse_args()

    # Set default directory list
    directory_list = [".", "*"]
    # directory_list = [".", "*", "../futile-1.8/*", "../futile-1.8/wrappers/*"]
    directory_list.extend(args.paths)

    # Build the graph
    graph = DependGraph()
    for directory in directory_list:
        graph.parse_source_tree(directory)
    graph.resolve_dependencies()

    # Save cache
    graph.save_cache(file_name=args.file_name)


if __name__ == "__main__":
    main()
