#! /usr/bin/env python
"""Script for generating cache for other scripts."""

import argparse

from fortlint import SourceTree


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

    # Build the source tree
    source_tree = SourceTree()
    for directory in directory_list:
        source_tree.parse_source_tree(directory)
    source_tree.resolve_dependencies()

    # Save cache
    source_tree.save_cache(file_name=args.file_name)


if __name__ == "__main__":
    main()
