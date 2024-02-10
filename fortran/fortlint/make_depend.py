#! /usr/bin/env python
"""Script for generating make.depend for Makefile."""

import argparse

from fortlint import DependGraph


def main():
    # Parse command-line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file-name", type=str, action="store",
                        default="sources.pkl")
    parser.add_argument("-o", "--output", type=str, action="store",
                        default="make.depend")
    args = parser.parse_args()

    graph = DependGraph()
    graph.load_cache(file_name=args.file_name)
    graph.write_make(args.output)


if __name__ == "__main__":
    main()
