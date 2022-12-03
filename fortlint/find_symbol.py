#! /usr/bin/env python
"""Script for finding the definitions or references of given symbol."""

import argparse

from fortlint import SourceTree


def main():
    # Parse command-line parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--ref", action="store_true", default=False)
    parser.add_argument("-f", "--file-name", type=str, action="store",
                        default="sources.pkl")
    parser.add_argument("symbol", type=str, action="store")
    args = parser.parse_args()

    # Load cache
    sources = SourceTree()
    sources.load_cache(file_name=args.file_name)

    # Search for the symbol and echo
    symbol = args.symbol
    if args.ref:
        candidates = sources.find_symbol(symbol, 'ref')
        print(f"References of '{symbol}' found in:")
        for item in candidates:
            print(f"{'':4s}{item}")
    else:
        candidates = sources.find_symbol(symbol, 'def')
        print(f"Definitions of '{symbol}' found in:")
        for item in candidates:
            print(f"{'':4s}{item}")


if __name__ == "__main__":
    main()
