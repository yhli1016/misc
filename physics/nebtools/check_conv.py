#! /usr/bin/env python
"""Check convergence of ionic steps from OSZICAR."""

import sys

from nebtools.base import get_int


def main():
    # Get arguments
    oszicar_name = sys.argv[1]

    # Get NELM from INCAR
    try:
        nelm = get_int("INCAR", "NELM")
    except RuntimeError:
        nelm = 60

    # Analyze OSZICAR
    with open(oszicar_name, "r") as oszicar:
        content = oszicar.readlines()

        # Get line numbers of starting and ending lines
        start_lines = []
        end_lines = []
        for i, line in enumerate(content):
            if line.find("rms") != -1:
                start_lines.append(i)
            elif line.find("E0") != -1:
                end_lines.append(i)

        # Check if the calculation is complete
        if len(start_lines) != len(end_lines):
            print(f"WARNING: ionic step {len(start_lines):4d} not completed")

        # Check convergence for each ionic step
        for i, nl in enumerate(end_lines):
            num_step = nl - start_lines[i] - 1
            if num_step >= nelm:
                print(f"WARNING: ionic step {i:4d} not converged")
            else:
                print(f"ionic step {i:4d} converged in {num_step:4d} iterations")


if __name__ == "__main__":
    main()
