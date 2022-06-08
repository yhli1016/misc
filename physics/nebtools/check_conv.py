#! /usr/bin/env python
"""Check convergence of ionic steps from OSZICAR."""

import re
import sys


def get_int(file_name, pattern):
    """
    Get value of pattern in form of 'pattern = value' from file.

    :param str file_name: name of the file
    :param str pattern: name of pattern
    :return: value: int
    """
    value = None
    with open(file_name, "r") as inf:
        content = inf.readlines()
        for line in content:
            result = re.search(r"%s( )?=( )?\d+" % pattern, line)
            if result is not None:
                value = int(result.group().split("=")[1])
    if value is None:
        raise RuntimeError(f"{pattern} not found")
    return value


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
            print("WARNING: ionic step %4d not completed" % (len(start_lines)))

        # Check convergence for each ionic step
        for i, nl in enumerate(end_lines):
            num_step = nl - start_lines[i] - 1
            if num_step >= nelm:
                print("WARNING: ionic step %4d not converged" % i)
            else:
                print("ionic step %4d converged in %4d iterations " % (i, num_step))


if __name__ == "__main__":
    main()
