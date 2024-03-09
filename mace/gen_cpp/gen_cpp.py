#! /usr/bin/env python
import re
from collections import defaultdict
from mace import Mace


def parse_template(file_name):
    # Read raw content
    with open(file_name, "r") as infile:
        content = infile.readlines()

    # Determine starting and ending line numbers
    begin_pattern = re.compile(r"^\s*\/\s*\*\s*begin\s+(\w+)\s*\*\s*\/$", re.I)
    end_pattern = re.compile(r"^\s*\/\s*\*\s*end\s+(\w+)\s*\*\s*\/$", re.I)
    entries = defaultdict(dict)
    for i, line in enumerate(content):
        result = re.search(begin_pattern, line)
        if result is not None:
            entries[result.group(1)]["begin"] = i
        result = re.search(end_pattern, line)
        if result is not None:
            entries[result.group(1)]["end"] = i

    # Check if any item has missing line numbers
    for key1, value in entries.items():
        for key2 in ("begin", "end"):
            if not key2 in value.keys():
                raise RuntimeError(f"segment {key1} has no {key2}")

    # Split content to yield the template
    template = dict()
    for key, value in entries.items():
        nl0, nl1 = value["begin"], value["end"]
        template[key] = content[nl0:nl1+1]
    return template


def main():
    out_file = open("out.cpp", "w")

    mace = Mace()
    template = parse_template("scalar.i")
    dtypes = [("int", "MPI_INT"), ("double", "MPI_DOUBLE"), ("std::complex<double>", "MPI_COMPLEX")]
    for value in template.values():
        for t in dtypes:
            mace["scalar_type"] = t[0]
            mace["mpi_type"] = t[1]
            out_file.writelines(mace.expand_lines(value))
            out_file.write("\n")

    mace.clear()
    template = parse_template("matrix.i")
    dtypes = [("MatrixXi", "MPI_INT"), ("MatrixXd", "MPI_DOUBLE"), ("MatrixXcd", "MPI_COMPLEX")]
    for value in template.values():
        for t in dtypes:
            mace["eigen_type"] = t[0]
            mace["mpi_type"] = t[1]
            out_file.writelines(mace.expand_lines(value))
            out_file.write("\n")

    out_file.close()


if __name__ == "__main__":
    main()
