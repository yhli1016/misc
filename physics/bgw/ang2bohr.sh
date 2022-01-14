#! /bin/bash

# This program converts lattice vectors from angstroms to bohrs.
#
# This program is useful when generating input for gsphere.inp.
#
# Usage: This program reads from stdin and writes to stdout.

awk '{printf "%14.9f%14.9f%14.9f\n", $1*1.889726125, $2*1.889726125, $3*1.889726125}'
