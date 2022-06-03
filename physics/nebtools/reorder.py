#! /usr/bin/env python
"""Reorder atoms in POSCAR according to reference structure."""

import argparse

import numpy as np
from ase import Atoms
from ase.io import read, write


# Parse command-line parameters
parser = argparse.ArgumentParser()
parser.add_argument("poscar_ref", type=str, action="store")
parser.add_argument("poscar_chk", type=str, action="store")
parser.add_argument("--tol", "-t",   type=float, action="store", default=0.5)
args = parser.parse_args()

# Read structures
struct_ref = read(args.poscar_ref, index=-1)
struct_chk = read(args.poscar_chk, index=-1)

# Match atoms
struct_new = []
for idx_ref, atom_ref in enumerate(struct_ref):
    # Search for matching atoms according to positions.
    matched = False
    for atom_chk in struct_chk:
        if atom_ref.symbol == atom_chk.symbol and \
                np.linalg.norm(atom_ref.position - atom_chk.position) <= args.tol:
            struct_new.append(atom_chk)
            matched = True

    # If not found, match them manually.
    if not matched:
        idx_match = int(input(f"{idx_ref}->: "))
        atom_chk = struct_chk[idx_match]
        if atom_ref.symbol == atom_chk.symbol:
            struct_new.append(atom_chk)
        else:
            print("ERROR: atomic species does not match!")

# Output
struct_new = Atoms(struct_new)
struct_new.cell = struct_chk.cell
write("POSCAR_new", struct_new, format="vasp", vasp5=True, direct=True)
