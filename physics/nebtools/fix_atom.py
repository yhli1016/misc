#! /usr/bin/env python
"""Fix given atoms."""

from ase.io import read, write
from ase.constraints import FixAtoms

from nebtools.neb import select_atoms


def main():
    poscar_in = "CONTCAR"
    poscar_out = "CONTCAR_fixed"
    selected_atoms = [94, 111, 41, 61, 60, 65, 24]
    fix_selected_atoms = True
    offset = 0

    image = read(poscar_in, index=-1)
    if fix_selected_atoms:
        fixed_atoms = selected_atoms
    else:
        fixed_atoms = select_atoms(image, exclude_atoms=selected_atoms,
                                   offset=offset)
    image.set_constraint(FixAtoms(fixed_atoms))
    write(poscar_out, image, format="vasp", vasp5=True, direct=True)


if __name__ == "__main__":
    main()
