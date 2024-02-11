"""Functions for manipulating cif files."""

import numpy as np
from toolkit.base.species import PERIODIC_TABLE


def get_lattice_constants(cif_content):
    """
    Extract lattice constants from content of cif file.

    :param cif_content: list of strings read from a cif file
    :return lattice_contants: list containg 6 real numbers
    """
    lattice_constants = []
    for const in ("a", "b", "c"):
        for line in cif_content:
            if line.find("_cell_length_%s" % const) != -1:
                lattice_constants.append(float(line.split()[1]))
                break
    for const in ("alpha", "beta", "gamma"):
        for line in cif_content:
            if line.find("_cell_angle_%s" % const) != -1:
                lattice_constants.append(float(line.split()[1]))
                break
    return lattice_constants


def get_atom_coordinates(cif_content, sort_atoms=False, sort_method="zval"):
    """
    Extract atomic symbols and coordinates from cif content.

    :param cif_content: list of strings read from a cif file
    :param sort_atoms: whether to sort atomic coordinates
    :param sort_method: sorting method
    :return atom_symbols: list of atomic symbols
    :return atom_coordinates: list of atomic coordinates
    """
    raw_atom_coordinates = []
    for line in cif_content:
        if line.find("Uiso") != -1:
            s = line.split()
            raw_atom_coordinates.append([s[1], float(s[2]), float(s[3]), float(s[4])])
    if sort_atoms:
        if sort_method == "zval":
            raw_atom_coordinates = sorted(raw_atom_coordinates, key=lambda x: PERIODIC_TABLE.index(x[0]))
        elif sort_method == "abc":
            raw_atom_coordinates = sorted(raw_atom_coordinates, key=lambda x: x[0])
        else:
            pass
    atom_symbols = [row[0] for row in raw_atom_coordinates]
    atom_coordinates = np.array([row[1:4] for row in raw_atom_coordinates])
    return atom_symbols, atom_coordinates


def get_atom_types(atom_symbols):
    """
    Get the list of unique atomic types from a list of atomic symbols.

    :param atom_symbols: list of atomic symbols
    :return atom_types: list of unique atomic types
    """
    atom_types = []
    for symbol in atom_symbols:
        if symbol not in atom_types:
            atom_types.append(symbol)
    return atom_types


def get_natom_per_type(atom_symbols):
    """
    Calculate the number of atoms of each type.

    :param atom_symbols: list of atomic symbols
    :return num_atom: dictionary with the key being atomic type and value being
                      the number of atoms of this type
    """
    atom_types = get_atom_types(atom_symbols)
    natom_per_type = dict()
    for atom_type in atom_types:
        natom_per_type[atom_type] = atom_symbols.count(atom_type)
    return natom_per_type