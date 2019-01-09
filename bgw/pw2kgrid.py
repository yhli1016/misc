#! /usr/bin/env python
"""
This program generates input files for kgrid.x from the input of pw.x.

Usage: pw2kgrid.py inp out nx ny nz dx dy dz qx qy qz gx gy gz
    inp: input of pwscf
    nx/ny/nz: dimension of the size of MP k-grid,
    dx/dy/dz: shift of MP k-grid
    qx/qy/qz: shift for epsilon and absorption calculation
    gx/gy/gz: dimension of FFT grid determined from pwscf calculation or
              gsphere.py

 Notes:
 (1) Only scf.in, nscf.in and bands.in are acceptable.
 (2) Lattice vectors must be specified manually, i.e. ibrav = 0.
 (3) Atomic positions must be fractional coordinates, i.e.
     ATOMIC_POSITIONS {crystal}.
 (4) Specifying FFT grid size is mandatory, otherwise the following GW/BSE
     calculations are likely to CRASH due to symmetry reasons.
"""

import sys


def read_inp(inp_name):
    with open(inp_name) as inp_file:
        raw_content = inp_file.readlines()
    content = [line for line in raw_content
               if line.find("#") == -1 and line.find("!") == -1]
    return content


def check_fractional(content):
    is_fractional = False
    for line in content:
        if (line.find("ATOMIC_POSITIONS") != -1
            and line.find("crystal") != -1
            and len(line.split()) == 2):
            is_fractional = True
            break
    return is_fractional


def check_ibrav(content):
    is_brav_zero = False
    for line in content:
        if (line.find("ibrav") != -1
            and len(line.split()) == 3
            and int(line.split()[2]) == 0):
            is_brav_zero = True
            break
    return is_brav_zero


def extract_nat(content):
    nat = 0
    for line in content:
        if (line.find("nat") != -1
            and len(line.split()) == 3):
            nat = int(line.split()[2])
            break
    return nat


def extract_species(content):
    species = []
    for line in content:
        if (line.find("UPF") != -1
            or line.find("upf") != -1
            and len(line.split()) == 3):
            species.append(line.split()[0])
    return species


def extract_vectors(content):
    vector = []
    for line in content:
        if line.find("CELL_PARAMETERS") != -1:
            nl_start = content.index(line) + 1
            for i in range(nl_start, nl_start+3):
                line_split = content[i].split()
                vector.append([float(line_split[0]),
                               float(line_split[1]),
                               float(line_split[2])])
    return vector


def extract_coordinates(content):
    symbols = []
    coordinates  = []
    nl_start = 0
    for line in content:
        if line.find("ATOMIC_POSITIONS") != -1:
            nl_start = content.index(line) + 1
            break
    nat = extract_nat(content)
    nl_end = nl_start + nat
    for i in range(nl_start, nl_end):
        line_split = content[i].split()
        symbols.append(line_split[0])
        coordinates.append([float(line_split[1]),
                            float(line_split[2]),
                            float(line_split[3])])
    return symbols, coordinates


def frac2cart(frac_coord, vectors):
    cart_coord = []
    for i in frac_coord:
        x0 = i[0]
        y0 = i[1]
        z0 = i[2]
        x1 = x0 * vectors[0][0] + y0 * vectors[1][0] + z0 * vectors[2][0]
        y1 = x0 * vectors[0][1] + y0 * vectors[1][1] + z0 * vectors[2][1]
        z1 = x0 * vectors[0][2] + y0 * vectors[1][2] + z0 * vectors[2][2]
        cart_coord.append([x1, y1, z1])
    return cart_coord


def main():
    argv = sys.argv
    inp_name, out_name = argv[1], argv[2]
    nx, ny, nz = int(argv[3]), int(argv[4]), int(argv[5])
    dx, dy, dz = float(argv[6]), float(argv[7]), float(argv[8])
    qx, qy, qz = float(argv[9]), float(argv[10]), float(argv[11])
    gx, gy, gz = int(argv[12]), int(argv[13]), int(argv[14])

    inp_content = read_inp(inp_name)

    # check PWSCF input
    if not check_fractional(inp_content):
        print("ERROR: Atomic positions are not fractional!")
        sys.exit(-1)

    if not check_ibrav(inp_content):
        print("ERROR: ibrav != 0")
        sys.exit(-1)

    # extract other information from PWSCF input
    nat = extract_nat(inp_content)
    species = extract_species(inp_content)
    vectors = extract_vectors(inp_content)
    symbols, frac_coord = extract_coordinates(inp_content)
    cart_coord = frac2cart(frac_coord, vectors)

    # Output
    with open(out_name, "w") as out_file:
        # write kgrid
        out_file.write("%8d%8d%8d\n" % (nx, ny, nz))
        out_file.write("%8.3f%8.3f%8.3f\n" % (dx, dy, dz))
        out_file.write("%8.3f%8.3f%8.3f\n" % (qx, qy, qz))
        out_file.write("\n")

        # write vectors
        for i in range(3):
            for j in range(3):
                out_file.write("%16.9f" % vectors[i][j])
            out_file.write("\n")
        out_file.write("\n")

        # write atoms
        out_file.write("%4d\n" % nat)
        for i in range(nat):
            out_file.write("%4d" % (species.index(symbols[i]) + 1))
            for j in range(3):
                out_file.write("%16.9f" % cart_coord[i][j])
            out_file.write("\n")
        out_file.write("\n")

        # write FFT grid and symmetry options
        out_file.write("%4d%4d%4d\n" % (gx, gy, gz))
        out_file.write("%8s\n%8s\n" % (".false.", ".false."))


if __name__ == "__main__":
   main()
