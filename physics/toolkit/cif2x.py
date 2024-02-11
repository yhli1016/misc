#! /usr/bin/env python
import argparse
import toolkit.base.cif as cif
import toolkit.base.lattice as lattice
import toolkit.base.interface as interface


def main():
    # Parse CLI parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("cif_name", type=str,  action="store")
    parser.add_argument("out_name", type=str, action="store")
    parser.add_argument("-f", "--format", type=str, action="store", default="qe")
    parser.add_argument("-s", "--shift-center", action="store_true", default=False)
    parser.add_argument("-d", "--dimension", type=int, action="store", default=3)
    parser.add_argument("-c", "--center", type=float,  action="store", default=0.0)
    parser.add_argument("-S", "--sort-atoms", action="store_true", default=False)
    parser.add_argument("-m", "--sort-method", type=str, action="store", default="zval")
    args = parser.parse_args()

    # Read and parse cif file
    with open (args.cif_name, "r") as cif_file:
        cif_content = cif_file.readlines()
    lattice_constants = cif.get_lattice_constants(cif_content)

    # Calculate lattice vectors and shift
    lattice_vectors = lattice.gen_lattice_vectors(lattice_constants)

    # Treat on atomic types and coordinates
    if args.format == "vasp":
        args.sort_atoms = True
    atom_symbols, atom_coordinates = cif.get_atom_coordinates(cif_content, sort_atoms=args.sort_atoms,
                                                              sort_method=args.sort_method)
    natom_total = len(atom_symbols)
    natom_per_type = cif.get_natom_per_type(atom_symbols)
    ntype = len(natom_per_type)
    if args.shift_center is True:
         atom_coordinates = lattice.shift_geo_center(atom_coordinates, args.dimension, args.center)

    # Output
    print("%-16s%4d" % ("Number of types:", ntype))
    print("%-16s%4d" % ("Number of atoms:", natom_total))
    print("Number of atoms for each type:")
    for key, value in natom_per_type.items():
        print("%4s%4d" % (key, value))
    with open(args.out_name, "w") as out_file:
        if args.format == "qe":
            interface.write_qe(lattice_vectors, atom_symbols, atom_coordinates, out_file)
        elif args.format == "siesta":
            interface.write_siesta(lattice_vectors, atom_symbols, atom_coordinates, out_file)
        elif args.format == "vasp":
            interface.write_vasp(lattice_vectors, atom_symbols, atom_coordinates, out_file)
        else:
            raise NotImplementedError("Format '%s' not implemented yet." % args.format)


if __name__ == "__main__":
   main()
