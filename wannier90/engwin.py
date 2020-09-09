#! /usr/bin/env python
import sys
import xml.dom.minidom


Har2eV = 13.60569253 * 2


def get_energies(xml_name):
    dom = xml.dom.minidom.parse(xml_name)
    root = dom.documentElement
    eng_full = []
    if root.nodeName == "fleurOutput":
        eigenvalues = root.getElementsByTagName("eigenvalues")[-1]
        eks = eigenvalues.getElementsByTagName("eigenvaluesAt")
        eng_full = [[float(f) for f in ek.childNodes[0].data.split()] for ek in eks]
    elif root.nodeName == "qes:espresso":
        eigenvalues = root.getElementsByTagName("eigenvalues")
        eng_full = [[float(f) for f in ek.childNodes[0].data.split()] for ek in eigenvalues]
    else:
        raise RuntimeError("Unknown xml output")
    return eng_full


def main():
    xml_name = sys.argv[1]
    option = sys.argv[2]
    if option == "e":
        # Band index is counted from 1, NOT 0.
        band_index = int(sys.argv[3])
        eng_full = get_energies(xml_name)
        eng_selected =[eng[band_index-1] for eng in eng_full]
        emin = min(eng_selected) * Har2eV
        emax = max(eng_selected) * Har2eV
        print("emin = %f" % emin)
        print("emax = %f" % emax)
    elif option == "n":
        # Energies are in eV, not Hartree.
        emin = float(sys.argv[3]) / Har2eV
        emax = float(sys.argv[4]) / Har2eV
        eng_full = get_energies(xml_name)
        for ik, ek in enumerate(eng_full):
            num_bands = 0
            for eng in ek:
                if eng >= emin and eng <= emax:
                    num_bands += 1
            print("ik = %d, nbnd = %d" % (ik+1, num_bands))
    else:
        print("ERROR: unknown option '%s'" % option)


if __name__ == "__main__":
    main()
