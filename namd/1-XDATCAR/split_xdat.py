#! /usr/bin/env python

"""
Create a directory named 'POS' under top directory and copy XDATCAR 
and this program there before invoking.
"""

# read XDATCAR
with open("XDATCAR", "r") as xdatcar_file:
    xdatcar = xdatcar_file.readlines()

# extract the entries of each configuration
conf_entries = [i for i in range(len(xdatcar)) \
                if xdatcar[i].find("configuration") != -1]

# extract shared header of POSCAR and number of atoms
poscar_header = xdatcar[:conf_entries[0]]
nline_per_conf= conf_entries[1] - conf_entries[0]


# determine the format of the postfix in POSCAR_$n
nconf = len(conf_entries)
ndigit = len(str(nconf))
fmt = "%0" + str(ndigit) + "d"

# extract atomic positions of each configuration and write to POSCAR_$n
for entry in conf_entries:
    conf_id = int(xdatcar[entry].split()[2])
    nl0 = entry
    nl1 = nl0 + nline_per_conf
    with open("POSCAR_" + str(fmt % conf_id), "w") as poscar_file:
        for line in poscar_header:
            poscar_file.write(line)
        for line in xdatcar[nl0:nl1]:
            poscar_file.write(line)
