#!/usr/bin/env python
import os


pot_dir = "/home/apps/vasp/potcar/paw_PBE"

with open("POSCAR", "r") as poscar:
    content = poscar.readlines()
    symbols = content[5].split()

# POTCAR
os.system(f"cp {pot_dir}/{symbols[0]}/POTCAR .")
for s in symbols[1:]:
    os.system(f"cat {pot_dir}/{s}/POTCAR >> POTCAR")

# vdw_kernel.bindat
os.system(f"cat {pot_dir}/../vdw_kernel.bindat > vdw_kernel.bindat")