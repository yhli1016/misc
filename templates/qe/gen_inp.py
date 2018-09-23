#! /usr/bin/env python
"""
Copied from bpo/vcrelax/bop/1.mace2, demonstrating how to run mace as a
library and how to extract text from file without calling 'include'.
"""

import sys
sys.path.append("/home/yhli/soft/bin")
from mace import main


m = macro = dict()

for i in range(1, 11):
    m["prefix"] = i
    m["nat"] = i + 4
    with open("../cfg/%d.cfg" % i, "r") as cfg:
        content = cfg.readlines()
    nl0 = content.index("begin_vectors_ang\n")
    nl1 = content.index("end_vectors_ang\n")
    nl2 = content.index("begin_coordinates_frac\n")
    nl3 = content.index("end_coordinates_frac\n")
    m["cell"] = "".join(content[nl0+1:nl1])
    m["pos"] = "".join(content[nl2+1:nl3])
    main(macro, "template.tpl", "%d.in" % i)
