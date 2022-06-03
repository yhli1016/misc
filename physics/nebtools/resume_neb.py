#! /usr/bin/env python
"""Resume NEB calculation from given step."""

import shutil

from ase.io import read, write

from check_conv import get_int


def backup(num_ts_images):
    """Back up results of last run."""
    last_run = get_int("run_neb.sh", "run")
    for i in range(1, num_ts_images+1):
        prefix = "%02d" % i
        shutil.copy(f"{prefix}/CONTCAR", f"{prefix}/CONTCAR_{last_run}")


def resume_from_traj():
    """Resume NEB from mep.traj."""
    mep_name = input("Input the name of mep to resume: ")
    ts_images = read(mep_name, index="1:-1")
    backup(len(ts_images))
    for i, image in enumerate(ts_images):
        prefix = "%02d" % (i + 1)
        write(f"{prefix}/CONTCAR", image, format="vasp", vasp5=True, direct=True)


def resume_from_poscar():
    """Resume NEB from CONTCAR"""
    num_ts_images = get_int("INCAR", "IMAGES")
    backup(num_ts_images)
    i_start = int(input("Input the number of run to resume: "))
    for i in range(1, num_ts_images+1):
        prefix = "%02d" % i
        shutil.copy(f"{prefix}/POSCAR_{i_start+1}", f"{prefix}/CONTCAR")


if __name__ == "__main__":
    resume_from_traj()
