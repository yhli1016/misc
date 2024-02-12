#! /usr/bin/env python
"""Copy results from scratch to submit directory for resuming the job."""

import shutil

from nebtools.base import get_int, get_option, get_job_name


def main():
    job_type = get_option("Input type of job", ("opt", "neb"))
    job_name = get_job_name(f"run_{job_type}.sh")
    host_name = get_option("Input host", ("euler", "csuc"))
    if host_name == "euler":
        scratch = f"/cluster/scratch/zhangwenj/{job_name}"
    else:
        scratch = f"/scratch/wzhang/{job_name}"
    results = ["OUTCAR", "OSZICAR", "CONTCAR"]

    if job_type == "opt":
        for r in results:
            shutil.copy(f"{scratch}/{r}", r)
    elif job_type == "neb":
        num_image = get_int("INCAR", "IMAGES")
        for i in range(1, num_image+1):
            prefix = f"{i:02d}"
            for r in results:
                try:
                    shutil.copy(f"{scratch}/{prefix}/{r}", f"{prefix}/{r}")
                except FileNotFoundError as err:
                    if r == "OUTCAR":
                        print(f"WARNING: {scratch}/{prefix}/OUTCAR not found."
                              f" Copying OUTCAR.gz instead.")
                        shutil.copy(f"{scratch}/{prefix}/OUTCAR.gz",
                                    f"{prefix}/OUTCAR.gz")
                    else:
                        raise err
    else:
        raise ValueError(f"Unknown job type: {job_type}")


if __name__ == "__main__":
    main()
