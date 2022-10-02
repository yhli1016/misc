#! /usr/bin/env python
"""Copy results from scratch to submit directory for resuming the job."""

import shutil

from nebtools.check_conv import get_int


def get_option(prompt, options):
    """
    Get options from stdin.

    :param str prompt: prompt to the user
    :param List[str] options: options from which the user choose
    :return: str, selected option
    """
    while True:
        option = input(prompt)
        if option in options:
            break
    return option


def get_job_name(script_name):
    """
    Get job name from job script.

    :param str script_name: file name of job script
    :return: str, job name
    """
    job_name = None
    with open(script_name, "r") as script_file:
        script = script_file.readlines()
        for line in script:
            if line.find("#BSUB -J") != -1 or line.find("#SBATCH -J") != -1:
                job_name = line.split()[2].rstrip("\n")
            elif line.find("#SBATCH --job-name") != -1:
                job_name = line.split('=')[1].rstrip("\n")
            else:
                pass
    if job_name is None:
        raise RuntimeError("Job name not found")
    return job_name


def main():
    job_type = get_option("Input type of job (opt/neb): ", ("opt", "neb"))
    job_name = get_job_name(f"run_{job_type}.sh")
    scratch = f"/cluster/scratch/zhangwenj/{job_name}"
    results = ["OUTCAR", "OSZICAR", "CONTCAR"]

    if job_type == "opt":
        for r in results:
            shutil.copy(f"{scratch}/{r}", r)
    elif job_type == "neb":
        num_image = get_int("INCAR", "IMAGES")
        for i in range(1, num_image+1):
            prefix = "%02d" % i
            for r in results:
                shutil.copy(f"{scratch}/{prefix}/{r}", f"{prefix}/{r}")
    else:
        raise ValueError(f"Unknown job type: {job_type}")


if __name__ == "__main__":
    main()
