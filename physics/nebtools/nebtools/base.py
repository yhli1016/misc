"""Common utilities."""

import re
from typing import Tuple


def get_int(file_name: str, key: str) -> int:
    """
    Get value of pattern in form of 'key = value' from file.

    :param file_name: name of the file
    :param key: name of pattern
    :return: value
    """
    value = None
    pattern = r"^\s*%s\s*=\s*[\.\-\d]+\s*$" % key
    pattern = re.compile(pattern, re.I)
    with open(file_name, "r") as inf:
        content = inf.readlines()
        for line in content:
            result = re.search(pattern, line)
            if result is not None:
                value = int(result.group().split("=")[1])
    if value is None:
        raise RuntimeError(f"{key} not found")
    return value


def get_input(prompt: str) -> str:
    """
    Get input from stdin.

    :param prompt: prompt to the user
    :return: input arguments
    """
    return input(f"\n{prompt}: ")


def get_option(prompt: str, options: Tuple[str, ...]) -> str:
    """
    Get one option from a list of options.

    :param prompt: prompt to the user
    :param options: available options
    :return: chosen option
    """
    full_prompt = f"{prompt} ("
    for s in options:
        full_prompt += f"{s}/"
    full_prompt = full_prompt[:-1] + ")"
    while True:
        option = get_input(full_prompt)
        if option in options:
            break
    return option


def get_job_name(script_name: str) -> str:
    """
    Get job name from job script.

    :param script_name: file name of job script
    :return: job name
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
