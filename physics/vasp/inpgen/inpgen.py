#! /usr/bin/env python
import os
from typing import Iterable, Dict

from mace import Mace


def get_input(prompt: str) -> str:
    """
    Get input from stdin.

    :param prompt: prompt to the user
    :return: input arguments
    """
    return input(f"\n{prompt}: ")


def get_option(prompt: str, options: Iterable[str]) -> str:
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


def get_macro() -> Dict[str, str]:
    """
    Get arguments as a macro.

    NOTE: keys in lowercase and starting with '.' are reserved for
    special purposes. DO NOT use them as macro names in templates.

    :return: the macros with keys being names and values being contents
    """
    macro = dict()

    # Get the type of job
    macro[".job"] = get_option("Input type of job",
                               options=("opt", "neb", "bader"))

    # Get specific job information
    if macro[".job"] == "neb":
        num_image = int(get_input("Input number of transition states"))
        macro["NIMAGE"] = num_image
        macro["DIR_TS"] = f"$(seq -f %02g 1 {num_image})"
        macro["DIR_TOT"] = f"$(seq -f %02g 0 {num_image + 1})"
        macro["NMAX"] = f"{num_image+1:02d}"

    # Get general job information
    macro["NAME"] = get_input("Input job name")
    macro["NCPU"] = get_input("Input number of cpu to use")
    macro["TIME"] = get_input("Input time limit (in hours)")
    macro["MEM"] = 2048  # maximum memory, temporarily fixed to 2048 MB
    restart = get_option("Input whether to restart",
                         options=("yes", "no"))
    if restart == "yes":
        macro["ISTART"] = 1
    else:
        macro["ISTART"] = 0
    macro["RUN"] = get_input("Input number of run")
    macro[".incar"] = get_option("Input whether to write INCAR",
                                 options=("yes", "no"))
    return macro


def main():
    # Location of template files
    home = os.environ["HOME"]
    root_dir = f"{home}/soft/inpgen"
    script_dir = f"{root_dir}/euler"
    incar_dir = f"{root_dir}/incar"

    # Get arguments
    macro = get_macro()
    m = Mace(macro)

    # Generate input file
    job = macro[".job"]
    m.expand_re(f"{script_dir}/run_{job}.sh", f"run_{job}.sh")
    if macro[".incar"] == "yes":
        Mace(macro).expand_re(f"{incar_dir}/INCAR_{job}", "INCAR")

    # Notify the user
    print(f"\nScript written to 'run_{job}.sh'")
    if macro[".incar"] == "yes":
        print("\nInput written to 'INCAR'")


if __name__ == "__main__":
    main()
