#! /usr/bin/env python
import os
from typing import Iterable, Dict

from mace import expand_file


def get_option(prompt: str, options: Iterable[str]) -> str:
    """
    Get one option from a list of options.

    :param prompt: prompt to the user
    :param options: available options
    :return: chosen option
    """
    while True:
        option = input(prompt)
        if option in options:
            break
    return option


def get_input() -> Dict[str, str]:
    """
    Get input parameters as a macro.

    NOTE: keys in lowercase and starting with '.' are reserved for
    special purposes. DO NOT use them as macro names in templates.

    :return: the macros with keys being names and values being contents
    """
    macro = dict()

    # Get the type of job
    macro[".job"] = get_option("\nInput type of job (opt/neb/bader): ",
                               ("opt", "neb", "bader"))

    # Get specific job information
    if macro[".job"] == "neb":
        num_image = int(input("\nInput number of transition states: "))
        macro["NIMAGE"] = num_image
        macro["DIR_TS"] = f"$(seq -f %02g 1 {num_image})"
        macro["DIR_TOT"] = f"$(seq -f %02g 0 {num_image + 1})"
        macro["NMAX"] = f"{num_image+1:02d}"

    # Get general job information
    macro["NAME"] = input("\nInput job name: ")
    macro["NCPU"] = input("\nInput number of cpu to use: ")
    macro["TIME"] = input("\nInput time limit (in hours): ")
    macro["MEM"] = 2048  # maximum memory, temporarily fixed to 2048 MB
    restart = get_option("\nInput whether to restart (yes/no): ",
                         ("yes", "no"))
    if restart == "yes":
        macro["ISTART"] = 1
    else:
        macro["ISTART"] = 0
    macro["RUN"] = input("\nInput number of run: ")
    macro[".incar"] = get_option("\nInput whether to write incar (yes/no): ",
                                 ("yes", "no"))
    return macro


def main():
    # Location of template files
    home = os.environ["HOME"]
    root_dir = f"{home}/soft/inpgen"
    script_dir = f"{root_dir}/euler"
    incar_dir = f"{root_dir}/incar"

    # Get arguments
    macro = get_input()

    # Generate input file
    job = macro[".job"]
    expand_file(macro, f"{script_dir}/run_{job}.m4", f"run_{job}.sh")
    if macro[".incar"] == "yes":
        expand_file(macro, f"{incar_dir}/INCAR_{job}.m4", "INCAR")

    # Notify the user
    print(f"\nScript written to 'run_{job}.sh'")
    if macro[".incar"] == "yes":
        print("\nInput written to 'INCAR'")


if __name__ == "__main__":
    main()
