#! /usr/bin/env python3
import os
import re


def get_versioned_packages():
    command = "zypper search -i oneapi | grep oneapi | awk '{print $3}'"
    pattern = re.compile(r"^intel\-oneapi\-([a-zA-Z]+\-)+[0-9\.]+$")
    full_packages = os.popen(command).readlines()
    versioned_packages = []
    for pkg in full_packages:
        result = re.search(pattern, pkg)
        if result is not None:
            versioned_packages.append(result.group().split("-"))
    return versioned_packages


def get_old_packages(versioned_packages):
    # Convert versioned packages into dictionary
    package_dict = {}
    for pkg in versioned_packages:
        name = tuple(pkg[:-1])
        ver = pkg[-1]
        try:
            package_dict[name].append(ver)
        except KeyError:
            package_dict[name] = [ver]

    # Get old packages
    old_packages = []
    for name, ver in package_dict.items():
        sorted_verions = sorted(ver)
        for old_ver in sorted_verions[:-1]:
            old_packages.append("".join([f"{_}-" for _ in name]) + old_ver)
    return old_packages


def main():
    versioned_packages = get_versioned_packages()
    old_packages = get_old_packages(versioned_packages)
    for pkg in old_packages:
        print(f"{pkg} ", end="")


if __name__ == "__main__":
    main()
