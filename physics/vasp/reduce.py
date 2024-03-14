#!/usr/bin/env python
import re
import os
from remove import get_all_files, FileSet, ProgressBar


def _reduce_outcar(outcar_name: str) -> None:
    with open(outcar_name, "r") as outcar:
        content = outcar.readlines()

    # Determine indices for each ionic and electronic step
    indices = dict()
    pattern = re.compile(r"^-+\s*Iteration\s*(\d+)\(\s*(\d+)\)")
    for nl, line in enumerate(content):
        result = re.search(pattern, line)
        if result is not None:
            ii, jj = int(result.group(1)), int(result.group(2))
            indices[(ii, jj)] = nl

    # Reduce OUTCAR if it has more than 1 ionic steps
    ionic_steps = set([_[0] for _ in indices.keys()])
    if len(ionic_steps) > 1:
        step_min = min(ionic_steps)
        step_max = max(ionic_steps)
        step_1st_nl = indices[(step_min, 1)]
        step_max_nl = indices[(step_max-1, 1)]
        outcar_reduced = f"{outcar_name}_re"
        with open(outcar_reduced, "w") as outcar:
            for line in content[:step_1st_nl]:
                outcar.write(line)
            for line in content[step_max_nl:]:
                outcar.write(line)
        os.system(f"rm -rf {outcar_name}")


def reduce_outcar(outcar_name: str) -> None:
    try:
        print(f"Reducing {outcar_name}")
        _reduce_outcar(outcar_name)
    except:
        print(f"Skipping {outcar_name}")


def reduce_all(file_set: FileSet) -> None:
    outcars = [_ for _ in get_all_files()
               if file_set.include(_.split("/")[-1])]
    prog_bar = ProgressBar(len(outcars), 100)
    for f in outcars:
        if f.split(".")[-1] == "gz":
            outcar_name = f.replace(".gz", "")
            if os.path.isfile(outcar_name):
                reduce_outcar(outcar_name)
                prog_bar.count()
            else:
                print(f"Unzipping {f}")
                os.system(f"gunzip {f}")
                reduce_outcar(outcar_name)
                prog_bar.count()
        else:
            outcar_name = f
            reduce_outcar(outcar_name)
            prog_bar.count()


def zip_all(file_set: FileSet) -> None:
    outcars = [_ for _ in get_all_files()
               if file_set.include(_.split("/")[-1])]
    prog_bar = ProgressBar(len(outcars), 100)
    for f in outcars:
        print(f"Zipping {f}")
        os.system(f"gzip {f}")
        prog_bar.count()


def main():
    heads = ["OUTCAR"]
    file_set = FileSet(head=heads)
    reduce_all(file_set)
    zip_all(file_set)


if __name__ == "__main__":
    main()
