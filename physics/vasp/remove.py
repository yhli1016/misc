#!/usr/bin/env python
import re
import os
from typing import List


class FileSet:
    def __init__(self, prefix: List[str] = None,
                 suffix: List[str] = None,
                 head: List[str] = None,
                 tail: List[str] = None,
                 body: List[str] = None,
                 word: List[str] = None) -> None:
        self.patterns = []
        if prefix is not None:  # prefix.*
            self.patterns.extend([re.compile(f"^{_}\.\w+$") for _ in prefix])
        if suffix is not None:  # *.suffix
            self.patterns.extend([re.compile(f"^\w+\.{_}$") for _ in suffix])
        if head is not None:  # head*
            self.patterns.extend([re.compile(f"^{_}.*$") for _ in head])
        if tail is not None:  # *tail
            self.patterns.extend([re.compile(f"^.*{_}$") for _ in tail])
        if body is not None:  # *body*
            self.patterns.extend([re.compile(f"^.*{_}.*$") for _ in body])
        if word is not None:  # exact match word
            self.patterns.extend([re.compile(f"^{_}$") for _ in word])

    def __add__(self, other):
        result = FileSet()
        result.patterns = self.patterns + other.patterns
        return result

    def include(self, file_name: str) -> bool:
        included = False
        for p in self.patterns:
            if re.search(p, file_name) is not None:
                included = True
                break
        return included


class ProgressBar:
    def __init__(self, num_tasks: int, num_scales: int = 10) -> None:
        self._num_tasks = num_tasks
        self._num_scales = num_scales
        self._scale_unit = self._num_tasks // self._num_scales
        self._next_scale = self._scale_unit
        self._counter = 0

    def count(self) -> None:
        self._counter += 1
        if self._counter >= self._next_scale:
            self._next_scale += self._scale_unit
            percent = 100.0 * self._counter / self._num_tasks
            print(f"finished {self._counter:6d} of {self._num_tasks:6d} "
                  f"({percent:6.2f}%)")


def get_all_files() -> List[str]:
    all_files = os.popen("find . -type f").readlines()
    all_files = [_.rstrip("\n") for _ in all_files]
    all_files = [_ for _ in all_files if os.path.isfile(_)]
    return all_files


def remove_files() -> None:
    words = ["INCAR", "KPOINTS"]
    heads = ["POSCAR", "CONTCAR", "OUTCAR"]
    suffixes = ["py", "txt"]
    bodies = ["README"]
    file_set = FileSet(word=words, head=heads, suffix=suffixes, body=bodies)

    all_files = get_all_files()
    files_to_remove = [_ for _ in all_files
                       if not file_set.include(_.split("/")[-1])]

    print(f"Total number of files: {len(all_files)}")
    print(f"Number of files to remove: {len(files_to_remove)}")
    prog_bar = ProgressBar(len(files_to_remove), 100)
    for f in files_to_remove:
        print(f"Removing {f}")
        os.system(f"rm -rf {f}")
        prog_bar.count()


if __name__ == "__main__":
    remove_files()
