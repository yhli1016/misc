#! /bin/bash

pkg=$(python3 ./get_old_oneapi.py)
zypper remove $pkg
