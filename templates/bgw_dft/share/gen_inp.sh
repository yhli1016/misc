#! /bin/bash

for directory in 01-scf 02-wfn 03-wfnq 04-wfn_path 05-wfn_fi; do
    cd ../$directory

    # Copy pseudo-potential files
    cp ../share/*.UPF .

    # Call mace.py to generate input files
    for prefix in scf bands p2b; do
        if [ -f $prefix.tpl ]; then
            mace.py -p ../share $prefix.tpl $prefix.in
        fi
    done
done
