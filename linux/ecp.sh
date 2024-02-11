#! /bin/bash

# Enhanced CoPy program
#
# This program is based on find and rsync, and acts as an enhanced cp program.
#
# Typically, this program is called before downloading data to localhost, in
# order to remove large files.
#
# Usage: ecp.sh SRC DESTINATION OPTIONS
# Example: ecp.sh 0 test "-size -10M"

if [ $# -ne 3 ]; then
    echo 
    echo "Usage: ecp.sh SRC DESTINATION OPTIONS"
    echo "Example: ecp.sh 0 test \"-size -10M\""
    echo
    exit -1
else
    src_dir=$1
    destination_dir=$2
    opt=$3

    find $src_dir $opt > filelist
    rsync -acP --files-from=filelist . $destination_dir
    wait
    rm filelist
fi
