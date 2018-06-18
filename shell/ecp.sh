#! /bin/bash

# Ehanced CoPy program
#
# This program is based on find and rsync, and acts as an enhanced cp program.
#
# Typically, this program is called before downloading data to localhost, in
# order to remove largefiles.
#
# Usage: ecp.sh SRC DEST OPTIONS
# Example: ecp.sh 0 test "-size -10M"

if [ $# -ne 3 ]; then
    echo 
    echo "Usage: ecp.sh SRC DEST OPTIONS"
    echo "Example: ecp.sh 0 test \"-size -10M\""
    echo
    exit -1
else
    srcdir=$1
    destdir=$2
    opt=$3

    find $srcdir $opt > filelist
    rsync -acP --files-from=filelist . $destdir
    wait
    rm filelist
fi
