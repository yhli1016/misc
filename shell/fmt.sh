#! /bin/bash

# This is a wrapper to the fmt program of GNU coreutils.

fmt -w 80 $1 > $1.tmp
mv $1.tmp $1
