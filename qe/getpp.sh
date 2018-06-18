#! /bin/bash

# Last modified on 2015-12-17
#
# This program downloads pseudo potential files from QE official host using either
# curl of wget.
#
# Part of this program are adopted from check-pw.x.j and environment_variables of
# QE source code.
#
# Usage: getpp.sh X.xxx-xxx.UPF

# remote host of QE
NETWORK_PSEUDO=http://www.quantum-espresso.org/wp-content/uploads/upf_files

# wget or curl needed to downloaded PP from web site
WGET="wget -O"
#WGET="curl -o"

# download the pseudo potential file
ppfile=$1
$WGET $ppfile $NETWORK_PSEUDO/$ppfile
