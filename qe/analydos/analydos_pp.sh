#! /bin/bash

if [ $# != 3 ]; then
    echo "USAGE: analydos_pp.sh dosfn bndfn outfn" 
    exit -1
else
    dosfn=$1
    bndfn=$2
    outfn=$3
   
    # The output of analy*dos.x (dosfn) may contain less kpoints and bands
    # than that in band structure data (bndstr.pd). So here we need to know
    # the lower and upper bounds of kpoint and band indices in dosfn.
    ikmin=$(awk 'NR>1 {print $1}' $dosfn | sort -g | uniq | head -1)
    ikmax=$(awk 'NR>1 {print $1}' $dosfn | sort -g | uniq | tail -1)
    ibmin=$(awk 'NR>1 {print $2}' $dosfn | sort -g | uniq | head -1)
    ibmax=$(awk 'NR>1 {print $2}' $dosfn | sort -g | uniq | tail -1)
    
    # clear outfn
    echo "# kpath energy projection" > $outfn

    # x.tmp and y.tmp are extracted from bndstr.pd which contains all the
    # kpoints and bands. So it is safe to filt data using NR.
    #
    # z.tmp is extracted from dosfn, where NR is not applicable. In stead
    # we use $1 and $2 to filt the data.
    awk '(NR>='$ikmin') && (NR<='$ikmax') {print $1}' $bndfn > x.tmp
    for ib in $(seq $ibmin $ibmax); do
        awk '(NR>='$ikmin') && (NR<='$ikmax') {print $('$ib'+1)}' $bndfn > y.tmp
        awk '(NR>1) && ($1>='$ikmin') && ($1<='$ikmax') && ($2=='$ib') {print $3}' $dosfn > z.tmp

        # Add to outfn; A blank line is essential for gnuplot.
        paste x.tmp y.tmp z.tmp | expand -t 8 >> $outfn
        echo >> $outfn
    done
fi 

rm *.tmp
