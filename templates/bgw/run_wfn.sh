#! /bin/bash

function run_job()
{
    # add your scipts here
    local ncore=$1
    local MPIRUN="mpirun -np $ncore"
    local topdir=$(pwd)

    cd $topdir/01-scf
    $MPIRUN pw.x -in scf.in &> scf.out
    wait

    for pref in 02-wfn 03-wfnq 04-wfn_path 05-wfn_fi; do
    	cd $topdir/$pref
	    cp -r ../01-scf/*.save .
	    wait
	    $MPIRUN pw.x -in bands.in &> bands.out
	    wait
	    $MPIRUN pw2bgw.x -in p2b.in &> p2b.out
	    wait
	    rm -r *.save *.wfc*
	    wait
    done
}

function get_ncore_tot()
{
    local nphys=$(grep "physical id" /proc/cpuinfo | sort | uniq | wc -l)
    local ncore=$(grep "cores" /proc/cpuinfo | sort | uniq | awk '{print $4}')
    local ncore_tot=$(echo $nphys $ncore | awk '{print $1 * $2}')
    echo $ncore_tot
}

function get_ncore_idle()
{
    local prog_list=$1
    local ncore_tot=$2
    local ncore_busy=0
    local ncore_idle=0
    for prog in $prog_list; do
        local ncore_i=$(ps -e | grep "$prog" | wc -l)
        ncore_busy=$(echo $ncore_i $ncore_busy | awk '{print $1 + $2}')
    done
    ncore_idle=$(echo $ncore_tot $ncore_busy | awk '{print $1 - $2}')
    echo $ncore_idle
}

function on_watch()
{
    # programs to watch on
    local prog_list="pw yambo vasp"

    # time interval and maximum number of trial
    # default is one trial per 60s in 3 days
    local dt=60
    local nmax=4320

    # total number of cores on this host and min/max number of cores for this job
    # default is to use all the cores
    local ncore_tot=$(get_ncore_tot)
    local ncore_min=$ncore_tot
    local ncore_max=$ncore_tot
    
    # main loop
    for i in $(seq $nmax); do
        local ncore_idle=$(get_ncore_idle "$prog_list" $ncore_tot)
        if [ $ncore_idle -ge $ncore_min  ]; then
            if [ $ncore_idle -le $ncore_max ]; then
                run_job $ncore_idle
                break
            else
                run_job $ncore_max
                break
            fi
        else
            sleep $dt
        fi
    done
}

on_watch

