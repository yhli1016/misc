#! /usr/bin/env python

import os
import time


def get_ncore_total():
    with open("/proc/cpuinfo", "r") as cpuinfo_file:
        cpuinfo = cpuinfo_file.readlines()
        physid = [int(line.split()[3]) for line in cpuinfo 
                  if line.find("physical id") != -1]
        cores = [int(line.split()[3]) for line in cpuinfo
                 if line.find("cores") != -1]
        ncpu_phys = len(set(physid))
        ncore_per_phys = cores[0]
    return ncpu_phys * ncore_per_phys


def run_job(ncore):
    #command = "bash run_job.sh ncore"
    #command = "mpirun -np ncore vasp.mpi.5.3.3 &> log"
    os.system(command.replace("ncore", str(ncore)))


def main():
    program_list = ["vasp"]
    dt = 60
    maxloop = 86400 / dt * 3
    ncore_total = get_ncore_total()
    ncore_min = ncore_total
    ncore_max = ncore_total
    
    for i in range(int(maxloop)):
        process_list = os.popen("ps -e").readlines()
        core_status = [True for program in program_list 
                       for process in process_list
                       if process.find(program) != -1]
        ncore_idle = ncore_total - len(core_status)

        if ncore_idle >= ncore_min:
            run_job(min([ncore_idle, ncore_max]))
            break
        else:
            time.sleep(dt)


if __name__ == "__main__":
    main()
