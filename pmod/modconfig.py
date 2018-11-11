#! /usr/bin/env python
import sys 


def add_mod(mod_name, mod_type, dest, dep=None, cmd=None):
    # Create the mod if not exist
    if mod_name not in avail_mods.keys():
        avail_mods[mod_name] = dict()

    # Add essential keys 
    mod = avail_mods[mod_name]
    if mod_type == "mod":
        mod["PATH"] = dest + "/bin"
        mod["LIBRARY_PATH"] = dest + "/lib"
        mod["LD_LIBRARY_PATH"] = dest + "/lib"
        mod["C_INCLUDE_PATH"] = dest + "/include"
        mod["CPLUS_INCLUDE_PATH"] = dest + "/include"
    elif mod_type == "path":
        mod["PATH"] = dest
    elif mod_type == "lib":
        mod["LIBRARY_PATH"] = dest
        mod["LD_LIBRARY_PATH"] = dest
    elif mod_type == "inc":
        mod["C_INCLUDE_PATH"] = dest
        mod["CPLUS_INCLUDE_PATH"] = dest
    elif mod_type == "py":
        mod["PYTHONPATH"] = dest
    else:
        print("echo \"Unknown mod_type %s\";" % mod_type)
        sys.stdout.flush()

    # Add optional keys
    if dep is not None:
        mod["DEPENDENCY"] = dep
    if cmd is not None:
        mod["CMD"] = cmd
        

def check_integrity():
    for mod_name, mod_env in avail_mods.items():
        if "DEPENDENCY" in mod_env.keys():
            for mod_dep in mod_env["DEPENDENCY"]:
                if mod_dep not in avail_mods.keys():
                    print("echo \"Module %s has unresolved dependency %s\";" %
                          (mod_name, mod_dep))


# Module configurations
avail_mods = dict()
prefix = "/home/yhli/soft"

# for shared libraries
add_mod("mkl-13.0.079", "mod", prefix+"/mkl-13.0.079")
add_mod("mkl-13.0.079", "lib", prefix+"/mkl-13.0.079/lib/intel64")
add_mod("fftw-3.3.4", "mod", prefix+"/fftw-3.3.4")
add_mod("hdf5-1.8.17", "mod", prefix+"/hdf5-1.8.17")
add_mod("libxc-4.2.3","mod", prefix+"/libxc-4.2.3")

# for openmpi
add_mod("openmpi-1.10.0", "mod", prefix+"/openmpi-1.10.0")

# for DFT and GW software
add_mod("qe-6.2", "path", prefix+"/qe-6.2/bin", 
        dep=["mkl-13.0.079", "openmpi-1.10.0", "hdf5-1.8.17"])

add_mod("qe-5.4.0", "path", prefix+"/espresso-5.4.0/bin", 
        dep=["mkl-13.0.079", "openmpi-1.10.0"])

add_mod("bgw-1.2.0", "path", prefix+"/bgw-1.2.0/bin", 
        dep=["mkl-13.0.079", "openmpi-1.10.0", "hdf5-1.8.17", "fftw-3.3.4"])

add_mod("elk-4.0.15", "path", prefix+"/elk-4.0.15/bin", 
        dep=["mkl-13.0.079", "openmpi-1.10.0"])

add_mod("wannier90-2.1.0", "path", prefix+"/wannier90-2.1.0",
        dep=["mkl-13.0.079", "openmpi-1.10.0"])

add_mod("vasp.5.4.1", "path", prefix+"/vasp.5.4.1/bin",
        dep=["mkl-13.0.079", "openmpi-1.10.0"])

add_mod("exciting.carbon", "path", prefix+"/exciting.carbon/bin",
        dep=["mkl-13.0.079", "openmpi-1.10.0"])

# for aipes
add_mod("anaconda3", "mod", prefix+"/anaconda3")
add_mod("amp-v0.6", "py", prefix+"/amp-v0.6", dep=["anaconda3"])
add_mod("aipes", "py", prefix+"/aipes", dep=["anaconda3", "amp-v0.6"])

# Final check
check_integrity()
