#! /bin/env python
import sys
import os
import argparse
from modconfig import avail_mods


def print_cmd(cmd):
    print(cmd)
    sys.stdout.flush()


def avail_mod():
    for mod_name in avail_mods.keys():
        fmt_str = "\\e[32m\\t%s\\e[0m" % mod_name
        print_cmd("echo -e \"%s\";" % fmt_str)


def list_mod():
    for mod_name, mod_env in avail_mods.items():
        # Determine the number of effective keys
        num_key = len(mod_env.keys())
        if "DEPENDENCY" in mod_env.keys():
            num_key -= 1
        if "CMD" in mod_env.keys():
            num_key -= 1

        # Check if the environmental variables have been correctly set
        num_set_env = 0
        for key, value in mod_env.items():
            if key not in ("DEPENDENCY", "CMD") and key in os.environ.keys():
                if os.environ[key].find(value) != -1:
                    num_set_env += 1
        if num_set_env == num_key:
            fmt_str = "\\e[32m\\t[L]\\t%s\\e[0m" % mod_name
            print_cmd("echo -e \"%s\";" % fmt_str)
        elif num_set_env in range(1, num_key):
            fmt_str = "\\e[33m\\t[B]\\t%s\\e[0m" % mod_name
            print_cmd("echo -e \"%s\";" % fmt_str)
        else:
            fmt_str = "\\e[31m\\t[U]\\t%s\\e[0m" % mod_name
            print_cmd("echo -e \"%s\";" % fmt_str)

                
def load_mod(mod_name):
    if mod_name in avail_mods.keys():
        mod_env = avail_mods[mod_name]
        # Execute initialization commands
        if "CMD" in mod_env.keys():
            for cmd in mod_env["CMD"]:
                print_cmd("%s;" % cmd)
        # Set environmetal variables
        for key in mod_env.keys():
            new_env = mod_env[key]
            if key not in ("DEPENDENCY", "CMD") and key in os.environ.keys():
                if os.environ[key].find(new_env) == -1:
                    new_env = new_env + ":" + os.environ[key]
                    print_cmd("export %s=%s;" % (key, new_env))
            elif key not in ("DEPENDENCY", "CMD"):
                print_cmd("export %s=%s;" % (key, new_env))
    else:
        print_cmd("echo \"Module %s not found\";" % mod_name)
        sys.exit(-1)


def unload_mod(mod_name):
    if mod_name in avail_mods.keys():
        mod_env = avail_mods[mod_name]
        for key in mod_env.keys():
            if key in os.environ.keys():
                new_env = os.environ[key]
                new_env = new_env.replace(mod_env[key] + ":", "")
                new_env = new_env.replace(mod_env[key], "")
                print_cmd("export %s=%s;" % (key, new_env))
    else:
        print_cmd("echo \"Module %s not found\";" % mod_name)
        sys.exit(-1)


def build_full_mods(mod_name):
    full_mods = [mod_name]
    recur_flag = True
    recur_count = 0
    while recur_flag == True and recur_count <= 10000:
        recur_flag = False

        # Loop over all the existing mods in full_mods
        for mod_name in full_mods:
            # Check if the current mod has been defined in avail_mods
            if mod_name in avail_mods.keys():
                mod_env = avail_mods[mod_name]

                # Check if the dependencies are already in full_mods
                if "DEPENDENCY" in mod_env.keys():
                    for mod_dep in mod_env["DEPENDENCY"]:
                        if mod_dep not in full_mods:
                            full_mods.append(mod_dep)
                            recur_flag = True
                    recur_count += 1
            else:
                print_cmd("echo \"Module %s not found\";" % mod_name)
                sys.exit(-1)
    return list(reversed(full_mods))


def load_mod_recur(mod_name):
    full_mods = build_full_mods(mod_name)
    for mod in full_mods:
        print_cmd("pmod load %s;" % mod)
                

def unload_mod_recur(mod_name):
    full_mods = build_full_mods(mod_name)
    for mod in full_mods:
        print_cmd("pmod unload %s;" % mod)


def load_mod_all():
    full_mods = avail_mods.keys()
    for mod in full_mods:
        print_cmd("pmod load %s;" % mod)


def unload_mod_all():
    full_mods = avail_mods.keys()
    for mod in full_mods:
        print_cmd("pmod unload %s;" % mod)
    

def main():
    # Parse cli-parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--all", default=False, action="store_true")
    parser.add_argument("-r", "--recursive", default=False, action="store_true")
    parser.add_argument("operation", type=str, action="store")
    parser.add_argument("mod_name", type=str, action="store", nargs="?")
    args = parser.parse_args()
    
    # Perform the required operation
    if args.operation == "avail":
        avail_mod()
    elif args.operation in ("list", "ls", "status", "stat"):
        list_mod()
    elif args.operation in ("load", "add"):
        if args.all == True:
            load_mod_all()
        elif args.recursive == True:
            load_mod_recur(args.mod_name)
        else:
            load_mod(args.mod_name)
    elif args.operation in ("unload", "remove", "rm"):
        if args.all == True:
            unload_mod_all()
        elif args.recursive == True:
            unload_mod_recur(args.mod_name)
        else:
            unload_mod(args.mod_name)
    else:
        print_cmd("echo \"Unkown operation %s\";" % args.operation)
        sys.exit(-1)
 

if __name__ == "__main__":
    main()
