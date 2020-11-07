# Set up environment variables used by bmod
export BMOD_ROOT=$HOME/soft/bmod.py
export BMOD_MOD=$BMOD_ROOT/modules
export BMOD_LOADED_MODS=""
if [[ ! "$PYTHONPATH" =~ "$BMOD_ROOT/base" ]]; then
    export PYTHONPATH=$BMOD_ROOT/base:$PYTHONPATH
fi

# Functions for manipulating environment variables and modules from CLI
reset_env ()
{
    local cmd=$1
    local envname=$2
    local pattern=$3
    eval $(python -m sandbox reset $cmd $envname $pattern)
}

set_env () {
    local cmd=$1
    local envname=$2
    local pattern=$3
    eval $(python -m sandbox set_env $cmd $envname $pattern)
}

set_mod () {
    local cmd=$1
    local presets=$2
    local dest=$3
    eval $(python -m sandbox set_mod $cmd $presets $dest)
}

# Top-level 'bmod' command
bmod () {
    local cmd=$1
    local modname=""
    local path=""
    local script=""
    if [[ "$cmd" == "add" || "$cmd" == "rm" ]]; then
        shift 
        for modname in $*; do
            # Locate the script file
            for path in $(echo $BMOD_MOD | sed 's/:/\n/g'); do
                if [[ -f $path/$modname.py ]]; then
                    script=$modname.py
                    break
                else
                    script=$(ls $path | egrep "^$modname[-/]+[0-9\.]+\.py$" | tail -1)
                    test $script && break
                fi
            done
            # Source the script file
            if [[ -f $path/$script ]]; then
                eval $(python $path/$script $cmd)
            else
                echo "ERROR: Module '$modname' not found"
                return -1
            fi
        done
    elif [[ "$cmd" == "ls" ]]; then
        echo "Currently loaded modules:"
        echo $BMOD_LOADED_MODS | sed 's/:/\n/g'| tac | \
            awk '{if(NF>0) printf "%4i) %s\n", NR, $1}'
    elif [[ "$cmd" == "av" ]]; then
        for path in $(echo $BMOD_MOD | sed 's/:/\n/g'); do
            echo "---- $path ----"
            ls $path | sed 's/.py//g' | sort | awk '{if(NF>0) printf "%4i) %s\n", NR, $1}'
        done
    elif [[ "$cmd" == "cl" ]]; then
        for path in $(echo $BMOD_MOD | sed 's/:/\n/g'); do
            for modname in $(echo $BMOD_LOADED_MODS | sed 's/:/\n/g'); do
                test -f $path/$modname.py && eval $(python $path/$modname.py rm)
            done
        done
    elif [[ "$cmd" == "pg" ]]; then
        for path in $(echo $BMOD_MOD | sed 's/:/\n/g'); do
            for script in $(ls $path); do
                test -f $path/$script && eval $(python $path/$script rm)
            done
        done
    else
        echo "ERROR: Illegal command '$cmd'"
    fi
}

complete -W "add rm ls av cl pg" bmod
