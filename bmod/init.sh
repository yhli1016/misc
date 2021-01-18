# Set up environment variables used by bmod
export BMOD_ROOT=$HOME/soft/bmod
export BMOD_MOD=$BMOD_ROOT/modules
export BMOD_LOADED_MODS=""

# Function for manipulating environment variables
reset_env () {
    local cmd=$1
    local envname=$2
    local pattern=$3
    if [[ "$cmd" == "add" ]]; then
        export $envname=$pattern
    elif [[ "$cmd" == "rm" ]]; then
        unset $envname
    else
        echo "ERROR: Illegal command '$cmd'"
    fi
}

set_env () {
    local cmd=$1
    local envname=$2
    local pattern=$3
    local envval=$(eval echo '$'"$envname")
    local flag=""
    local newval=""
    if [[ "$cmd" == "add" ]]; then
        if [[ -z "$envval" ]]; then
            export $envname=$pattern
        else
            flag=$(echo $envval | awk -F ':' \
                   '{for(i=1;i<=NF;i++) if($i=="'$pattern'") {print $i; break}}')
            test -z "$flag" && export $envname=$pattern:$envval
        fi
    elif [[ "$cmd" == "rm" ]]; then
        if [[ ! -z "$envval" ]]; then
            newval=$(echo $envval | awk -F ':' \
                     '{for(i=1;i<=NF;i++) if($i!="'$pattern'") printf "%s ", $i}' | \
                     awk '{for(i=1;i<NF;i++) printf "%s:", $i; print $NF}')
            test "$newval" != "$envval" && export $envname=$newval
        fi
    else
        echo "ERROR: Illegal command '$cmd'"
    fi
}

# Wrapper for set_env with presets of environment variables
set_mod () {
    local cmd=$1
    local presets=$(echo $2 | sed 's/+/\n/g')
    local pattern=$3
    local path=""
    for preset in $presets; do
        case $preset in
        bin)
            set_env $cmd "PATH" $pattern
            ;;
        lib)
            set_env $cmd "LIBRARY_PATH" $pattern
            set_env $cmd "LD_RUN_PATH" $pattern
            set_env $cmd "LD_LIBRARY_PATH" $pattern
            ;;
        config)
            set_env $cmd "PKG_CONFIG_PATH" $pattern
            ;;
        inc)
            set_env $cmd "C_INCLUDE_PATH" $pattern
            set_env $cmd "CPLUS_INCLUDE_PATH" $pattern
            ;;
        py)
            set_env $cmd "PYTHONPATH" $pattern
            ;;
        pkg)
            test -d $pattern/bin && set_env $cmd "PATH" $pattern/bin
            for path in lib lib64; do
                if [[ -d $pattern/$path ]]; then
                    set_env $cmd "LIBRARY_PATH" $pattern/$path
                    set_env $cmd "LD_RUN_PATH" $pattern/$path
                    set_env $cmd "LD_LIBRARY_PATH" $pattern/$path
                    test -d $pattern/$path/pkgconfig && \
                    set_env $cmd "PKG_CONFIG_PATH" $pattern/$path/pkgconfig
                fi
            done
            if [[ -d $pattern/include ]]; then
                set_env $cmd "C_INCLUDE_PATH" $pattern/include
                set_env $cmd "CPLUS_INCLUDE_PATH" $pattern/include
            fi
            ;;
        *)
            echo "ERROR: Illegal preset '$preset'"
            ;;
        esac
    done
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
                if [[ -f $path/$modname.sh ]]; then
                    script=$modname.sh
                    break
                else
                    script=$(ls $path | egrep "^$modname[-/]+[0-9\.]+\.sh$" | tail -1)
                    test $script && break
                fi
            done
            # Source the script file and update BMOD_LOADED_MODS
            if [[ -f $path/$script ]]; then
                source $path/$script $cmd
                if [[ $? == 0 ]]; then
                    set_env $cmd "BMOD_LOADED_MODS" $script
                fi
            else
                echo "ERROR: Module '$modname' not found"
                return -1
            fi
        done
    elif [[ "$cmd" == "ls" ]]; then
        echo "Currently loaded modules:"
        echo $BMOD_LOADED_MODS | sed -e 's/:/\n/g' -e 's/.sh//g' | tac | \
            awk '{if(NF>0) printf "%4i) %s\n", NR, $1}'
    elif [[ "$cmd" == "av" ]]; then
        for path in $(echo $BMOD_MOD | sed 's/:/\n/g'); do
            echo "---- $path ----"
            ls $path | sed 's/.sh//g' | sort | awk '{if(NF>0) printf "%4i) %s\n", NR, $1}'
        done
    elif [[ "$cmd" == "cl" ]]; then
        for path in $(echo $BMOD_MOD | sed 's/:/\n/g'); do
            for script in $(echo $BMOD_LOADED_MODS | sed 's/:/\n/g'); do
                if [[ -f $path/$script ]]; then
                    source $path/$script rm
                    set_env rm "BMOD_LOADED_MODS" $script
                fi
            done
        done
    elif [[ "$cmd" == "pg" ]]; then
        for path in $(echo $BMOD_MOD | sed 's/:/\n/g'); do
            for script in $(ls $path); do
                source $path/$script rm
                set_env rm "BMOD_LOADED_MODS" $script
            done
        done
    else
        echo "ERROR: Illegal command '$cmd'"
    fi
}

complete -W "add rm ls av cl pg" bmod
