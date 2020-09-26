# Set up environment variables
export BMOD_ROOT=$HOME/proj/misc/bmod
export BMOD_MOD=$BMOD_ROOT/modules
export BMOD_LOADED_MODS=""

# Function for setting environment variables
set_env () {
    local cmd=$1
    local envname=$2
    local pattern=$3
    local envval=$(eval echo '$'"$envname")
    local flag=""
    local newval=""
    if [[ "$cmd" == "add" ]]; then
        if [ -z "$envval" ]; then
            export $envname=$pattern
        else
            flag=$(echo $envval | awk -F ':' 'BEGIN{flag=0} \
                   {for(i=1;i<=NF;i++) if($i=="'$pattern'") flag=1} \
                   END{print flag}')
            if [[ "$flag" == "0" ]]; then
                export $envname=$pattern:$envval
            fi
        fi
    elif [[ "$cmd" == "rm" ]]; then
        if [[ ! -z "$envval" ]]; then
            newval=$(echo $envval | awk -F ':' '{for(i=1;i<=NF;i++) \
                     if($i!="'$pattern'") printf "%s ", $i}' | \
                     awk '{for(i=1;i<NF;i++) printf "%s:", $i; print $NF}')
            if [[ ! "$newval" == "$envval" ]]; then
                export $envname=$newval
            fi
        fi
    else
        echo "ERROR: Illegal command '$cmd'"
    fi
}

# Wrapper for set_env with presets of environment variables
set_mod () {
    local cmd=$1
    local presets=$(echo $2 | sed 's/+/ /g')
    local pattern=$3
    for preset in $presets; do
        case $preset in
        bin)
            set_env $cmd "PATH" $pattern
            ;;
        lib)
            set_env $cmd "LIBRARY_PATH" $pattern
            set_env $cmd "LD_LIBRARY_PATH" $pattern
            ;;
        inc)
            set_env $cmd "C_INCLUDE_PATH" $pattern
            set_env $cmd "CXX_INCLUDE_PATH" $pattern
            ;;
        config)
            set_env $cmd "PKG_CONFIG_PATH" $pattern
            ;;
        py)
            set_env $cmd "PYTHONPATH" $pattern
            ;;
        pkg)
            test -d $pattern/bin && set_mod ${cmd} bin $pattern/bin
            test -d $pattern/lib && set_mod ${cmd} lib $pattern/lib
            test -d $pattern/lib64 && set_mod ${cmd} lib $pattern/lib64
            test -d $pattern/include && set_mod ${cmd} inc $pattern/include
            test -d $pattern/lib/pkgconfig && set_mod ${cmd} config $pattern/lib/pkgconfig
            test -d $pattern/lib64/pkgconfig && set_mod ${cmd} config $pattern/lib64/pkgconfig
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
    local script=""
    if [[ "$cmd" == "add" || "$cmd" == "rm" ]]; then
        shift 
        for modname in $*; do
            # Get the full name of script file
            if [ -f $BMOD_MOD/$modname.sh ]; then
                script=$modname.sh
            else
                script=$(ls $BMOD_MOD | egrep "^$modname[-/]+[0-9\.]+\.sh" | tail -1)
            fi
            # Source the script file and update BMOD_LOADED_MODS
            if [ -f $BMOD_MOD/$script ]; then
                source $BMOD_MOD/$script $cmd
                if [[ $? == 0 ]]; then
                    modname=$(echo $script | awk -F '.sh' '{print $1}')
                    set_env $cmd "BMOD_LOADED_MODS" $modname
                fi
            else
                echo "ERROR: Module '$modname' not found"
                return -1
            fi
        done
     elif [[ "$cmd" == "ls" ]]; then
        echo "Loaded modules:"
        echo $BMOD_LOADED_MODS | \
            awk -F ':' '{for(i=1;i<=NF;i++) print $i}' | \
            sort | awk '{printf "%4i) %s\n", NR, $1}'
     elif [[ "$cmd" == "av" ]]; then
        echo "Available modules:"
        ls $BMOD_MOD | awk -F '.sh' '{print $1}' | sort | \
            awk '{printf "%4i) %s\n", NR, $1}'
     else
        echo "ERROR: Illegal command '$cmd'"
     fi
}

complete -W "add rm ls av" bmod
