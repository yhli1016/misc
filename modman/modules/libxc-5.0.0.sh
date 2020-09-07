if [[ "$1" == "add" ]]; then
    modadd pkg $HOME/soft/libxc-5.0.0
elif [[ "$1" == "rm" ]]; then
    modrm pkg $HOME/soft/libxc-5.0.0
else
    echo "Illegal command: $cmd"
fi
