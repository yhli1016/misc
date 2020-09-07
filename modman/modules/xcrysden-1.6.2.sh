if [[ "$1" == "add" ]]; then
    modadd bin+lib $HOME/soft/xcrysden-1.6.2-bin-shared
elif [[ "$1" == "rm" ]]; then
    modrm bin+lib $HOME/soft/xcrysden-1.6.2-bin-shared
else
    echo "Illegal command: $cmd"
fi