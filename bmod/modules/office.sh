set_mod $1 bin+lib $HOME/soft/doublecmd
set_mod $1 bin $HOME/soft/Typora-linux-x64
set_mod $1 bin $HOME/soft/TeXmacs/bin

if [[ "$1" == "add" ]]; then
    export TEXMACS_PATH=$HOME/soft/TeXmacs
else
    unset TEXMACS_PATH
fi
