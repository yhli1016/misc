#! /bin/bash

# Back up configuration
mkdir -p bak/config
itemlist=".bash* .dircolors .inputrc .octaverc .vimrc \
.ssh"
for item in $itemlist; do
	cp -r ~/$item bak/config
done

# Back up software
mkdir -p bak/soft
itemlist="bin octave templates"
for item in $itemlist; do
    cp -r ~/soft/$item bak/soft
done

# Archiving
tar -cjf bak.tar.bz2 bak && rm -rf bak
