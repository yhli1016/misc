#! /bin/bash

# Back up configuration
itemlist=".bash_profile .bashrc .bash_logout \
.dircolors .inputrc .vimrc .ssh"

mkdir rc_bak
for item in $itemlist; do
	cp -r ~/$item rc_bak
done

# Archiving
tar -cjf rc_bak.tar.bz2 rc_bak && rm -rf rc_bak
