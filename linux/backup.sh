#! /bin/bash

# Back up configuration
itemlist="\
.alias \
.bash_logout \
.bashrc \
.dircolors \
.gitconfig \
.inputrc \
.profile \
.ssh \
.vim \
.vimrc"


mkdir rc_bak
for item in $itemlist; do
	cp -r ~/$item rc_bak
done

# Archiving
tar -cjf rc_bak_$(date +%F).tar.bz2 rc_bak && rm -rf rc_bak
