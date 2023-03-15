#! /bin/bash

# Back up configuration
itemlist="\
.alias \
.bash_profile \
.bash_login \
.profile \
.bashrc \
.bash_logout \
.dircolors \
.gitconfig \
.inputrc \
.ssh \
.vim \
.vimrc \
soft/bmod"


mkdir rc_bak
for item in $itemlist; do
	cp -r ~/$item rc_bak
done

# Archiving
tar -cjf rc_bak_$(date +%F).tar.bz2 rc_bak && rm -rf rc_bak
