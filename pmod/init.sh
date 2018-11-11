PMODPATH=$HOME/soft/pmod
export PATH=$PMODPATH:$PATH
if [ ! -x "$PMODPATH/modcmd.py" ]; then
    chmod +x $PMODPATH/modcmd.py
fi

function pmod ()
{
    cmds=`modcmd.py $*`
    eval $cmds
}
alias module=pmod
