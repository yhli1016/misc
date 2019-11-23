#! /bin/bash

envman () {
	local cmd=$1
	local envname=$2
    local envval=$(eval echo '$'"$envname")
	local token=$3
	case $cmd in
	add)
		if [ -z "$envval" ]; then
			export $envname=$token
		elif [[ ! "$envval" =~ "$token" ]]; then
			export $envname=$token:$envval
		fi
		;;
	rm)
		if [[ "$envval" =~ "$token" ]]; then
			export $envname=$(echo $envval | sed 's,'$token:',,g' | sed 's,'$token',,g')
		fi
		;;
	*)
		echo "Illegal command $cmd"
		;;
	esac
}

modman () {
	local cmd=$1
	local preset=$2
	local token=$3
	case $preset in
	bin)
		envman $cmd "PATH" $token 
		;;
	lib)
		envman $cmd "LIBRARY_PATH" $token 
		envman $cmd "LD_LIBRARY_PATH" $token 
		;;
	inc)
		envman $cmd "C_INCLUDE_PATH" $token 
		envman $cmd "CXX_INCLUDE_PATH" $token 
		;;
	py)
		envman $cmd "PYTHONPATH" $token 
		;;
	pkg)
		test -d $2/bin && modman ${cmd} bin $token/bin
		test -d $2/lib && modman ${cmd} lib $token/lib
		test -d $2/lib64 && modman ${cmd} lib $token/lib64
		test -d $2/include && modman ${cmd} inc $token/include
		;;
	*)
		echo "Illegal preset $preset"
		;;
	esac
}

# for programs
modman add bin $HOME/software/bin
modman add bin $HOME/software/pycharm-community-2019.2.3/bin
modman add lib $HOME/software/xcrysden-1.6.2-bin-shared
modman add bin $HOME/software/xcrysden-1.6.2-bin-shared
modman add pkg $HOME/software/elpa-2019.05.002
modman add lib $HOME/software/wannier90-3.0.0
modman add bin $HOME/software/wannier90-3.0.0
modman add bin $HOME/software/fleurMaXR3/build
#modman add bin $HOME/software/fleurMaXR3.1/build
modman add bin $HOME/software/elk-6.3.2/bin
