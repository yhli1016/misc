#!/bin/sh
# compute dependencies for the PWscf directory tree

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# run from directory where this script is
echo "The script for adding new dependencies is running"
echo "Usage: $0 DIR DEPENDENCY_DIRS"
echo "dependencies in DEPENDENCY_DIRS will be searched for DIR"
shift
DEPENDS=$*

# generate dependencies file (only for directories that are present)
moduledep.sh $DEPENDS > make.depend
includedep.sh $DEPENDS >> make.depend

# handle special cases: modules for C-fortran binding, hdf5, MPI
sed '/@mpi@/d' make.depend > make.depend.tmp
mv make.depend.tmp make.depend

# remove cyclic dependencies
sed 's/://g' make.depend | \
awk '!match($2, $1) {printf "%s: %s\n", $1, $2}' > make.depend.tmp
mv make.depend.tmp make.depend

# check for missing dependencies 
if grep @ make.depend
then
   notfound=1
   echo WARNING: dependencies not found in directory $DIR
else
    echo directory $DIR : ok
fi

if test "$notfound" = ""
then
    echo all dependencies updated successfully
fi
