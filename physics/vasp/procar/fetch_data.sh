#! /bin/bash

direct="r_c"
for i in 0.98 0.99 1.00 1.01 1.02; do
    mkdir $i
    scp h3clogin:~/vo2/pdos/$direct/$i/{OUTCAR.*,PROCAR.*} $i
done
