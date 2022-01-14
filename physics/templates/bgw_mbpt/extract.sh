#! /bin/bash

for i in $(seq 200 200 2000); do
    vbm=$(awk '$1==9 {print $10}' $i.log)
    cbm=$(awk '$1==10 {print $10}' $i.log)
    echo $i $vbm $cbm | awk '{print $1, $2, $3, $3 - $2}'
done
