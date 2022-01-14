#! /bin/bash

## read valence bands maxima of MF(LDA) and GW
#echo "Give VBM of DFT energies"
#read vbm_mf
#echo "Give VBM of GW energies"
#read vbm_qp

vbm_mf=$1
vbm_qp=$2

# remove the head lines of bandstructure.dat
awk 'NR>2 {print}' bandstructure.dat > effdata.dat

# determine band_index_min and band_index_max
band_index_min=$(awk '{print $2}' effdata.dat | sort -g | uniq | head -1)
band_index_max=$(awk '{print $2}' effdata.dat | sort -g | uniq | tail -1)

# extract band structure for every band within [ band_index_min, band_index_max ]
for i in $(seq $band_index_min $band_index_max); do
    awk '$2=='$i' {printf "%15.9f\n", $6 - '$vbm_mf'}' effdata.dat > mf_$i.dat
    awk '$2=='$i' {printf "%15.9f\n", $7 - '$vbm_qp'}' effdata.dat > qp_$i.dat
done

# generate kpath.dat
awk '$2=='$band_index_min' {printf "%12.5f%12.5f%12.5f\n", $3, $4, $5}' effdata.dat > kpoints.dat
kpath.py kpoints.dat kpath.dat

# assemble all data files
for flavor in mf qp; do
    if [ -s bnd_$flavor.dat ]; then
        rm bnd_$flavor.dat
    fi
    touch bnd_$flavor.dat
    for i in $(seq $band_index_min $band_index_max); do
        paste bnd_$flavor.dat "$flavor"_"$i".dat | expand -t 2 > temp.dat
        mv temp.dat bnd_$flavor.dat
    done
    paste kpath.dat bnd_$flavor.dat | expand -t 2 > temp.dat
    mv temp.dat bnd_$flavor.dat
done

# sort energies
sorteqp.py bnd_mf.dat bnd_mf.dat
sorteqp.py bnd_qp.dat bnd_qp.dat

# clean up
rm mf_*.dat qp_*.dat kpoints.dat kpath.dat effdata.dat
