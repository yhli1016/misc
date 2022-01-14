#! /bin/bash

runtest()
{
	cat > sigma.inp << EOF
screened_coulomb_cutoff 12
bare_coulomb_cutoff     $n

number_bands 1280

screening_semiconductor

cell_slab_truncation

band_index_min 9
band_index_max 10

begin kpoints
  0.000000000  0.000000000  0.000000000   1.0
end
EOF
	mpirun -np 20 sigma.cplx.x &> sigma.out
	wait
	mv sigma.inp $n.inp
	mv sigma.out $n.out
	mv sigma_hp.log $n.log
}

for n in $(seq 15 5 70); do
	runtest
done
