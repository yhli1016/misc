#! /bin/bash

runtest()
{
	cat > sigma.inp << EOF
screened_coulomb_cutoff 16.0
bare_coulomb_cutoff     70.0

number_bands $n

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

for n in $(seq 200 200 2000); do
	runtest
done
