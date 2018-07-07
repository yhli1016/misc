#! /bin/bash

runtest()
{
	cat > epsilon.inp << EOF
epsilon_cutoff 12.0

number_bands $n

cell_slab_truncation

begin qpoints
  0.001000000  0.000000000  0.000000000   1.0   1
  0.000000000  0.166666667  0.000000000   1.0   0
  0.000000000  0.333333333  0.000000000   1.0   0
  0.000000000  0.500000000  0.000000000   1.0   0
  0.000000000  0.666666667  0.000000000   1.0   0
  0.000000000  0.833333333  0.000000000   1.0   0
  0.166666667  0.333333333  0.000000000   1.0   0
  0.166666667  0.500000000  0.000000000   1.0   0
  0.166666667  0.666666667  0.000000000   1.0   0
  0.333333333  0.666666667  0.000000000   1.0   0
end
EOF

	cat > sigma.inp << EOF
screened_coulomb_cutoff 12
bare_coulomb_cutoff     30

number_bands 1280

screening_semiconductor

cell_slab_truncation

band_index_min 9
band_index_max 10

begin kpoints
  0.000000000  0.000000000  0.000000000   1.0
end
EOF
	mpirun -np 20 epsilon.cplx.x &> epsilon.out
	wait
	mpirun -np 20 sigma.cplx.x &> sigma.out
	wait

	mv epsilon.inp $n.epsilon.inp
	mv epsilon.out $n.epsilon.out
	mv epsilon.log $n.epsilon.log
	mv sigma.inp $n.sigma.inp
	mv sigma.out $n.sigma.out
	mv sigma_hp.log $n.sigma.log
}

for n in $(seq 209 200 2009); do
	runtest
done
