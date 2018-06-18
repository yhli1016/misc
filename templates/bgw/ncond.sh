#! /bin/bash

runtest()
{
	cat > absorption.inp << EOF
number_val_bands_coarse  9
number_cond_bands_coarse 11

number_val_bands_fine  5
number_cond_bands_fine $i

use_symmetries_coarse_grid
use_symmetries_fine_grid
use_symmetries_shifted_grid

eqp_co_corrections

diagonalization

spin_singlet

screening_semiconductor

cell_slab_truncation

use_momentum
polarization 1.0 0.0 0.0

gaussian_broadening
energy_resolution 0.05

write_eigenvectors -1
EOF

	mpirun -np 20 absorption.cplx.x &> absorption.out
	wait

	mv absorption.inp $i.inp
	mv absorption.out $i.out
	mv absorption_eh.dat $i.abs.dat
	mv eigenvalues.dat $i.eig.dat

}

for i in $(seq 2 2 10); do
	runtest
done
