#! /usr/bin/env python

# input file for hopterm.py

# lattice vectors
lat = [
[ 2.027733378,   3.512137236,   0.000000000],
[-2.027733378,   3.512137236,   0.000000000],
[ 0.000000000,   0.000000000,  20.500000000]]

# fractional coordinates of orbitals
orb = [
[0.333334671,   0.333334671,   0.431781565],
[0.333334671,   0.333334671,   0.431781565],
[0.333334671,   0.333334671,   0.431781565],
[0.333334671,   0.333334671,   0.431781565],
[0.333334812,   0.333334812,   0.568224866],
[0.333334812,   0.333334812,   0.568224866],
[0.333334812,   0.333334812,   0.568224866],
[0.333334812,   0.333334812,   0.568224866],
[0.666665091,   0.666665091,   0.369744698],
[0.666665091,   0.666665091,   0.369744698],
[0.666665091,   0.666665091,   0.369744698],
[0.666665091,   0.666665091,   0.369744698],
[0.666665427,   0.666665427,   0.630258871],
[0.666665427,   0.666665427,   0.630258871],
[0.666665427,   0.666665427,   0.630258871],
[0.666665427,   0.666665427,   0.630258871]]

# types of orbitals
orbtype = [
"s", "px", "py", "pz",
"s", "px", "py", "pz",
"s", "px", "py", "pz",
"s", "px", "py", "pz"]

# types of atoms
atomtype = [
"In1", "In1", "In1", "In1",
"In2", "In2", "In2", "In2",
"Se1", "Se1", "Se1", "Se1",
"Se2", "Se2", "Se2", "Se2"]

# lower and upper bond for Rn and cutoff distance
namin = -10; namax = 10;
nbmin = -10; nbmax = 10;
ncmin =  0; ncmax = 0;
dmax = 10.0;

# flags controls how the hopping terms are reduced
reduce_onsite = True
reduce_conjugate = True

# flags controls the output
write_lat = True
write_orb = True
output_flavor = "matlab"
