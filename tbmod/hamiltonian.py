"""Functions for setting up the Hamiltonian and diagonalization."""

import numpy as np


class TBModel(object):
    def __init__(self, num_orbital):
        self.hr_terms = []
        self.hamiltonian = np.zeros((num_orbital, num_orbital), dtype="complex")

    def read_hr(self, filename, threshold=1.0e-4):
        """
        Read in the the matrix elements of Hamiltonian in LCAO basis
        (seedname_hr.dat) produced by Wannier90.x.

        :param filename: filename of seedname_hr.dat
        :param threshold: float, threshold for the matrix elements
        :return: None.
        """
        with open(filename, "r") as hr_file:
            hr_content = hr_file.readlines()
        for line in hr_content:
            line_split = line.split()
            Rn = np.array([float(line_split[0]), float(line_split[1]),
                           float(line_split[2])])
            ij = (int(line_split[3])-1, int(line_split[4])-1)
            tij = float(line_split[5]) + 1j * float(line_split[6])
            if np.abs(tij) >= threshold:
                self.hr_terms.append([Rn, ij, tij])

    def set_hamiltonian(self, kpoint):
        """
        Set up the Hamiltonian in Bloch basis set for given kpoint.

        :param kpoint: 1*3 array
        :return: None.
        """
        self.hamiltonian *= 0.0
        for hr in self.hr_terms:
            Rn = hr[0]
            ij = hr[1]
            tij = hr[2]
            hij = tij * np.exp(1j * 2 * np.pi * np.dot(kpoint, Rn))
            self.hamiltonian[ij] += hij

    def eval_energies(self, kpoints, orbital_list=None):
        """
        Get the energies for each kpoint.

        :param kpoints: N*3 array
        :param orbital_list: indices of orbitals on which the projection will
                be evaluated
        :return: energies, Nk * N_orbital array
        :return: projection, Nk * N_orbital array
        """
        orbital_indices = set([hr[1][0] for hr in self.hr_terms])
        assert len(orbital_indices) == self.hamiltonian.shape[0]
        num_k_point = len(kpoints)
        energies = np.zeros((num_k_point, self.hamiltonian.shape[0]),
                            dtype="float")
        if orbital_list is not None:
            projection = np.zeros((num_k_point, self.hamiltonian.shape[0]),
                                  dtype="float")
            mask = np.array([True if i in orbital_list else False
                             for i in range(self.hamiltonian.shape[0])])

        for ik, kpoint in enumerate(kpoints):
            self.set_hamiltonian(kpoint)
            if orbital_list is not None:
                eigval, eigvec = np.linalg.eig(self.hamiltonian)
                energies[ik] = eigval.real
                for ib in range(self.hamiltonian.shape[0]):
                    projection[ik, ib] += np.sum(np.abs(eigvec[mask, ib])**2)
                    projection[ik, ib] /= np.sum(np.abs(eigvec[:, ib])**2)
            else:
                energies[ik] = np.sort(np.linalg.eigvals(self.hamiltonian)).real

        return energies, projection
