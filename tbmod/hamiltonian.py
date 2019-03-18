"""Class representing the Hamiltonian."""

import numpy as np
from core import set_ham


class TBModel(object):
    """A simple yet efficient tight-binding model."""

    def __init__(self, num_orbital):
        self.Rn = None
        self.ij = None
        self.tij = None
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
        Rn_list, ij_list, tij_list = [], [], []
        for line in hr_content:
            line_split = line.split()
            Rn = [float(line_split[0]), float(line_split[1]),
                  float(line_split[2])]
            ij = [int(line_split[3]) - 1, int(line_split[4]) - 1]
            tij = float(line_split[5]) + 1j * float(line_split[6])
            if np.abs(tij) >= threshold:
                Rn_list.append(Rn)
                ij_list.append(ij)
                tij_list.append(tij)
        self.Rn = np.array(Rn_list, dtype="float")
        self.ij = np.array(ij_list, dtype="int32")
        self.tij = np.array(tij_list, dtype="complex")
        orbital_indices = set([ij[0] for ij in self.ij])
        assert len(orbital_indices) == self.hamiltonian.shape[0]

    def eval_energies(self, kpoints, orbital_list=None):
        """
        Get the energies for each kpoint.

        :param kpoints: N*3 array
        :param orbital_list: indices of orbitals on which the projection will
                be evaluated, counted from 1
        :return: energies, Nk * N_orbital array
        :return: projection, Nk * N_orbital array
        """
        num_k_point = len(kpoints)
        energies = np.zeros((num_k_point, self.hamiltonian.shape[0]),
                            dtype="float")
        projection = np.zeros((num_k_point, self.hamiltonian.shape[0]),
                              dtype="float")
        if orbital_list is not None:
            mask = np.array([True if i+1 in orbital_list else False
                             for i in range(self.hamiltonian.shape[0])])

        for ik, kpoint in enumerate(kpoints):
            set_ham(self.hamiltonian, self.ij, self.Rn, self.tij, kpoint)
            if orbital_list is not None:
                eigval, eigvec = np.linalg.eig(self.hamiltonian)
                energies[ik] = eigval.real
                for ib in range(self.hamiltonian.shape[0]):
                    projection[ik, ib] += np.sum(np.abs(eigvec[mask, ib])**2)
                    projection[ik, ib] /= np.sum(np.abs(eigvec[:, ib])**2)
            else:
                energies[ik] = np.sort(np.linalg.eigvals(self.hamiltonian)).real

        return energies, projection
