#! /usr/bin/env python
import re
import math
from typing import Tuple, List, Dict

import numpy as np
import matplotlib.pyplot as plt


def gaussian(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    """
    Gaussian type broadening function.

    :param x: incoming x
    :param mu: center of the Gaussian function
    :param sigma: half-width of the Gaussian function
    :return: normalized Gaussian function value at each x
    """
    part_a = 1.0 / (sigma * math.sqrt(2 * math.pi))
    part_b = np.exp(-(x - mu)**2 / (2 * sigma**2))
    return part_a * part_b


def get_fermi(outcar_name: str = "OUTCAR") -> float:
    """
    Get Fermi energy from OUTCAR.

    :param outcar_name: file name of OUTCAR
    :return: fermi energy
    """
    with open(outcar_name, "r") as outcar:
        content = outcar.readlines()
    fermi = 0.0
    for line in content:
        result = re.search(r"\s*Fermi energy:\s*(-*\d*\.*\d*)", line)
        if result is not None:
            fermi = float(result.group(1))
    return fermi


class StateProjection:
    """
    Container for holding the projection of single state on atomic orbitals.

    Attributes
    ----------
    energy: float
        energy of the state
    projection: List[Dict[str, float]]
        projection of the state of atomic orbitals, each item corresponds to an
        ion, with keys being orbital symbols and values being projections
    """
    def __init__(self, energy: float,
                 projection: List[Dict[str, float]]) -> None:
        self.energy = energy
        self.projection = projection


class Procar:
    """
    Container for holding the projections of all states in PROCAR.

    Attributes
    ----------
    _data: List[StateProjection]
        projections of all states
    """
    def __init__(self, procar_name: str = None) -> None:
        """
        :param procar_name: name of PROCAR
        """
        self._data = []
        if procar_name is not None:
            self.parse(procar_name)

    def parse(self, procar_name: str = "PROCAR") -> None:
        """
        Load and parse PROCAR.

        :param procar_name: name of PROCAR
        :return: None
        """
        # Clear existing data
        self._data = []

        # Load data and trim blank lines
        with open(procar_name, "r") as procar_file:
            procar = procar_file.readlines()
        procar = [_ for _ in procar if len(_.split()) != 0]

        # Get number of ions
        result = re.search(r"# of ions:\s*(\d*)", procar[1])
        num_ion = int(result.group(1))

        # Determine header lines
        nl_start = [_ for _, line in enumerate(procar)
                    if re.search(r"^band", line) is not None]

        # Extract data
        for nl0 in nl_start:
            # Get energy
            result = re.search(r"# energy\s*(-*\d*\.*\d*)", procar[nl0])
            energy = float(result.group(1))

            # Assemble projection
            projection = []
            symbols = procar[nl0+1].split()
            for line in procar[nl0+2:nl0+2+num_ion]:
                proj_ion = dict()
                data = line.split()
                for i, s in enumerate(symbols):
                    if i != 0:
                        proj_ion[s] = float(data[i])
                projection.append(proj_ion)
            self._data.append(StateProjection(energy, projection))

    def sum_proj(self, atom_ids: List[int],
                 symbols: List[str]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Sum over specified atom indices and orbital symbols.

        :param atom_ids: atom indices
        :param symbols: orbital symbols
        :return: (energy, projection)
        """
        energy_array, proj_array = [], []
        for state in self._data:
            proj_sum = 0.0
            for i in atom_ids:
                for s in symbols:
                    proj_sum += state.projection[i][s]
            energy_array.append(state.energy)
            proj_array.append(proj_sum)
        return np.array(energy_array), np.array(proj_array)

    def eval_pdos(self, atom_ids: List[int],
                  symbols: List[str],
                  e_min: float = -2.0,
                  e_max: float = 2.0,
                  e_step: float = 0.01,
                  sigma: float  = 0.1) -> Tuple[np.ndarray, np.ndarray]:
        """
        Evaluate PDOS.

        :param atom_ids: atom indices
        :param symbols: orbital symbols
        :param e_min: lower bound of energies
        :param e_max: upper bound of energies
        :param e_step: resolution of energies
        :param sigma: broadening parameter
        :return: (energy, pdos)
        """
        num_grid = int((e_max - e_min) / e_step)
        energy_fi = np.linspace(e_min, e_max, num_grid + 1)
        pdos_fi = np.zeros_like(energy_fi)
        energy_co, pdos_co = self.sum_proj(atom_ids, symbols)
        for i, eng in enumerate(energy_co):
            pdos_fi += gaussian(energy_fi, eng, sigma) * pdos_co[i]
        pdos_fi /= energy_co.size
        return energy_fi, pdos_fi


def main():
    atom_ids = [1]
    symbols = [["x2-y2"], ["dxz", "dyz"]]
    e_step = 0.001
    sigma = 0.05
    strain = ["0.98", "0.99", "1.00", "1.01", "1.02"]
    direct = "r_c"

    for r in ["rot"]:
        for s in strain:
            fermi = get_fermi(f"{direct}/{s}/OUTCAR.{r}")
            e_min = fermi - 1.5
            e_max = fermi
            procar = Procar(f"{direct}/{s}/PROCAR.{r}")
            # energy, pdos0 = procar.eval_pdos(atom_ids, symbols[0],
            #                                  e_min, e_max, e_step, sigma)
            energy, pdos1 = procar.eval_pdos(atom_ids, symbols[1],
                                             e_min, e_max, e_step, sigma)
            overlap = np.sum(pdos1) * e_step
            print(r, s, f"{overlap:.3e}")

    for r in ["rot"]:
        for s in strain:
            fermi = get_fermi(f"{direct}/{s}/OUTCAR.{r}")
            e_min = fermi - 1.5
            e_max = fermi + 1.5
            procar = Procar(f"{direct}/{s}/PROCAR.{r}")
            energy, pdos0 = procar.eval_pdos(atom_ids, symbols[0],
                                             e_min, e_max, e_step, sigma)
            energy, pdos1 = procar.eval_pdos(atom_ids, symbols[1],
                                             e_min, e_max, e_step, sigma)
            energy -= fermi
            plt.plot(energy, pdos0, label="d//")
            plt.plot(energy, pdos1, label="$\pi*$")
            plt.axvline(color="k", linestyle="--")
            plt.xlim(np.min(energy), np.max(energy))
            plt.title(f"{r}, {s}")
            plt.legend()
            plt.show()



if __name__ == "__main__":
    main()
