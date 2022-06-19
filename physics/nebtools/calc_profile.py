#! /usr/bin/env python
"""Calculate energy profile of reaction paths."""

import numpy as np

# Energies of shared reactants and products
CO2 = -18.40232529  # Energy of CO2 molecule
CO = -12.06750907   # Energy of CO molecule
H2 = -7.16771221    # Energy of H2 molecule
H2O = -12.80802156  # Energy of H2O molecule


class Path:
    """
    Class for evaluating the energy profile of a single reaction path.

    Attributes
    ----------
    label: list of strings
        labels for reactants and products
    energy: list of floats
        energies of reactants and products
    unit: string
        unit of energy
    """
    def __init__(self, unit="ev") -> None:
        self.label = []
        self.energy = []
        if unit not in ("kjm", "ev"):
            raise ValueError(f"Illegal unit: {unit}")
        self.unit = unit

    def scale_energy(self, unit="ev"):
        """Get scaled energy in given unit."""
        ev2kjm = 96.4916
        if unit not in ("kjm", "ev"):
            raise ValueError(f"Illegal unit: {unit}")
        if unit == "kjm":
            if self.unit == "kjm":
                scale_factor = 1.0
            else:
                scale_factor = ev2kjm
        else:
            if self.unit == "ev":
                scale_factor = 1.0
            else:
                scale_factor = 1.0 / ev2kjm
        energy = np.array(self.energy) * scale_factor
        return energy

    def add_eng(self, label, energy):
        """
        Add an energy level in the reaction path.

        :param string label: label for the state
        :param float energy: energy of the state
        :return: None
        """
        if energy is not None:  # DO NOT DELETE THIS LINE!
            self.label.append(label)
            self.energy.append(energy)

    def eval_eng(self, unit="eV"):
        """
        Print energy levels and differences of the reaction path.

        :param string unit: unit of energies for output
        :return: None
        """
        energy = self.scale_energy(unit=unit)
        for i, label in enumerate(self.label):
            eng = energy[i]
            eng_align = eng - energy[0]
            eng_delta = eng - energy[i-1] if i > 0 else 0
            if i > 0:
                print("%16s : %8.2f%8.2f" % (label, eng_align, eng_delta))
            else:
                print("%16s : %8.2f%8s" % (label, eng_align, "diff"))


class MultiPath:
    """
    Base class for evaluating the energy profile of multiple reaction paths.

    Attributes
    ----------
    paths: List[Path]
        list of single reaction paths
    """
    def __init__(self) -> None:
        self.paths = []

    def gen_paths(self):
        """To be implemented in derived classes."""
        pass

    def eval_eng(self, unit="ev"):
        """
        Print energy levels and differences of the reaction paths.

        :param string unit: unit of energies for output
        :return: None
        """
        self.gen_paths()
        for path in self.paths:
            path.eval_eng(unit=unit)
            print()


class RWGS(MultiPath):
    """Class for evaluating the energy profile of RWGS reaction."""
    def __init__(self) -> None:
        super().__init__()

        # Energy of substrate
        self.sub = None
    
        # Energies of CO formation
        self.co2_fe_linear = None  # CO2-Fe* in linear configuration
        self.co2_l2b_ts = None     # Transition state from linear to bending
        self.co2_fe_bend = None    # CO2-Fe* in bent configuration
        self.co_form_ts = None     # Transition state of co formation
        self.co_fe_o_sub = None    # Co-Fe* and O-*
        self.o_sub = None          # O-* after release of CO gas

        # Energies of H2O formation
        self.h2_fe = None         # H2-Fe*
        self.h2o_form1_ts = None  # Transition state of migration of 1st H atom
        self.h_fe_h_o = None      # H-Fe* and H-O*
        self.h2o_form2_ts = None  # Transition state of migration of 2nd H atom
        self.h2o_sub = None       # H2O-*

    def gen_paths(self):
        co_path = Path()
        co_path.add_eng('co2 + sub', CO2+self.sub)
        co_path.add_eng('co2_fe_linear', self.co2_fe_linear)
        co_path.add_eng('co2_l2b_ts', self.co2_l2b_ts)
        co_path.add_eng('co2_fe_bend', self.co2_fe_bend)
        co_path.add_eng('co_form_ts', self.co_form_ts)
        co_path.add_eng('co_fe_o_sub', self.co_fe_o_sub)
        co_path.add_eng('co + o_sub', CO+self.o_sub)
        self.paths = [co_path]

        h2o_path = Path()
        h2o_path.add_eng('h2 + o_sub', H2+self.o_sub)
        h2o_path.add_eng('h2_fe', self.h2_fe)
        h2o_path.add_eng('h2o_form1_ts', self.h2o_form1_ts)
        h2o_path.add_eng('h_fe_h_o', self.h_fe_h_o)
        h2o_path.add_eng('h2o_form2_ts', self.h2o_form2_ts)
        h2o_path.add_eng('h2o_sub', self.h2o_sub)
        h2o_path.add_eng('h2o + sub', H2O+self.sub)
        self.paths.append(h2o_path)


def main():
    # 0 coverage
    print("Reaction Energies of RWGS on Fe@Mo2C with 0 coverage")
    path = RWGS()
    path.sub = -198.12728341
    path.co2_fe_bend = -217.79137763
    path.co_form_ts = -217.401000
    path.co_fe_o_sub = -218.62729779
    path.o_sub = -204.99528597
    path.h2_fe = -212.40798162
    path.h2o_form1_ts = -212.002300
    path.h_fe_h_o = -212.31001109
    path.h2o_form2_ts = -211.133700
    path.h2o_sub = -211.62499841
    path.eval_eng()

    # 0.33 coverage
    print("Reaction Energies of RWGS on Fe@Mo2C with 0.33 coverage")
    path = RWGS()
    path.sub = -218.87184700
    path.co2_fe_bend = -238.67602257
    path.co_form_ts = -237.759400
    path.co_fe_o_sub = -239.33374512
    path.o_sub = -225.50533791
    path.h2_fe = -233.25968018
    path.h2o_form1_ts = -232.698800
    path.h_fe_h_o = -233.61507818
    path.h2o_sub = -232.48088446
    path.h2o_form2_ts = -231.723
    path.h2o_sub = -232.48088446
    path.eval_eng()

    # 0.78
    print("Reaction Energies of RWGS on Fe@Mo2C with 0.78 coverage")
    path = RWGS()
    path.sub = -245.94153180
    path.co2_fe_bend = -265.30342210
    path.co_form_ts = -264.997200
    path.co_fe_o_sub = -266.0301423
    path.o_sub = -252.68471571
    path.h2_fe = -260.11341903
    path.h2o_form1_ts = -259.498800
    path.h_fe_h_o = -259.86970027
    path.h2o_form2_ts = -259.021500
    path.h2o_sub = -259.72247630
    path.eval_eng()


if __name__ == "__main__":
    main()
