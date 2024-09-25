#! /usr/bin/env python

from nebtools.profile import ReactionPath, PathCollection


# Energies of shared reactants and products
CO2 = -18.40232529  # Energy of CO2 molecule
CO = -12.06750907   # Energy of CO molecule
H2 = -7.16771221    # Energy of H2 molecule
H2O = -12.80802156  # Energy of H2O molecule


class RWGS(PathCollection):
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

    def build_paths(self):
        co_path = ReactionPath()
        co_path.add_state('co2 + sub', CO2 + self.sub)
        co_path.add_state('co2_fe_linear', self.co2_fe_linear)
        co_path.add_state('co2_l2b_ts', self.co2_l2b_ts)
        co_path.add_state('co2_fe_bend', self.co2_fe_bend)
        co_path.add_state('co_form_ts', self.co_form_ts)
        co_path.add_state('co_fe_o_sub', self.co_fe_o_sub)
        co_path.add_state('co + o_sub', CO + self.o_sub)
        self.paths = [co_path]

        h2o_path = ReactionPath()
        h2o_path.add_state('h2 + o_sub', H2 + self.o_sub)
        h2o_path.add_state('h2_fe', self.h2_fe)
        h2o_path.add_state('h2o_form1_ts', self.h2o_form1_ts)
        h2o_path.add_state('h_fe_h_o', self.h_fe_h_o)
        h2o_path.add_state('h2o_form2_ts', self.h2o_form2_ts)
        h2o_path.add_state('h2o_sub', self.h2o_sub)
        h2o_path.add_state('h2o + sub', H2O + self.sub)
        self.paths.append(h2o_path)


def main():
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
    path.eval_energy_differences()


if __name__ == "__main__":
    main()
