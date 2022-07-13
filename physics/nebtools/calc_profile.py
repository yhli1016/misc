#! /usr/bin/env python

from nebtools.profile import Path, MultiPath


# Energies of shared reactants and products
CO2 = -18.40232529  # Energy of CO2 molecule
CO = -12.06750907   # Energy of CO molecule
H2 = -7.16771221    # Energy of H2 molecule
H2O = -12.80802156  # Energy of H2O molecule


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
