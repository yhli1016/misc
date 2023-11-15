from typing import Tuple, Dict, Iterable, Union
from copy import deepcopy

import sympy as sp


c_type = Union[int, float, complex, sp.Basic]


class AtomicOrbital:
    """
    Class representing an atomic orbital composed of multiple |l,m,s> states.

    Attributes
    ----------
    _l_max: int
        maximum angular quantum number l
    _coeff: Dict[Tuple[int, int, int], c_type]
        keys: (l, m, s), values: coefficients on |l,m,s>, where s is the quantum
        number of sigma_z which should be either -1 or 1
    _allowed_qn: Set[Tuple[int, int, int]]
        set of allowed (l, m, s) pairs
    """
    def __init__(self, l_max: int = 2) -> None:
        """
        :param l_max: maximum angular quantum number l
        """
        self._l_max = l_max
        self._coeff = self._init_coeff()
        self._allowed_qn = set(self._coeff.keys())

    def __setitem__(self, key: Tuple[int, int, int],
                    value: c_type = 1) -> None:
        """
        Set the coefficient on state |l,m,s>.

        :param key: (l, m, s) of the state
        :param value: the new coefficient
        :return: None
        """
        if key not in self._allowed_qn:
            raise KeyError(f"Undefined key {key}")
        self._coeff[key] = value

    def __getitem__(self, key: Tuple[int, int, int]) -> c_type:
        """
        Get the coefficient on state |l,m,s>.

        :param key: (l, m, s) of the state
        :return: the coefficient
        """
        if key not in self._allowed_qn:
            raise KeyError(f"Undefined key {key}")
        return self._coeff[key]

    def _init_coeff(self) -> Dict[Tuple[int, int, int], c_type]:
        """
        Build initial coefficients.

        :return: dictionary with keys being (l, m, s) and values being zero
        """
        coefficients = dict()
        for l_i in range(self._l_max + 1):
            for m in range(-l_i, l_i+1):
                for s in (-1, 1):
                    coefficients[(l_i, m, s)] = 0
        return coefficients

    def l_plus(self) -> None:
        """
        Apply the l+ operator on the atomic orbital.

        As the formula is
            l+|l,m> = sqrt((l - m) * (l + m + 1)) * h_bar * |l,m+1>
        the actual coefficient should be multiplied by a factor of h_bar.

        :return: None
        """
        new_coefficients = self._init_coeff()
        for key, value in self._coeff.items():
            l, m, s = key
            key_new = (l, m+1, s)
            if key_new in self._allowed_qn:
                factor = sp.sqrt((l - m) * (l + m + 1))
                new_coefficients[key_new] = value * factor
        self._coeff = new_coefficients

    def l_minus(self) -> None:
        """
        Apply the l- operator on the atomic orbital.

        As the formula is
            l-|l,m> = sqrt((l + m) * (l - m + 1)) * h_bar * |l,m-1>.
        the actual coefficient should be multiplied by a factor of h_bar.

        :return: None
        """
        new_coefficients = self._init_coeff()
        for key, value in self._coeff.items():
            l, m, s = key
            key_new = (l, m-1, s)
            if key_new in self._allowed_qn:
                factor = sp.sqrt((l + m) * (l - m + 1))
                new_coefficients[key_new] = value * factor
        self._coeff = new_coefficients

    def l_z(self) -> None:
        """
        Apply the lz operator on the atomic orbital.

        As the formula is
            lz|l,m> = m * h_bar * |l,m>
        the actual coefficient should be multiplied by a factor of h_bar.

        :return: None
        """
        for key in self._coeff.keys():
            m = key[1]
            self._coeff[key] *= m

    def l_square(self) -> None:
        """
        Apply l2 operator on the atomic orbital.

        As the formula is
            l2|l,m> = l * (l + 1) * h_bar**2 * |l,m>
        the actual coefficient should be multiplied by a factor of h_bar**2.

        :return: None
        """
        for key in self._coeff.keys():
            l_i = key[0]
            self._coeff[key] *= (l_i * (l_i + 1))

    def s_plus(self) -> None:
        """
        Apply the s+ (not sigma+) operator on the atomic orbital.

        As the formula is
            s+|s> = h_bar * |s+2>
        the actual coefficient should be multiplied by a factor of h_bar. Since
        s is the quantum number of sigma_z rather than s_z, we increase it by 2,
        not 1.

        :return: None
        """
        new_coefficients = self._init_coeff()
        for key, value in self._coeff.items():
            l, m, s = key
            key_new = (l, m, s+2)
            if key_new in self._allowed_qn:
                new_coefficients[key_new] = value
        self._coeff = new_coefficients

    def s_minus(self) -> None:
        """
        Apply the s- (not sigma-) operator on the atomic orbital.

        As the formula is
            s-|s> = h_bar * |s-2>
        the actual coefficient should be multiplied by a factor of h_bar. Since
        s is the quantum number of sigma_z rather than s_z, we decrease it by 2,
        not 1.

        :return: None
        """
        new_coefficients = self._init_coeff()
        for key, value in self._coeff.items():
            l, m, s = key
            key_new = (l, m, s-2)
            if key_new in self._allowed_qn:
                new_coefficients[key_new] = value
        self._coeff = new_coefficients

    def s_z(self) -> None:
        """
        Apply the s_z (not sigma_z) operator on the atomic orbital.

        As the formula is
            s_z |s> = (+-0.5) * h_bar * |s>
        the actual coefficient should be multiplied by a factor of h_bar.

        :return: None
        """
        for key in self._coeff.keys():
            s = key[2]
            self._coeff[key] *= (s * sp.Rational(1, 2))

    def s_square(self) -> None:
        """
        Apply the s2 (not sigma2) operator on the atomic orbital.

        As the formula is
            s2 |s> = 0.75 * h_bar**2 * |s>
        the actual coefficient should be multiplied by a factor of h_bar**2.

        :return: None
        """
        for key in self._coeff.keys():
            self._coeff[key] *= sp.Rational(3, 4)

    def product(self, ket) -> c_type:
        """
        Evaluate the inner product <self|ket>.

        :param ket: the ket vector
        :return: the inner product
        """
        product = 0
        for key, value in self._coeff.items():
            product += value.conjugate() * ket[key]
        return product

    def mtxel(self, ket, operators: Iterable[str]) -> c_type:
        """
        Evaluate the matrix element <self|operators|ket>.

        :param ket: the ket vector
        :param operators: the operators
        :return: the matrix element
        :raises ValueError: if any operator is not in l+, l-, lz, l2 or their
            spin operator counterparts
        """
        ket_copy = deepcopy(ket)
        for op in operators:
            if op == "l+":
                ket_copy.l_plus()
            elif op == "l-":
                ket_copy.l_minus()
            elif op == "lz":
                ket_copy.l_z()
            elif op == "l2":
                ket_copy.l_square()
            elif op == "s+":
                ket_copy.s_plus()
            elif op == "s-":
                ket_copy.s_minus()
            elif op == "sz":
                ket_copy.s_z()
            elif op == "s2":
                ket_copy.s_square()
            else:
                raise ValueError(f"Illegal operator {op}")
        return self.product(ket_copy)


class SOC2:
    """
    Class for evaluating spin-orbital coupling terms.

    Attributes
    ----------
    _orbital_labels: Set[str]
        labels of atomic orbitals
    _spin_labels: Set[str]
        directions of spins
    _orbital_basis: Dict[str, AtomicOrbital]
        collection of atomic orbitals s, px, py, pz. etc
    """
    def __init__(self) -> None:
        self._orbital_labels = ("s", "px", "py", "pz",
                                "dxy", "dx2-y2", "dyz", "dzx", "dz2")
        self._spin_labels = ("up", "down")

        # Initialize atomic orbitals
        self._orbital_basis = dict()
        for orbital in self._orbital_labels:
            for spin in self._spin_labels:
                label = (orbital, spin)
                self._orbital_basis[label] = AtomicOrbital(l_max=2)

        # Reference:
        # https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
        c = sp.sqrt(sp.Rational(1, 2))
        ci = c * sp.I
        for spin in self._spin_labels:
            s = 1 if spin == "up" else -1
            # s state
            self._orbital_basis[("s", spin)][(0, 0, s)] = 1

            # p states
            self._orbital_basis[("py", spin)][(1, -1, s)] = ci
            self._orbital_basis[("py", spin)][(1, 1, s)] = ci
            self._orbital_basis[("pz", spin)][(1, 0, s)] = 1
            self._orbital_basis[("px", spin)][(1, -1, s)] = c
            self._orbital_basis[("px", spin)][(1, 1, s)] = -c

            # d states
            self._orbital_basis[("dxy", spin)][(2, -2, s)] = ci
            self._orbital_basis[("dxy", spin)][(2, 2, s)] = -ci
            self._orbital_basis[("dyz", spin)][(2, -1, s)] = ci
            self._orbital_basis[("dyz", spin)][(2, 1, s)] = ci
            self._orbital_basis[("dz2", spin)][(2, 0, s)] = 1
            self._orbital_basis[("dzx", spin)][(2, -1, s)] = c
            self._orbital_basis[("dzx", spin)][(2, 1, s)] = -c
            self._orbital_basis[("dx2-y2", spin)][(2, -2, s)] = c
            self._orbital_basis[("dx2-y2", spin)][(2, 2, s)] = c

    @staticmethod
    def _eval_soc(bra: AtomicOrbital, ket: AtomicOrbital) -> c_type:
        # l_dos_s = 0.5 * (l+*s- + l-*s+) + lz*sz
        p1 = bra.mtxel(ket, ["l+", "s-"])
        p2 = bra.mtxel(ket, ["l-", "s+"])
        p3 = bra.mtxel(ket, ["lz", "sz"])
        product = sp.Rational(1, 2) * (p1 + p2) + p3
        return product

    def print_soc_table(self, spin_i: str = "up", spin_j: str = "up") -> None:
        """
        Print SOC terms between orbital basis functions, for generating the
        coefficients in 'SOCTable' class.

        :param spin_i: spin direction of bra
        :param spin_j: spin direction kf ket
        :return: None
        :raises ValueError: if spin directions are illegal
        """
        for spin in (spin_i, spin_j):
            if spin not in self._spin_labels:
                raise ValueError(f"Illegal spin direction {spin}")
        soc_table = dict()
        for label_i in self._orbital_labels:
            for label_j in self._orbital_labels:
                bra = self._orbital_basis[(label_i, spin_i)]
                ket = self._orbital_basis[(label_j, spin_j)]
                soc = self._eval_soc(bra, ket)
                if abs(soc) > 1.0e-5:
                    soc_table[(label_i, label_j)] = soc
        print("Factor: h_bar**2")
        print(soc_table)

    def eval(self, label_i: str = "s",
             spin_i: str = "up",
             label_j: str = "s",
             spin_j: str = "down") -> c_type:
        """
        Evaluate the matrix element <i,s_i|l*s|j,s_j>.

        :param label_i: orbital label of bra
        :param spin_i: spin direction of bra
        :param label_j: orbital label of ket
        :param spin_j: spin direction of ket
        :return: matrix element in h_bar**2
        :raises ValueError: if orbital labels or spin directions are illegal
        """
        for label in (label_i, label_j):
            if label not in self._orbital_labels:
                raise ValueError(f"Illegal orbital label {label}")
        for spin in (spin_i, spin_j):
            if spin not in self._spin_labels:
                raise ValueError(f"Illegal spin direction {spin}")
        bra = self._orbital_basis[(label_i, spin_i)]
        ket = self._orbital_basis[(label_j, spin_j)]
        soc = self._eval_soc(bra, ket)
        return soc


def main():
    import tbplas as tb
    soc = tb.SOC()
    soc2 = SOC2()
    orbital_labels = ("s", "px", "py", "pz",
                      "dxy", "dx2-y2", "dyz", "dzx", "dz2")
    spin_labels = ("up", "down")

    # Check accuracy and speed
    timer = tb.Timer()
    for o1 in orbital_labels:
        for s1 in spin_labels:
            for o2 in orbital_labels:
                for s2 in spin_labels:
                    timer.tic("soc1")
                    v1 = soc.eval(o1, s1, o2, s2)
                    timer.toc("soc1")
                    timer.tic("soc2")
                    v2 = soc2.eval(o1, s1, o2, s2)
                    timer.toc("soc2")
                    if abs(v1 - v2) >= 1.0e-5:
                        print(v1 - v2)
    timer.report_total_time()

    # Print SOC table
    soc2.print_soc_table("up", "up")
    soc2.print_soc_table("up", "down")
    soc2.print_soc_table("down", "up")


if __name__ == "__main__":
    main()
