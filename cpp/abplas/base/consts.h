///
/// CODATA recommended values of fundamental physical constants: 2018
/// Review of Modern Physics, 93, 025010 (2021).
///
/// @author: Yunhai Li (yhli)
///

#ifndef ABPLAS_BASE_CONSTS_H
#define ABPLAS_BASE_CONSTS_H

namespace abplas {
namespace base {

/// @brief PI
constexpr double PI = 3.14159265358979323846;

/// @brief twice of PI
constexpr double TWOPI = 2.0 * PI;

/// @brief Scaling factor from bohr to angstrom
/// Table XXXI, pp. 45
constexpr double BOHR2ANG = 0.529177210903;

/// @brief Scaling factor from bohr to nm
constexpr double BOHR2NM = 0.1 * BOHR2ANG;

/// @brief Scaling factor from angstrom to bohr
constexpr double ANG2BOHR = 1.0 / BOHR2ANG;

/// @brief Scaling factor from nm to bohr
constexpr double NM2BOHR = 1.0 / BOHR2NM;

/// @brief Scaling factor from Hartree to eV
/// Table XXXI, pp. 45
constexpr double HAR2EV = 27.211386245988;

/// @brief Scaling factor from eV to Hartree
constexpr double EV2HAR = 1.0 / HAR2EV;

/// @brief Scaling factor from Hartree to Joule
/// Table XXXI, pp. 45
constexpr double HAR2J = 4.3597447222071e-18;

/// @brief Scaling factor from Joule to Hartree
constexpr double J2HAR = 1.0 / HAR2J;

/// @brief Reduced Planck constant in J*s
/// Table XXXI, pp. 45
constexpr double H_BAR_JS = 1.054571817e-34;

/// @brief Reduced Planck constant in eV*s
/// Table XXXI, pp. 45
constexpr double H_BAR_EVS = 6.582119569e-16;

/// @brief Reduced Planck constant in Hartree*s
constexpr double H_BAR_HARS = H_BAR_EVS * EV2HAR;

/// @brief Boltzmann constant in J/K
/// Table XXXI, pp. 49
constexpr double KB_JK = 1.380649e-23;

/// @brief Boltzmann constant in eV/K
/// Table XXXI, pp. 49
constexpr double KB_EVK = 8.617333262e-5;

/// @brief Boltzmann const in Hartree/K
constexpr double KB_HARK = KB_EVK * EV2HAR;

/// @brief Scaling factor from a.u. time to femto-second
/// Table XXXIV, pp.51
constexpr double AU2FS = H_BAR_JS / HAR2J * 1.0e15;

/// @brief Scaling factor from femto-second to a.u. time
constexpr double FS2AU = 1.0 / AU2FS;

}
}

#endif