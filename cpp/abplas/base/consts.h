#ifndef ABPLAS_BASE_CONSTS_H
#define ABPLAS_BASE_CONSTS_H

namespace abplas {
namespace base {

/// @brief PI
constexpr double PI = 3.14159265358979323846;

/// @brief twice of PI
constexpr double TWOPI = 2.0 * PI;

/// @brief Scaling factor from bohr to angstrom
constexpr double BOHR2ANG = 0.529177210671212;

/// @brief Scaling factor from bohr to nm
constexpr double BOHR2NM = 0.1 * BOHR2ANG;

/// @brief Scaling factor from angstrom to bohr
constexpr double ANG2BOHR = 1.0 / BOHR2ANG;

/// @brief Scaling factor from nm to bohr
constexpr double NM2BOHR = 1.0 / BOHR2NM;

/// @brief Scaling factor from Hartree to eV
constexpr double HAR2EV = 27.21138602;

/// @brief Scaling factor from eV to Hartree
constexpr double EV2HAR = 1.0 / HAR2EV;

/// @brief Scaling factor from Hartree to Joule
constexpr double HAR2J = 4.3597447222071e-18;

/// @brief Scaling factor from Joule to Hartree
constexpr double J2HAR = 1.0 / HAR2J;

/// @brief Reduced Planck constant in J*s
constexpr double H_BAR_JS = 1.054571817e-34;

/// @brief Reduced Planck constant in eV*s
constexpr double H_BAR_EVS = 6.582119569e-16;

/// @brief Reduced Planck constant in Hartree*s
constexpr double H_BAR_HARS = H_BAR_EVS * EV2HAR;

/// @brief Boltzmann constant in J/K
constexpr double KB_JK = 1.380649e-23;

/// @brief Boltzmann constant in eV/K
constexpr double KB_EVK = 8.617333262e-5;

/// @brief Boltzmann const in Hartree/K
constexpr double KB_HARK = KB_EVK * EV2HAR;

/// @brief Scaling factor from a.u. time to femto-second
constexpr double AU2FS = H_BAR_JS / HAR2J * 1.0e15;

/// @brief Scaling factor from femto-second to a.u. time
constexpr double FS2AU = 1.0 / AU2FS;

}
}

#endif