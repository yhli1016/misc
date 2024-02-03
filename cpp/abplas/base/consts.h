#ifndef ABPLAS_BASE_CONSTS_H
#define ABPLAS_BASE_CONSTS_H

namespace abplas {
namespace base {

/// @brief multiplication factor from bohr to angstrom
constexpr double BOHR2ANG = 0.5291772109253;

/// @brief multiplication factor from bohr to nm
constexpr double BOHR2NM = 0.1 * BOHR2ANG;

/// @brief multiplication factor from angstrom to bohr
constexpr double ANG2BOHR = 1.0 / BOHR2ANG;

/// @brief multiplcation factor from nm to bohr
constexpr double NM2BOHR = 1.0 / BOHR2NM;

}
}

#endif