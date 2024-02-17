#ifndef ABPLAS_BASE_LATTICE_H
#define ABPLAS_BASE_LATTICE_H

#include <eigen3/Eigen/Dense>

namespace abplas {
namespace base {

/// @brief Convert fractional coordinates to Cartesian coordinates
/// @param[in] latticeVectors Cartesian coordinates of lattice vectors in columnwise-order
/// @param[in] fracCoordinates Fractional coordinates to convert in columnwise-order
/// @param[out] cartCoordinates Cartesian coordinates after conversion in columnwise-order
/// @note cartCoordinates will be resized to that of fracCoordinates automatically
void frac2cart(const Eigen::Matrix3d &latticeVectors, const Eigen::Matrix3Xd &fracCoordinates,
               Eigen::Matrix3Xd &cartCoordinates);

/// @brief Convert Cartesian coordinates to fractional coordinates
/// @param[in] latticeVectors Cartesian coordinates of lattice vectors in columnwise-order
/// @param[in] cartCoordinates Cartesian coordinates to convert in columnwise-order
/// @param[out] fracCoordinates Fractional coordinates after conversion in columnwise-order
/// @note fracCoordinates will be resized to that of cartCoordinates automatically
void cart2frac(const Eigen::Matrix3d &latticeVectors, const Eigen::Matrix3Xd &cartCoordinates,
               Eigen::Matrix3Xd &fracCoordinates);

/// @brief Calculate Cartesian coordinates of reciprocal lattice vectors
/// @param[in] realLatticeVectors Cartesian coordinates of real lattice in columnwise-order
/// @param[out] recipLatticeVectors Cartesian coordinates of reciprocal lattice in columnwise-order
void real2recip(const Eigen::Matrix3d &realLatticeVectors, Eigen::Matrix3d &recipLatticeVectors);

/// @brief Calculate Cartesian coordinates of reciprocal lattice vectors using another algorithm
/// @param[in] realLatticeVectors Cartesian coordinates of real lattice in columnwise-order
/// @param[out] recipLatticeVectors Cartesian coordinates of reciprocal lattice in columnwise-order
void real2recip2(const Eigen::Matrix3d &realLatticeVectors, Eigen::Matrix3d &recipLatticeVectors);

/// @brief Calculate lattice volume
/// @param[in] latticeVectors Cartesian coordinates of lattice vectors in columnwise-order
/// @return the volume
double calcLatticeVolume(const Eigen::Matrix3d &latticeVectors);

}
}

#endif