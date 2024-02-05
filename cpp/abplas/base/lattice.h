#ifndef ABPLAS_BASE_LATTICE_H
#define ABPLAS_BASE_LATTICE_H

#include <eigen3/Eigen/Dense>

namespace abplas {
namespace base {

void frac2cart(const Eigen::Matrix3d &latticeVectors, const Eigen::Matrix3Xd &fracCoordinates,
               Eigen::Matrix3Xd &cartCoordinates);

void cart2frac(const Eigen::Matrix3d &latticeVectors, const Eigen::Matrix3Xd &cartCoordinates,
               Eigen::Matrix3Xd &fracCoordinates);

void real2recip(const Eigen::Matrix3d &realLatticeVectors, Eigen::Matrix3d &recipLatticeVectors);

void real2recip2(const Eigen::Matrix3d &realLatticeVectors, Eigen::Matrix3d &recipLatticeVectors);

double calcLatticeVolume(const Eigen::Matrix3d &latticeVectors);

}
}

#endif