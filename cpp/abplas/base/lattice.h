#ifndef ABPLAS_BASE_LATTICE_H
#define ABPLAS_BASE_LATTICE_H

#include <eigen3/Eigen/Dense>

namespace abplas {
namespace base {

void frac2cart(const Eigen::Matrix3d &latticeVectors, const Eigen::Matrix3Xd &fracCoordinates,
               Eigen::Matrix3Xd &cartCoordinates);

void cart2frac(const Eigen::Matrix3d &latticeVectors, const Eigen::Matrix3Xd &cartCoordinates,
               Eigen::Matrix3Xd &fracCoordinates);

}
}

#endif