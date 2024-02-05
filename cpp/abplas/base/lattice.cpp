#include "lattice.h"

#include <cmath>

#include "consts.h"


namespace abplas{
namespace base {

void frac2cart(const Eigen::Matrix3d &latticeVectors, const Eigen::Matrix3Xd &fracCoordinates,
               Eigen::Matrix3Xd &cartCoordinates) {
    int num_col = fracCoordinates.cols();
    cartCoordinates = Eigen::Matrix3Xd(3, num_col);
    for(int idx = 0; idx < num_col; idx++) {
        cartCoordinates.col(idx) = latticeVectors * fracCoordinates.col(idx);
    }
}

void cart2frac(const Eigen::Matrix3d &latticeVectors, const Eigen::Matrix3Xd &cartCoordinates,
               Eigen::Matrix3Xd &fracCoordinates) {
    int num_col = cartCoordinates.cols();
    fracCoordinates = Eigen::Matrix3Xd(3, num_col);
    Eigen::Matrix3d conv_mat = latticeVectors.inverse();
    for(int idx = 0; idx < num_col; idx++) {
        fracCoordinates.col(idx) = conv_mat * cartCoordinates.col(idx);
    }
}

void real2recip(const Eigen::Matrix3d &realLatticeVectors, Eigen::Matrix3d &recipLatticeVectors) {
    recipLatticeVectors = Eigen::Matrix3d();
    Eigen::Matrix3d coeff = realLatticeVectors.transpose();
    Eigen::Matrix3d product = abplas::base::TWOPI * Eigen::Matrix3d::Identity(3, 3);
    for (int idx = 0; idx < 3; idx++) {
        recipLatticeVectors.col(idx) = coeff.colPivHouseholderQr().solve(product.col(idx));
    }
}

void real2recip2(const Eigen::Matrix3d &realLatticeVectors, Eigen::Matrix3d &recipLatticeVectors) {
    recipLatticeVectors = Eigen::Matrix3d();
    Eigen::Vector3d a0 = realLatticeVectors.col(0);
    Eigen::Vector3d a1 = realLatticeVectors.col(1);
    Eigen::Vector3d a2 = realLatticeVectors.col(2);
    double volume = a0.cross(a1).dot(a2);
    recipLatticeVectors.col(0) = abplas::base::TWOPI / volume * a1.cross(a2);
    recipLatticeVectors.col(1) = abplas::base::TWOPI / volume * a2.cross(a0);
    recipLatticeVectors.col(2) = abplas::base::TWOPI / volume * a0.cross(a1);
}

double calcLatticeVolume(const Eigen::Matrix3d &latticeVectors) {
    Eigen::Vector3d a0 = latticeVectors.col(0);
    Eigen::Vector3d a1 = latticeVectors.col(1);
    Eigen::Vector3d a2 = latticeVectors.col(2);
    double volume = std::abs(a0.cross(a1).dot(a2));
    return volume;
}

}
}