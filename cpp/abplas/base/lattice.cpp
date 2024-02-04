#include "lattice.h"


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

}
}