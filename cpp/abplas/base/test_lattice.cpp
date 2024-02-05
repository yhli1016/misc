#include <iostream>

#include "lattice.h"
using namespace abplas::base;

void testConvert() {
    std::cout << "-------- Coordinates conversion --------\n";
    Eigen::Matrix3d lat_vec {{1.0, 1.0, -1.0},
                             {2.0, -1.5, 1.0},
                             {0.0, 0.1, 2.1}};
    Eigen::Matrix3Xd frac_coord {{0.5, 0.0, 0.0, 0.7, 6.0},
                                {1.0, 1.0, 0.0, 0.1, -0.1},
                                {0.0, 1.0, 1.0, -1.2, 0.7}};
    // Print origin matrices
    Eigen::Matrix3Xd cart_coord;
    std::cout << "lat_vec:\n";
    std::cout << lat_vec << "\n";
    std::cout << "frac_coord:\n";
    std::cout << frac_coord << "\n";

    // Print cartesian coordinates
    frac2cart(lat_vec, frac_coord, cart_coord);
    std::cout << "cart_coord:\n";
    std::cout << cart_coord << "\n";

    // Transform back to fractional coordinates and compare
    Eigen::Matrix3Xd frac_coord2;
    cart2frac(lat_vec, cart_coord, frac_coord2);
    Eigen::Matrix3Xd delta = frac_coord - frac_coord2;
    std::cout << "abs(delta):\n";
    std::cout << delta.array().abs() << "\n";
}

void testRealRecip() {
    std::cout << "-------- Coordinates conversion --------\n";
    Eigen::Matrix3d lat_vec {{1.0, 1.0, -1.0},
                             {2.0, -1.5, 1.0},
                             {0.0, 0.1, 2.1}};

    // volume of lattice
    double volume = abplas::base::calcLatticeVolume(lat_vec);
    std::cout << "Volume = " << volume << "\n";

    // reciprocal lattice vectors
    Eigen::Matrix3d recip_vec;
    abplas::base::real2recip(lat_vec, recip_vec);
    std::cout << "algo1:\n";
    std::cout << recip_vec << "\n";
    std::cout << "algo2:\n";
    abplas::base::real2recip2(lat_vec, recip_vec);
    std::cout << recip_vec << "\n";
}

int main() {
    testConvert();
    testRealRecip();
}
