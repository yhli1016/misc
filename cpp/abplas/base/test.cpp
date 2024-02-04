#include "input.h"
#include "lattice.h"
using namespace abplas::base;

template<typename T>
void print(T &elements) {
    for (const auto &i: elements) {
        std::cout << i << "\n";
    }
}

void testValue(InputFile &inp) {
    std::cout << "-------- Key value --------\n";
    std::stringstream ss;
    int nx, ny, nz;
    if (inp.getValue("^\\s*k_grid(\\s+\\d+){3}\\s*$", ss)) {
        ss >> nx >> ny >> nz;
        std::cout << "k_grid = " << nx << " " << ny << " " << nz << "\n";
    } else {
        std::cout << "k_grid not found\n";
    }
    std::string mix_algo;
    if (inp.getValue("^\\s*mix_algo\\s+\\w+\\s*$", ss)) {
        ss >> mix_algo;
        std::cout << "mix_algo = " << mix_algo << "\n";
    }
}

void testNum(StructFile &sf) {
    std::cout << "-------- Numbers --------\n";
    std::cout << "Number of species " << sf.getNumSpecies() << "\n";
    std::cout << "Number of lattice vectors " << sf.getNumLattice() << "\n";
    std::cout << "Number of atoms " << sf.getNumAtoms() << "\n";
}

void testSpecies(StructFile &sf) {
    std::cout << "-------- Species --------\n";
    std::vector<std::string> elements;
    Eigen::VectorXd mass;
    std::vector<std::string> pseudo_pots;
    sf.getSpecies(elements, mass, pseudo_pots);
    print(elements);
    std::cout << mass << "\n";
    print(pseudo_pots);
}

void testLattice(StructFile &sf) {
    std::cout << "-------- Lattice vectors --------\n";
    Eigen::Matrix3d lattice;
    sf.getLattice(lattice);
    std::cout << lattice.transpose() << "\n";
}

void testPositions(StructFile &sf) {
    std::cout << "-------- Positions --------\n";
    std::vector<std::string> elements;
    Eigen::Matrix3Xd positions;
    sf.getPositions(elements, positions);
    print(elements);
    std::cout << positions.transpose() << "\n";
}

void testConvert() {
    std::cout << "-------- Coordinates conversion --------\n";
    Eigen::Matrix3d lat_vec {{1.0, 1.0, -1.0},
                             {2.0, -1.5, 1.0},
                             {0.0, 0.1, 2.1}};
    Eigen::Matrix3Xd frac_coord {{0.5, 0.0, 0.0, 0.7, 6.0},
                                {1.0, 1.0, 0.0, 0.1, -0.1},
                                {0.0, 1.0, 1.0, -1.2, 0.7}};
    Eigen::Matrix3Xd cart_coord;
    std::cout << lat_vec << "\n\n";
    std::cout << frac_coord << "\n\n";
    frac2cart(lat_vec, frac_coord, cart_coord);
    std::cout << cart_coord << "\n\n";
    cart2frac(lat_vec, cart_coord, frac_coord);
    std::cout << frac_coord << "\n\n";
}

int main() {
    InputFile inp("INPUT");
    StructFile sf("STRUCT");
    testValue(inp);
    testNum(sf);
    testSpecies(sf);
    testLattice(sf);
    testPositions(sf);
    testConvert();
}