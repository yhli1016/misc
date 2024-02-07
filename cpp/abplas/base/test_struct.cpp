#include "input.h"
#include "lattice.h"
using namespace abplas::base;

template<typename T>
void print(T &elements) {
    for (const auto &i: elements) {
        std::cout << i << "\n";
    }
}

void testNum(StructFile &sf) {
    std::cout << "-------- Numbers --------\n";
    std::cout << "Number of species " << sf.getNumSpecies() << "\n";
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

int main() {
    StructFile sf("STRUCT");
    testNum(sf);
    testSpecies(sf);
    testLattice(sf);
    testPositions(sf);
}