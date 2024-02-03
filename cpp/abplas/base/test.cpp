#include "input.h"
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
    Eigen::MatrixXd positions;
    sf.getPositions(elements, positions);
    print(elements);
    std::cout << positions.transpose() << "\n";
}

int main() {
    InputFile inp("INPUT");
    StructFile sf("STRUCT");
    testValue(inp);
    testNum(sf);
    testSpecies(sf);
    testLattice(sf);
    testPositions(sf);
}