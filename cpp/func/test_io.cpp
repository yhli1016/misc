#include <fstream>
#include <iostream>
#include <sstream>
#include <regex>
#include <cstdlib>
#include <assert.h>
#include <eigen3/Eigen/Dense>


class StructParser {
    private:
        std::ifstream m_InFile;
        std::regex m_SpeciesHead;
        std::regex m_SpeciesBody;
        std::regex m_SpeciesTail;
        std::regex m_LatticeHead;
        std::regex m_LatticeBody;
        std::regex m_LatticeTail;
        std::regex m_PositionsHead;
        std::regex m_PositionsBody;
        std::regex m_PositionsTail;
    public:
        StructParser();
        ~StructParser();
        void checkSanity();
        int getNumSpecies();
        int getNumLattice();
        int getNumAtoms();
        void parseSpecies(std::vector<std::string> &elements, Eigen::VectorXd &mass, std::vector<std::string> &pseudo_pots);
        void parseLattice(Eigen::Matrix3d &lattice);
        void parsePositions(std::vector<std::string> &elements, Eigen::MatrixXd &positions);
    private:
        void rewind();
        bool findPattern(const std::regex &pattern);
        int indexPattern(const std::regex &pattern);
};

StructParser::StructParser() {
    m_InFile.open("STRUCT", std::ios::in);
    if (!m_InFile.is_open()) {
        std::cout << "Failed to open STRUCT\n";
        std::exit(-1);
    }
    m_SpeciesHead = std::regex("^\\s*begin\\s+species\\s*$", std::regex_constants::icase);
    m_SpeciesBody = std::regex("^\\s*\\w+\\s+[\\d\\.]+\\s+[\\w\\.\\-]+\\s*$", std::regex_constants::icase);
    m_SpeciesTail = std::regex("^\\s*end\\s+species\\s*$", std::regex_constants::icase);
    m_LatticeHead = std::regex("^\\s*begin\\s+lattice\\s+\\w+\\s*$", std::regex_constants::icase);
    m_LatticeBody = std::regex("^\\s*[\\d\\.\\-]+(\\s+[\\d\\.\\-]+){2}\\s*$", std::regex_constants::icase);
    m_LatticeTail = std::regex("^\\s*end\\s+lattice\\s*$", std::regex_constants::icase);
    m_PositionsHead = std::regex("^\\s*begin\\s+positions\\s+\\w+\\s*$", std::regex_constants::icase);
    m_PositionsBody = std::regex("^\\s*\\w+(\\s+[\\d\\.\\-]+){3}\\s*$", std::regex_constants::icase);
    m_PositionsTail = std::regex("^\\s*end\\s+positions\\s*$", std::regex_constants::icase);
    checkSanity();
}

StructParser::~StructParser() {
    if (m_InFile.is_open()) {
        m_InFile.close();
    }
}

void StructParser::rewind() {
    m_InFile.clear();
    m_InFile.seekg(0, std::ios::beg);
}

bool StructParser::findPattern(const std::regex &pattern) {
    rewind();
    bool status = false;
    std::string buffer;
    while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, pattern)) {
            status = true;
            break;
        }
    }
    return status;
}

int StructParser::indexPattern(const std::regex &pattern) {
    rewind();
    int idx = 0;
    std::string buffer;
    while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, pattern)) {
            break;
        } else {
            idx++;
        }
    }
    return idx;
}

void StructParser::checkSanity() {
    // atomic species
    if (!findPattern(m_SpeciesHead)) {
        std::cout << "Begin species not found!\n";
        std::exit(-1);
    }
    if (!findPattern(m_SpeciesTail)) {
        std::cout << "End species not found!\n";
        std::exit(-1);
    }

    // lattice vectors
    if (!findPattern(m_LatticeHead)) {
        std::cout << "Begin lattice not found!\n";
        std::exit(-1);
    }
    if (!findPattern(m_LatticeTail)) {
        std::cout << "End lattice not found!\n";
        std::exit(-1);
    }

    // atomic positions
    if (!findPattern(m_PositionsHead)) {
        std::cout << "Begin positions not found!\n";
        std::exit(-1);
    }
    if (!findPattern(m_PositionsTail)) {
        std::cout << "End positions not found!\n";
        std::exit(-1);
    }
}

int StructParser::getNumSpecies() {
    int nl_start = indexPattern(m_SpeciesHead);
    int nl_end = indexPattern(m_SpeciesTail);
    int num_species = nl_end - nl_start - 1;
    return num_species;
}

int StructParser::getNumLattice() {
    int nl_start = indexPattern(m_LatticeHead);
    int nl_end = indexPattern(m_LatticeTail);
    int num_lattice = nl_end - nl_start - 1;
    return num_lattice;
}

int StructParser::getNumAtoms() {
    int nl_start = indexPattern(m_PositionsHead);
    int nl_end = indexPattern(m_PositionsTail);
    int num_atoms = nl_end - nl_start - 1;
    return num_atoms;
}

void StructParser::parseSpecies(std::vector<std::string> &elements,
                                Eigen::VectorXd &mass,
                                std::vector<std::string> &pseudo_pots) {
    // Further check struct sanity
    int num_species = getNumSpecies();
    if (num_species < 1) {
        std::cout << "Wrong number of atomic species!\n";
        std::exit(-1);
    }

    // Initialize parameters
    elements = std::vector<std::string>(num_species);
    mass = Eigen::VectorXd(num_species);
    pseudo_pots = std::vector<std::string>(num_species);

    // Parsing
    rewind();
    std::string buffer;
    bool in_block = false;
    int row = 0;
    while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, m_SpeciesHead)) {
            in_block = true;
            continue;
        }
        if (std::regex_match(buffer, m_SpeciesTail)) {
            in_block = false;
            break;
        }
        if (in_block) {
            if (!std::regex_match(buffer, m_SpeciesBody)) {
                std::cout << "Illegal line in species: " << buffer << "\n";
                std::exit(-1);
            }
            std::stringstream ss(buffer);
            ss >> elements[row] >> mass(row) >> pseudo_pots[row];
            row++;
        }
    }
}

void StructParser::parseLattice(Eigen::Matrix3d &lattice) {
    // Further check struct sanity
    int num_lattice = getNumLattice();
    if (num_lattice != 3) {
        std::cout << "Wrong number of lattice vectors!\n";
        std::exit(-1);
    }

    // Initialize parameters
    lattice = Eigen::Matrix3d();

    // Parsing
    rewind();
    std::string buffer;
    bool in_block = false;
    int row = 0;
    while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, m_LatticeHead)) {
            in_block = true;
            continue;
        }
        if (std::regex_match(buffer, m_LatticeTail)) {
            in_block = false;
            break;
        }
        if (in_block) {
            if (!std::regex_match(buffer, m_LatticeBody)) {
                std::cout << "Illegal line in lattice: " << buffer << "\n";
                std::exit(-1);
            }
            std::stringstream ss(buffer);
            ss >> lattice(row, 0) >> lattice(row, 1) >> lattice(row, 2);
            row++;
        }
    }
}

void StructParser::parsePositions(std::vector<std::string> &elements,
                                  Eigen::MatrixXd &positions) {
    // Further check struct sanity
    int num_atoms = getNumAtoms();
    if (num_atoms < 1) {
        std::cout << "Wrong number of atoms!\n";
        std::exit(-1);
    }

    // Initialize parameters
    elements = std::vector<std::string>(num_atoms);
    positions = Eigen::MatrixX3d(num_atoms, 3);

    // Parsing
    rewind();
    std::string buffer;
    bool in_block = false;
    int row = 0;
        while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, m_PositionsHead)) {
            in_block = true;
            continue;
        }
        if (std::regex_match(buffer, m_PositionsTail)) {
            in_block = false;
            break;
        }
        if (in_block) {
            if (!std::regex_match(buffer, m_PositionsBody)) {
                std::cout << "Illegal line in positions: " << buffer << "\n";
                std::exit(-1);
            }
            std::stringstream ss(buffer);
            ss >> elements[row] >> positions(row, 0) >> positions(row, 1) >> positions(row, 2);
            row++;
        }
    }
}

template<typename T>
void print(T &elements) {
    for (const auto &i: elements) {
        std::cout << i << "\n";
    }
}

void testNum(StructParser &sp) {
    std::cout << "-------- Numbers --------\n";
    std::cout << "Number of species " << sp.getNumSpecies() << "\n";
    std::cout << "Number of lattice vectors " << sp.getNumLattice() << "\n";
    std::cout << "Number of atoms " << sp.getNumAtoms() << "\n";
}

void testSpecies(StructParser &sp) {
    std::cout << "-------- Species --------\n";
    std::vector<std::string> elements;
    Eigen::VectorXd mass;
    std::vector<std::string> pseudo_pots;
    sp.parseSpecies(elements, mass, pseudo_pots);
    print(elements);
    std::cout << mass << "\n";
    print(pseudo_pots);
}

void testLattice(StructParser &sp) {
    std::cout << "-------- Lattice vectors --------\n";
    Eigen::Matrix3d lattice;
    sp.parseLattice(lattice);
    std::cout << lattice << "\n";
}

void testPositions(StructParser &sp) {
    std::cout << "-------- Positions --------\n";
    std::vector<std::string> elements;
    Eigen::MatrixXd positions;
    sp.parsePositions(elements, positions);
    print(elements);
    std::cout << positions << "\n";
}

int main() {
    StructParser sp;
    testNum(sp);
    testSpecies(sp);
    testLattice(sp);
    testPositions(sp);
}