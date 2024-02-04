#include "input.h"
#include "consts.h"
#include "lattice.h"

namespace {

void resetStringStream(std::stringstream &ss) {
    ss.clear();
    ss.str("");
}

double getScaleFactorLength(const std::string &currentUnit) {
    std::regex bohr("^\\s*b(ohr)?\\s*$", std::regex_constants::icase);
    std::regex ang("^\\s*a(ng)?(strom)?\\s*$", std::regex_constants::icase);
    std::regex nm("^\\s*n(ano)?m?(eter)?\\s*$", std::regex_constants::icase);
    double scale_factor = 1.0;
    if (std::regex_match(currentUnit, bohr)) {
        scale_factor = 1.0;
    } else if (std::regex_match(currentUnit, ang)) {
        scale_factor = abplas::base::ANG2BOHR;
    } else if (std::regex_match(currentUnit, nm)){
        scale_factor = abplas::base::NM2BOHR;
    } else {
        std::cout << "Unknown length unit " << currentUnit << "\n";
        std::exit(-1);
    }
    return scale_factor;
}

} // namespace

namespace abplas {
namespace base {

InputFile::InputFile(const std::string &fileName) {
    m_InFile.open(fileName, std::ios::in);
    if (!m_InFile.is_open()) {
        std::cout << "Failed to open " << fileName << "\n";
        std::exit(-1);
    }
}

InputFile::~InputFile() {
    if (m_InFile.is_open()) {
        m_InFile.close();
    }
}

void InputFile::rewind() {
    m_InFile.clear();
    m_InFile.seekg(0, std::ios::beg);
}

bool InputFile::findPattern(const std::regex &pattern) {
    rewind();
    bool status = false;
    std::string buffer = "\0";
    while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, pattern)) {
            status = true;
            break;
        }
    }
    return status;
}

int InputFile::indexPattern(const std::regex &pattern) {
    rewind();
    int idx = 0;
    std::string buffer = "\0";
    while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, pattern)) {
            break;
        } else {
            idx++;
        }
    }
    return idx;
}

bool InputFile::getValue(const std::string &keyPattern, std::stringstream &valueStream) {
    ::resetStringStream(valueStream);
    rewind();
    std::string buffer = "\0", key = "\0";
    std::regex pattern(keyPattern, std::regex_constants::icase);
    bool matched = false;
    while(std::getline(m_InFile, buffer)) {
        if(std::regex_match(buffer, pattern)) {
            valueStream << buffer;
            valueStream >> key;
            matched = true;
            break;
        }
    }
    return matched;
}

StructFile::StructFile(const std::string &fileName): InputFile(fileName) {
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

void StructFile::checkSanity() {
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

int StructFile::getNumSpecies() {
    int nl_start = indexPattern(m_SpeciesHead);
    int nl_end = indexPattern(m_SpeciesTail);
    int num_species = nl_end - nl_start - 1;
    return num_species;
}

int StructFile::getNumLattice() {
    int nl_start = indexPattern(m_LatticeHead);
    int nl_end = indexPattern(m_LatticeTail);
    int num_lattice = nl_end - nl_start - 1;
    return num_lattice;
}

int StructFile::getNumAtoms() {
    int nl_start = indexPattern(m_PositionsHead);
    int nl_end = indexPattern(m_PositionsTail);
    int num_atoms = nl_end - nl_start - 1;
    return num_atoms;
}

void StructFile::getSpecies(std::vector<std::string> &elements,
                                Eigen::VectorXd &mass,
                                std::vector<std::string> &pseudoPots) {
    // Further check struct sanity
    int num_species = getNumSpecies();
    if (num_species < 1) {
        std::cout << "Wrong number of atomic species!\n";
        std::exit(-1);
    }

    // Initialize parameters
    elements = std::vector<std::string>(num_species);
    mass = Eigen::VectorXd(num_species);
    pseudoPots = std::vector<std::string>(num_species);

    // Parsing
    rewind();
    std::string buffer = "\0";
    std::stringstream ss(buffer);
    bool in_block = false;
    int idx = 0;
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
            ::resetStringStream(ss);
            ss << buffer;
            ss >> elements[idx] >> mass(idx) >> pseudoPots[idx];
            idx++;
        }
    }
}

void StructFile::getLattice(Eigen::Matrix3d &lattice) {
    // Further check struct sanity
    int num_lattice = getNumLattice();
    if (num_lattice != 3) {
        std::cout << "Wrong number of lattice vectors!\n";
        std::exit(-1);
    }

    // Initialize parameters
    lattice = Eigen::Matrix3d();
    std::string length_unit = "bohr";

    // Parsing
    rewind();
    std::string buffer = "\0";
    std::stringstream ss(buffer);
    bool in_block = false;
    int idx = 0;
    while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, m_LatticeHead)) {
            ::resetStringStream(ss);
            ss << buffer;
            ss >> length_unit >> length_unit >> length_unit;
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
            ::resetStringStream(ss);
            ss << buffer;
            ss >> lattice(0, idx) >> lattice(1, idx) >> lattice(2, idx);
            idx++;
        }
    }

    // Unit and coordinates conversion
    lattice *= ::getScaleFactorLength(length_unit);
}

void StructFile::getPositions(std::vector<std::string> &elements,
                                  Eigen::Matrix3Xd &positions) {
    // Further check struct sanity
    int num_atoms = getNumAtoms();
    if (num_atoms < 1) {
        std::cout << "Wrong number of atoms!\n";
        std::exit(-1);
    }

    // Initialize parameters
    elements = std::vector<std::string>(num_atoms);
    positions = Eigen::Matrix3Xd(3, num_atoms);
    std::string length_unit = "bohr";

    // Parsing
    rewind();
    std::string buffer = "\0";
    std::stringstream ss(buffer);
    bool in_block = false;
    int idx = 0;
        while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, m_PositionsHead)) {
            ::resetStringStream(ss);
            ss << buffer;
            ss >> length_unit >> length_unit >> length_unit;
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
            ::resetStringStream(ss);
            ss << buffer;
            ss >> elements[idx] >> positions(0, idx) >> positions(1, idx) >> positions(2, idx);
            idx++;
        }
    }

    // Unit and coordinates conversion
    std::regex crystal("^\\s*c(rystal)?\\s*$", std::regex_constants::icase);
    if (std::regex_match(length_unit, crystal)) {
        Eigen::Matrix3d lattice;
        Eigen::Matrix3Xd positions_cart(3, num_atoms);
        getLattice(lattice);
        frac2cart(lattice, positions, positions_cart);
        positions = positions_cart;
    } else {
        positions *= ::getScaleFactorLength(length_unit);
    }
}

} // namespace
} // namespace