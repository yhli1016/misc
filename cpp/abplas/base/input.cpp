#include "input.h"

#include <iostream>
#include <cstdlib>

#include "consts.h"
#include "lattice.h"

namespace {

/// @brief Reset the stringstream object
/// @param[in,out] ss stringstream object to reset 
void resetStringStream(std::stringstream &ss) {
    ss.clear();
    ss.str("");
}

/// @brief Count the number of lines between header and tail
/// @param[in] idxStart line number of header, starting from either 0 or 1
/// @param[in] idxEnd line number of tail, starting from either 0 or 1
/// @return number of lines, with header and tail excluded
/// @note The header and tail lines are excluded, otherwise it should be +1 instead of -1.
int countLine(int idxStart, int idxEnd) {
    return idxEnd - idxStart - 1;
}

}

namespace abplas {
namespace base {

double getScaleFactorLength(const std::string &currentUnit) {
    std::regex bohr(R"(^\s*b(ohr)?\s*$)", std::regex_constants::icase);
    std::regex ang(R"(^\s*a(ng)?(strom)?\s*$)", std::regex_constants::icase);
    std::regex nm(R"(^\s*n(ano)?m?(eter)?\s*$)", std::regex_constants::icase);
    double scale_factor;
    if (std::regex_match(currentUnit, bohr)) {
        scale_factor = 1.0;
    } else if (std::regex_match(currentUnit, ang)) {
        scale_factor = abplas::base::ANG2BOHR;
    } else if (std::regex_match(currentUnit, nm)){
        scale_factor = abplas::base::NM2BOHR;
    } else {
        std::cout << "WARNING: unknown length unit: " << currentUnit << "\n";
        std::cout << "INFO: treating as Bohr\n";
        scale_factor = 1.0;
    }
    return scale_factor;
}

double getScaleFactorEnergy(const std::string &currentUnit) {
    std::regex har(R"(^\s*h(ar)?(tree)?\s*$)", std::regex_constants::icase);
    std::regex ev(R"(^\s*e(lectron)?v(olt)?\s*$)", std::regex_constants::icase);
    double scale_factor;
    if (std::regex_match(currentUnit, har)) {
        scale_factor = 1.0;
    } else if (std::regex_match(currentUnit, ev)) {
        scale_factor = abplas::base::EV2HAR;
    } else {
        std::cout << "WARNING: unknown energy unit: " << currentUnit << "\n";
        std::cout << "INFO: treating as Hartree\n";
        scale_factor = 1.0;
    }
    return scale_factor;
}

double getScaleFactorTime(const std::string &currentUnit) {
    std::regex au(R"(^\s*au\s*$)", std::regex_constants::icase);
    std::regex fs(R"(^\s*fs\s*$)", std::regex_constants::icase);
    double scale_factor;
    if (std::regex_match(currentUnit, au)) {
        scale_factor = 1.0;
    } else if (std::regex_match(currentUnit, fs)) {
        scale_factor = abplas::base::FS2AU;
    } else {
        std::cout << "WARNING: unknown time unit: " << currentUnit << "\n";
        std::cout << "INFO: treating as a.u.\n";
        scale_factor = 1.0;
    }
    return scale_factor;
}

ParameterNotFound::ParameterNotFound(const char * msg) {
    m_Message = msg;
}

const char * ParameterNotFound::what() const noexcept {
    return m_Message;
}

// ParameterNotFound::ParameterNotFound(const std::string &msg) {
//     m_Message = msg;
// }

// const char * ParameterNotFound::what() const noexcept {
//     return m_Message.data();
// }

InputFile::InputFile(const std::string &fileName) {
    m_InFile.open(fileName, std::ios::in);
    if (!m_InFile.is_open()) {
        std::cout << "ERROR: failed to open: " << fileName << "\n";
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

int InputFile::indexPattern(const std::regex &pattern) {
    rewind();
    bool found_pattern = false;
    int idx = 0;
    std::string buffer;
    while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, pattern)) {
            found_pattern = true;
            break;
        } else {
            idx++;
        }
    }
    if (!found_pattern) {
        idx = -1;
    }
    return idx;
}

bool InputFile::getValue(const std::string &keyPattern, std::stringstream &valueStream) {
    ::resetStringStream(valueStream);
    rewind();
    std::string buffer, key;
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
    m_SpeciesHead = std::regex(R"(^\s*begin\s+species\s*$)", std::regex_constants::icase);
    m_SpeciesBody = std::regex(R"(^\s*\w+\s+[\d\.]+\s+[\w\.\-]+\s*$)", std::regex_constants::icase);
    m_SpeciesTail = std::regex(R"(^\s*end\s+species\s*$)", std::regex_constants::icase);
    m_LatticeHead = std::regex(R"(^\s*begin\s+lattice\s+\w+\s*$)", std::regex_constants::icase);
    m_LatticeBody = std::regex(R"(^\s*[\d\.\-]+(\s+[\d\.\-]+){2}\s*$)", std::regex_constants::icase);
    m_LatticeTail = std::regex(R"(^\s*end\s+lattice\s*$)", std::regex_constants::icase);
    m_PositionsHead = std::regex(R"(^\s*begin\s+positions\s+\w+\s*$)", std::regex_constants::icase);
    m_PositionsBody = std::regex(R"(^\s*\w+(\s+[\d\.\-]+){3}\s*$)", std::regex_constants::icase);
    m_PositionsTail = std::regex(R"(^\s*end\s+positions\s*$)", std::regex_constants::icase);
    m_NumSpecies = 0;
    m_NumAtoms = 0;
    checkSanity();
}

void StructFile::checkSanity() {
    // atomic species
    int idx_start = indexPattern(m_SpeciesHead);
    int idx_end = indexPattern(m_SpeciesTail);
    m_NumSpecies = ::countLine(idx_start, idx_end);
    if (idx_start == -1) {
        std::cout << "ERROR: begin species not found\n";
        std::exit(-1);
    }
    if (idx_end == -1) {
        std::cout << "ERROR: end species not found\n";
        std::exit(-1);
    }
    if (m_NumSpecies < 1) {
        std::cout << "ERROR: wrong number of species\n";
        std::exit(-1);
    }

    // lattice vectors
    idx_start = indexPattern(m_LatticeHead);
    idx_end = indexPattern(m_LatticeTail);
    int num_lattice = ::countLine(idx_start, idx_end);
    if (idx_start == -1) {
        std::cout << "ERROR: begin lattice not found\n";
        std::exit(-1);
    }
    if (idx_end == -1) {
        std::cout << "ERROR: end lattice not found\n";
        std::exit(-1);
    } 
    if (num_lattice != 3) {
        std::cout << "ERROR: wrong number of lattice vectors\n";
        std::exit(-1);
    }

    // atomic positions
    idx_start = indexPattern(m_PositionsHead);
    idx_end = indexPattern(m_PositionsTail);
    m_NumAtoms = ::countLine(idx_start, idx_end);
    if (idx_start == -1) {
        std::cout << "ERROR: begin positions not found\n";
        std::exit(-1);
    }
    if (idx_end == -1) {
        std::cout << "ERROR: end positions not found\n";
        std::exit(-1);
    }
    if (m_NumAtoms < 1) {
        std::cout << "ERROR: wrong number of atoms\n";
        std::exit(-1);
    }
}

int StructFile::getNumSpecies() const {
    return m_NumSpecies;
}

int StructFile::getNumAtoms() const {
    return m_NumAtoms;
}

void StructFile::getSpecies(std::vector<std::string> &elements,
                                Eigen::VectorXd &mass,
                                std::vector<std::string> &pseudoPots) {
    // Initialize parameters
    elements = std::vector<std::string>(m_NumSpecies);
    mass = Eigen::VectorXd(m_NumSpecies);
    pseudoPots = std::vector<std::string>(m_NumSpecies);

    // Parsing
    rewind();
    std::string buffer;
    std::stringstream ss;
    bool in_block = false;
    int idx = 0;
    while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, m_SpeciesHead)) {
            in_block = true;
            continue;
        } else if (std::regex_match(buffer, m_SpeciesTail)) {
            break;
        } else if (in_block) {
            if (!std::regex_match(buffer, m_SpeciesBody)) {
                std::cout << "ERROR: illegal line in species: " << buffer << "\n";
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
    // Initialize parameters
    lattice = Eigen::Matrix3d();
    std::string length_unit = "bohr";

    // Parsing
    rewind();
    std::string buffer;
    std::stringstream ss;
    bool in_block = false;
    int idx = 0;
    while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, m_LatticeHead)) {
            ::resetStringStream(ss);
            ss << buffer;
            ss >> length_unit >> length_unit >> length_unit;
            in_block = true;
            continue;
        } else if (std::regex_match(buffer, m_LatticeTail)) {
            break;
        } else if (in_block) {
            if (!std::regex_match(buffer, m_LatticeBody)) {
                std::cout << "ERROR: illegal line in lattice: " << buffer << "\n";
                std::exit(-1);
            }
            ::resetStringStream(ss);
            ss << buffer;
            ss >> lattice(0, idx) >> lattice(1, idx) >> lattice(2, idx);
            idx++;
        }
    }

    // Unit and coordinates conversion
    lattice *= getScaleFactorLength(length_unit);
}

void StructFile::getPositions(std::vector<std::string> &elements,
                                  Eigen::Matrix3Xd &positions) {
    // Initialize parameters
    elements = std::vector<std::string>(m_NumAtoms);
    positions = Eigen::Matrix3Xd(3, m_NumAtoms);
    std::string length_unit = "bohr";

    // Parsing
    rewind();
    std::string buffer;
    std::stringstream ss;
    bool in_block = false;
    int idx = 0;
        while (std::getline(m_InFile, buffer)) {
        if (std::regex_match(buffer, m_PositionsHead)) {
            ::resetStringStream(ss);
            ss << buffer;
            ss >> length_unit >> length_unit >> length_unit;
            in_block = true;
            continue;
        } else if (std::regex_match(buffer, m_PositionsTail)) {
            break;
        } else if (in_block) {
            if (!std::regex_match(buffer, m_PositionsBody)) {
                std::cout << "ERROR: illegal line in positions: " << buffer << "\n";
                std::exit(-1);
            }
            ::resetStringStream(ss);
            ss << buffer;
            ss >> elements[idx] >> positions(0, idx) >> positions(1, idx) >> positions(2, idx);
            idx++;
        }
    }

    // Unit and coordinates conversion
    std::regex crystal(R"(^\s*c(rystal)?\s*$)", std::regex_constants::icase);
    if (std::regex_match(length_unit, crystal)) {
        Eigen::Matrix3d lattice;
        Eigen::Matrix3Xd positions_cart(3, m_NumAtoms);
        getLattice(lattice);
        frac2cart(lattice, positions, positions_cart);
        positions = positions_cart;
    } else {
        positions *= getScaleFactorLength(length_unit);
    }
}


}
}