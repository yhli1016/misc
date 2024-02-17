#ifndef ABPLAS_BASE_INPUT_H
#define ABPLAS_BASE_INPUT_H

#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>

#include <eigen3/Eigen/Dense>

namespace abplas {
namespace base {

/// @brief Get the scaling factor from current unit to Bohr
/// @param[in] currentUnit kind of current unit, should be bohr, angstrom or nm,
///     see the regular expressions in the definition for more details
/// @return the scaling factor, fall back to 1.0 if the current unit cannot be
///     recognized
double getScaleFactorLength(const std::string &currentUnit);

/// @brief Get the scaling factor from current unit to Hartree
/// @param[in] currentUnit kind of current unit, should be hartree or ev,
///     see the regular expressions in the definition for more details
/// @return the scaling factor, fall back to 1.0 if the current unit cannot be
///     recognized
double getScaleFactorEnergy(const std::string &currentUnit);

/// @brief Get the scaling factor from current unit to a.u. time
/// @param[in] currentUnit kind of current unit, should be au or fs,
///     see the regular expressions in the definition for more details
/// @return the scaling factor, fall back to 1.0 if the current unit cannot be
///     recognized
double getScaleFactorTime(const std::string &currentUnit);

/// @brief Exception class for testing purpose
class ParameterNotFound: public std::exception {
    private:
        const char *m_Message;
    public:
        ParameterNotFound(const char *msg);
        virtual const char * what() const noexcept;
};

// Exception based on std::string
// class ParameterNotFound: public std::exception {
//     private:
//         std::string m_Message;
//     public:
//         ParameterNotFound(const std::string &msg);
//         virtual const char * what() const noexcept;
// };

/// @brief Base class for parsing input and struct files
class InputFile {
    public:
        explicit InputFile(const std::string &fileName);
        ~InputFile();

        /// @brief Get the value of key and put it in stringstream
        /// @param[in] keyPattern pattern of the key
        /// @param[out] valueStream stringstream object for holding the value
        /// @return whether the key is defined in input
        bool getValue(const std::string &keyPattern, std::stringstream &valueStream);

    protected:
        /// @brief Go to the first line of input file after parsing an item
        void rewind();

        /// @brief Get the index of line satisfying the pattern
        /// @param[in] pattern the pattern to search
        /// @return line number counted from 0, -1 if not found
        int indexPattern(const std::regex &pattern);

    protected:
        /// @brief file to parse
        std::ifstream m_InFile;
};

/// @brief Class for parsing struct files
class StructFile: public InputFile {
    public:
        explicit StructFile(const std::string &fileName);
        ~StructFile() = default;

        /// @brief Get the number of species in struct file
        /// @return number of species
        int getNumSpecies() const;

        /// @brief Get the number of atoms in struct file
        /// @return number of atoms
        int getNumAtoms() const;

        /// @brief Get the list of species in struct file
        /// @param[out] elements labels of species 
        /// @param[out] mass atomic masses of species
        /// @param[out] pseudoPots file names of pseudo-potentials
        /// @note There is not need to pre-allocate the lists. They will be resized automatically.
        ///     Also, there is not need to pre-allocate the lists. They will be resized automatically.
        void getSpecies(std::vector<std::string> &elements, Eigen::VectorXd &mass, std::vector<std::string> &pseudoPots);
        
        
        /// @brief Get the Cartesian coordinates of lattice vectors in bohr in struct file
        /// @param[out] lattice 3*3 matrix for holding the lattice vectors
        /// @note The coordinates are in columnwise-order, i.e., lattice(:,0) is the 0th lattice vector.
        void getLattice(Eigen::Matrix3d &lattice);
    
        /// @brief Get the symbols and Cartesian coordinates of atoms in bohr in struct file 
        /// @param[out] elements atomic symbols
        /// @param[out] positions atomic coordinates
        /// @note The coordinates are in columnwise-order, i.e., positions(:,0) corresponds to the 0th atom.
        ///     Also, there is not need to pre-allocate the lists and matrices. They will be resized automatically.
        void getPositions(std::vector<std::string> &elements, Eigen::Matrix3Xd &positions);

    private:
        /// @brief Check for format errors in struct file
        /// @note If any errors have been detected, the program will be terminated by calling exit,
        ///     since it makes no sense if the struct file is wrong.
        void checkSanity();

    private:
        // Regex patterns for detecting the header, body, and tail lines of species, lattice and
        // positions blocks in struct file
        std::regex m_SpeciesHead;
        std::regex m_SpeciesBody;
        std::regex m_SpeciesTail;
        std::regex m_LatticeHead;
        std::regex m_LatticeBody;
        std::regex m_LatticeTail;
        std::regex m_PositionsHead;
        std::regex m_PositionsBody;
        std::regex m_PositionsTail;

        // Number of species and atoms in struct file
        int m_NumSpecies;
        int m_NumAtoms;
};

}
}
#endif