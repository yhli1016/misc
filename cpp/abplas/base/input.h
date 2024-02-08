#ifndef ABPLAS_BASE_INPUT_H
#define ABPLAS_BASE_INPUT_H

#include <fstream>
#include <sstream>
#include <regex>
#include <eigen3/Eigen/Dense>

namespace abplas {
namespace base {

double getScaleFactorLength(const std::string &currentUnit);

double getScaleFactorEnergy(const std::string &currentUnit);

double getScaleFactorTime(const std::string &currentUnit);

class InputFile {
    public:
        explicit InputFile(const std::string &fileName);
        ~InputFile();
        bool getValue(const std::string &keyPattern, std::stringstream &valueStream);
    protected:
        void rewind();
        int indexPattern(const std::regex &pattern);
    protected:
        std::ifstream m_InFile;
};

class StructFile: public InputFile {
    public:
        explicit StructFile(const std::string &fileName);
        ~StructFile() = default;
        int getNumSpecies() const;
        int getNumAtoms() const;
        void getSpecies(std::vector<std::string> &elements, Eigen::VectorXd &mass, std::vector<std::string> &pseudoPots);
        void getLattice(Eigen::Matrix3d &lattice);
        void getPositions(std::vector<std::string> &elements, Eigen::Matrix3Xd &positions);
    private:
        void checkSanity();
    private:
        std::regex m_SpeciesHead;
        std::regex m_SpeciesBody;
        std::regex m_SpeciesTail;
        std::regex m_LatticeHead;
        std::regex m_LatticeBody;
        std::regex m_LatticeTail;
        std::regex m_PositionsHead;
        std::regex m_PositionsBody;
        std::regex m_PositionsTail;
        int m_NumSpecies;
        int m_NumAtoms;
};

}
}
#endif