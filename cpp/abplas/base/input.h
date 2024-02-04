#ifndef ABPLAS_BASE_INPUT_H
#define ABPLAS_BASE_INPUT_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <regex>
#include <cstdlib>
#include <eigen3/Eigen/Dense>

namespace abplas {
namespace base {

class InputFile {
    public:
        InputFile(const std::string &fileName);
        ~InputFile();
        bool getValue(const std::string &keyPattern, std::stringstream &valueStream);
    protected:
        void rewind();
        bool findPattern(const std::regex &pattern);
        int indexPattern(const std::regex &pattern);
    protected:
        std::ifstream m_InFile;
};

class StructFile: public InputFile {
    public:
        StructFile(const std::string &fileName);
        ~StructFile(){};
        int getNumSpecies();
        int getNumLattice();
        int getNumAtoms();
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
};

} // namespace
} // namespace
#endif