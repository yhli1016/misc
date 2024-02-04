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
    int nx = 0, ny = 0, nz = 0;
    if (inp.getValue("^\\s*k_grid(\\s+\\d+){3}\\s*$", ss)) {
        ss >> nx >> ny >> nz;
        std::cout << "k_grid = " << nx << " " << ny << " " << nz << "\n";
    } else {
        std::cout << "k_grid not found\n";
    }
    std::string mix_algo = "\0";
    if (inp.getValue("^\\s*mix_algo\\s+\\w+\\s*$", ss)) {
        ss >> mix_algo;
        std::cout << "mix_algo = " << mix_algo << "\n";
    } else {
        std::cout << "mix_algo not found\n";
    }
    double sigma = 0.0;
    if (inp.getValue("^\\s*sigma\\s+[\\d\\.\\-]+\\s*$", ss)) {
        ss >> sigma;
        std::cout << "sigma = " << sigma << "\n";
    } else {
        std::cout << "sigma not found\n";
    }
}

int main() {
    InputFile inp("INPUT");
    testValue(inp);
}