#include "input.h"
using namespace abplas::base;

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

void testUnits(InputFile &inp) {
    std::cout << "-------- Length --------\n";
    std::stringstream ss;
    double test = 0.0;
    std::string unit = "\0";
    std::string head = "^\\s*", tail = "\\s+[\\d\\.\\-]+\\s+\\w+\\s*$";

    // Test length units
    std::vector<std::string> vec = {"len_a", "len_b", "len_c", "len_d"};
    for (const auto &s: vec) {
        if (inp.getValue(head+s+tail, ss)) {
            ss >> test >> unit;
            std::cout << s << " = " << test * getScaleFactorLength(unit) << "\n";
        } else {
            std::cout << s << " not found\n";
        }
    }

    // Test energy units
    vec = {"energy_a", "energy_b", "energy_c"};
    for (const auto &s: vec) {
        if (inp.getValue(head+s+tail, ss)) {
            ss >> test >> unit;
            std::cout << s << " = " << test * getScaleFactorEnergy(unit) << "\n";
        } else {
            std::cout << s << " not found\n";
        }
    }

    // Test time units
    vec = {"time_a", "time_b", "time_c"};
    for (const auto &s: vec) {
        if (inp.getValue(head+s+tail, ss)) {
            ss >> test >> unit;
            std::cout << s << " = " << test * getScaleFactorTime(unit) << "\n";
        } else {
            std::cout << s << " not found\n";
        }
    }
}

int main() {
    InputFile inp("INPUT");
    testValue(inp);
    testUnits(inp);
}