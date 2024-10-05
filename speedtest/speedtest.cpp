#include <iostream>

#include <base/parallel.h>

#include "diag.h"
#include "tbpm.h"

int main(int argc, char* argv[])
{
    tbplas::base::MPIENV_INIT(argc, argv);
    Eigen::initParallel();
    std::string test = "dos";
    if (argc >= 2) {
        test = argv[1];
    }
    if (test == "diag_all") {
        my::tests::test_diag_bands();
        my::tests::test_diag_dos();
        my::tests::test_diag_states();
        my::tests::test_spin();
        my::tests::test_z2();
        my::tests::test_diag_dyn_pol();
        my::tests::test_diag_ac_cond();
    } else if (test == "haydock") {
        my::tests::test_haydock();
    } else if (test == "dos") {
        my::tests::test_dos();
    } else if (test == "ac_cond") {
        my::tests::test_ac_cond();
    } else if (test == "dyn_pol") {
        my::tests::test_dyn_pol();
    } else if (test == "dc_cond") {
        my::tests::test_dc_cond();
    } else if (test == "hall_cond") {
        my::tests::test_hall_cond();
    } else if (test == "qe") {
        my::tests::test_qe();
    } else if (test == "wft") {
        my::tests::test_wft();
    } else if (test == "qe_speed") {
        my::tests::test_qe_speed();
    } else if (test == "wft_speed") {
        my::tests::test_wft_speed();
    } else if (test == "bands") {
        my::tests::test_bands();
    } else {
        std::cout << "Unknown test: " << test << "\n";
    }
    tbplas::base::MPIENV_FINALIZE();
}
