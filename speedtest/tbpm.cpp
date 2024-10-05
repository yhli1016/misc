#include <base/kpoints.h>
#include <base/utils.h>
#include <builder/advanced.h>
#include <tbpm/solver.h>

#include "model.h"

namespace my::tests {

using tbplas::base::gen_kpath;
using tbplas::base::Timer;
using tbplas::builder::dim_t;
using tbplas::tbpm::TBPMSolver;

#define FAST_ALGO false
#define WITH_ONSITE false
#define QUICK_TEST true

void test_haydock() // Passed
{
    // Model
    double t = -2.7;
    double a = 0.142;
    dim_t dim;
    if (QUICK_TEST) {
        dim = { 1024, 1024, 1 };
    } else {
        dim = { 4096, 4096, 1 };
    }
    model_t model = make_graphene_sc(t, a, WITH_ONSITE, dim);

    // Parameters
    TBPMSolver<model_t> solver(model);
    solver.config.rescale = 9.0;
    solver.config.ldos_recursion_depth = 2000;
    solver.config.dimension = 2;

    // Calculation
    Timer timer;
    timer.tic("ldos_haydock");
    solver.calc_ldos_haydock();
    timer.toc("ldos_haydock");
    timer.report_total_time();
}

void test_dos() // Passed
{
    // Model
    double t = -2.7;
    double a = 0.142;
    dim_t dim;
    if (QUICK_TEST) {
        dim = { 1024, 1024, 1 };
    } else {
        dim = { 4096, 4096, 1 };
    }
    model_t model = make_graphene_sc(t, a, WITH_ONSITE, dim);

    // Parameters
    TBPMSolver<model_t> solver(model);
    solver.config.num_random_samples = 1;
    solver.config.rescale = 9.0;
    solver.config.ldos = false;
    solver.config.ldos_orbital_indices = { 1 };
    solver.config.num_time_steps = 1024;
    solver.config.dimension = 2;
    solver.config.use_fast_algo = FAST_ALGO;

    // Calculation
    Timer timer;
    timer.tic("dos");
    solver.calc_corr_dos();
    timer.toc("dos");
    timer.report_total_time();
}

void test_ac_cond() // Passed
{
    // Model
    // Reference: tbplas/tests/test_units/test_ac.py
    double t = 3.0;
    double a = 0.142;
    dim_t dim;
    if (QUICK_TEST) {
        dim = { 512, 512, 1 };
    } else {
        dim = { 4096, 4096, 1 };
    }
    model_t model = make_graphene_sc(t, a, WITH_ONSITE, dim);

    // Parameters
    TBPMSolver<model_t> solver(model);
    solver.config.rescale = 9.0;
    solver.config.num_random_samples = 1;
    solver.config.num_time_steps = 1024;
    solver.config.num_spin = 2;
    solver.config.dimension = 2;
    solver.config.use_fast_algo = FAST_ALGO;

    // Calculation
    Timer timer;
    timer.tic("ac_cond");
    solver.calc_corr_ac_cond();
    timer.toc("ac_cond");
    timer.report_total_time();
}

void test_dyn_pol() // Passed
{
    // Model
    // Reference: tbplas/tests/test_units/test_dp.py and test_eps.py
    double t = 3.0;
    double a = 0.142;
    dim_t dim;
    if (QUICK_TEST) {
        dim = { 512, 512, 1 };
    } else {
        dim = { 4096, 4096, 1 };
    }
    model_t model = make_graphene_sc(t, a, WITH_ONSITE, dim);

    // Parameters
    TBPMSolver<model_t> solver(model);
    solver.config.rescale = 9.0;
    solver.config.num_random_samples = 1;
    solver.config.num_time_steps = 1024;
    solver.config.num_spin = 1;
    solver.config.dyn_pol_q_points = Eigen::Matrix3Xd(3, 2);
    solver.config.dyn_pol_q_points.col(0) = Eigen::Vector3d(0.86602540, 0.5, 0.0) / a;
    solver.config.dyn_pol_q_points.col(1) = Eigen::Vector3d(4.122280922013927, 2.38, 0.0);
    solver.config.dimension = 2;
    solver.config.use_fast_algo = FAST_ALGO;

    // Calculation
    Timer timer;
    timer.tic("dyn_pol");
    solver.calc_corr_dyn_pol();
    timer.toc("dyn_pol");
    timer.report_total_time();
}

void test_dc_cond() // Passed
{
    // Model
    double t = -2.7;
    double a = 0.142;
    dim_t dim;
    if (QUICK_TEST) {
        dim = { 256, 256, 1 };
    } else {
        dim = { 1024, 1024, 1 };
    }
    model_t model = make_graphene_sc(t, a, WITH_ONSITE, dim);

    // Parameters
    TBPMSolver<model_t> solver(model);
    solver.config.rescale = 9.0;
    solver.config.num_random_samples = 1;
    solver.config.num_time_steps = 1024;
    solver.config.num_spin = 1;
    solver.config.dc_cond_energy_limits = { { -2.5, 2.5 } };
    solver.config.dc_cond_num_time_steps = 1024;
    solver.config.dimension = 2;
    solver.config.use_fast_algo = FAST_ALGO;

    // Calculation
    Timer timer;
    timer.tic("dc_cond");
    solver.calc_corr_dc_cond();
    timer.toc("dc_cond");
    timer.report_total_time();
}

void test_hall_cond() // Passed
{
    // Model
    // Reference: tbplas/tests/test_units/test_hall.py
    double t = 3.0;
    double a = 0.142;
    dim_t dim;
    if (QUICK_TEST) {
        dim = { 256, 256, 1 };
    } else {
        dim = { 1024, 1024, 1 };
    }
    model_t model = make_graphene_sc(t, a, WITH_ONSITE, dim);
    model.apply_magnetic_field(100.0);

    // Parameters
    TBPMSolver<model_t> solver(model);
    solver.config.rescale = 9.0;
    solver.config.num_random_samples = 1;
    solver.config.num_time_steps = 1024;
    solver.config.num_spin = 2;
    solver.config.temperature = 10.0;
    solver.config.dckb_component = 2;
    solver.config.dckb_energies = Eigen::VectorXd::LinSpaced(1000, -1.0, 1.0);
    solver.config.dimension = 2;
    solver.config.dckb_unit = "h";

    // Calculation
    Timer timer;
    timer.tic("hall_cond");
    solver.calc_hall_cond();
    timer.toc("hall_cond");
    timer.report_total_time();
}

void test_qe() // Passed
{
    // Model
    double t = -2.7;
    double a = 0.142;
    model_t model = make_graphene_sc(t, a, WITH_ONSITE, { 17, 17, 1 });
    model.remove_orbitals({ 288 });

    // Parameters
    TBPMSolver<model_t> solver(model);
    solver.config.dimension = 2;
    solver.config.num_random_samples = 1;
    solver.config.rescale = 9.0;
    solver.config.num_time_steps = 1024;
    solver.config.qe_energies = Eigen::Vector3d(-1.0, 0.0, 1.0);
    solver.config.use_fast_algo = FAST_ALGO;

    // Calculation
    Timer timer;
    timer.tic("qe");
    solver.calc_quasi_eigenstates();
    timer.toc("qe");
    timer.report_total_time();
}

void test_wft() // Passed
{
    // Model
    double t = -2.7;
    double a = 0.142;
    model_t model = make_graphene_sc(t, a, WITH_ONSITE, { 17, 17, 1 });

    // Parameters
    TBPMSolver<model_t> solver(model);
    solver.config.dimension = 2;
    solver.config.num_random_samples = 1;
    solver.config.rescale = 9.0;
    solver.config.num_time_steps = 1024;
    solver.config.wft_time_save = { 0, 1, 1024 };
    solver.config.wft_wf0 = Eigen::VectorXcd::Zero(model.get_num_orbitals());
    for (int i = 0; i < model.get_num_orbitals(); ++i) {
        if (i % 2 == 0) {
            solver.config.wft_wf0[i] = 1.0;
        }
    }
    solver.config.use_fast_algo = FAST_ALGO;

    // Calculation
    Timer timer;
    timer.tic("wft");
    solver.calc_wft();
    timer.toc("wft");
    timer.report_total_time();
}

void test_qe_speed()
{
    // Model
    double t = -2.7;
    double a = 0.142;
    model_t model = make_graphene_sc(t, a, WITH_ONSITE, { 4096, 4096, 1 });

    // Parameters
    TBPMSolver<model_t> solver(model);
    solver.config.dimension = 2;
    solver.config.num_random_samples = 1;
    solver.config.rescale = 9.0;
    solver.config.num_time_steps = 1024;
    solver.config.qe_energies = Eigen::Vector3d(-1.0, 0.0, 1.0);
    solver.config.use_fast_algo = FAST_ALGO;
    solver.config.save_data = false;

    // Calculation
    if (!QUICK_TEST) {
        Timer timer;
        timer.tic("qs_speed");
        solver.calc_quasi_eigenstates();
        timer.toc("qs_speed");
        timer.report_total_time();
    }
}

void test_wft_speed()
{
    // Model
    double t = -2.7;
    double a = 0.142;
    model_t model = make_graphene_sc(t, a, WITH_ONSITE, { 4096, 4096, 1 });

    // Parameters
    TBPMSolver<model_t> solver(model);
    solver.config.dimension = 2;
    solver.config.num_random_samples = 1;
    solver.config.rescale = 9.0;
    solver.config.num_time_steps = 1024;
    solver.config.wft_time_save = { 0, 1, 1024 };
    solver.config.wft_wf0 = Eigen::VectorXcd::Zero(model.get_num_orbitals());
    for (int i = 0; i < model.get_num_orbitals(); ++i) {
        if (i % 2 == 0) {
            solver.config.wft_wf0[i] = 1.0;
        }
    }
    solver.config.use_fast_algo = FAST_ALGO;
    solver.config.save_data = false;

    // Calculation
    if (!QUICK_TEST) {
        Timer timer;
        timer.tic("wft_speed");
        solver.calc_wft();
        timer.toc("wft_speed");
        timer.report_total_time();
    }
}

void test_bands()
{
    // Model
    double t = -2.7;
    double a = 0.142;
    model_t model = make_graphene_sc(t, a, WITH_ONSITE, { 3, 3, 1 });

    // Parameters
    TBPMSolver<model_t> solver(model);
    solver.config.dimension = 2;
    solver.config.num_random_samples = 1;
    solver.config.rescale = 9.0;
    solver.config.num_time_steps = 1024;
    solver.config.use_fast_algo = FAST_ALGO;

    // K-points
    Eigen::MatrixX3d k_points { { 0.0, 0.0, 0.0 }, { 2. / 3, 1. / 3, 0.0 }, { 0.5, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
    Eigen::Matrix3Xd k_path;
    Eigen::VectorXi k_idx;
    std::tie(k_path, k_idx) = gen_kpath(k_points.transpose(), { 40, 40, 40 });
    solver.config.k_points = k_path;

    // Calculation
    Timer timer;
    timer.tic("bands");
    solver.calc_corr_bands();
    timer.toc("bands");
    timer.report_total_time();
    std::cout << k_idx << "\n";
}
} // namespace my::tests
