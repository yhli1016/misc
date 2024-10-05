#include <base/kpoints.h>
#include <base/utils.h>
#include <builder/materials.h>
#include <diag/base.h>
#include <diag/lindhard.h>
#include <diag/magnetic.h>
#include <diag/z2.h>

#include "model.h"

namespace my::tests {

using tbplas::base::gen_kmesh;
using tbplas::base::gen_kpath;
using tbplas::base::Timer;
using tbplas::builder::make_bismuth;
using tbplas::builder::make_graphene_diamond;
using tbplas::builder::make_graphene_soc;
using tbplas::diag::DiagSolver;
using tbplas::diag::Lindhard;
using tbplas::diag::SpinTexture;
using tbplas::diag::Z2;

void test_diag_bands()
{
    model_t model = make_graphene_diamond();

    Eigen::MatrixX3d k_points { { 0.0, 0.0, 0.0 }, { 2. / 3, 1. / 3, 0.0 }, { 0.5, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
    Eigen::Matrix3Xd k_path;
    Eigen::VectorXi k_idx;
    std::tie(k_path, k_idx) = gen_kpath(k_points.transpose(), { 1000, 1000, 1000 });
    DiagSolver<model_t> solver(model);
    solver.config.prefix = "graphene";
    solver.config.k_points = k_path;

    Timer timer;
    timer.tic("bands");
    auto data = solver.calc_bands();
    timer.toc("bands");
    if (solver.is_master()) {
        timer.report_total_time();
    }
}

void test_diag_dos()
{
    model_t model = make_graphene_diamond();

    DiagSolver<model_t> solver(model);
    solver.config.prefix = "graphene";
    solver.config.k_points = gen_kmesh({ 1000, 1000, 1 });
    solver.config.e_min = -10.0;
    solver.config.e_max = 10.0;

    Timer timer;
    timer.tic("dos");
    auto data = solver.calc_dos();
    timer.toc("dos");
    if (solver.is_master()) {
        timer.report_total_time();
    }
}

void test_diag_states()
{
    model_t model = make_graphene_diamond();

    DiagSolver<model_t> solver(model);
    solver.config.prefix = "graphene";
    solver.config.k_points = gen_kmesh({ 3, 3, 1 });

    Timer timer;
    timer.tic("states");
    solver.calc_states();
    timer.toc("states");
    if (solver.is_master()) {
        timer.report_total_time();
    }
}

void test_spin()
{
    model_t model_graph = make_graphene_soc();
    //
    SpinTexture<model_t> solver_graph(model_graph);
    solver_graph.config.prefix = "kane_mele";
    solver_graph.config.k_points = 2 * (tbplas::base::gen_kmesh({ 640, 640, 1 }).array() - 0.5);
    solver_graph.config.k_points.row(2).array() = 0.0;
    solver_graph.config.spin_major = false;
    //
    Timer timer;
    timer.tic("spin_kane_mele");
    auto data = solver_graph.calc_spin_texture();
    timer.toc("spin_kane_mele");

    model_t model_bismuth = make_bismuth();
    //
    SpinTexture<model_t> solver_bismuth(model_bismuth);
    solver_bismuth.config.prefix = "bismuth";
    solver_bismuth.config.k_points = 2 * (tbplas::base::gen_kmesh({ 640, 640, 1 }).array() - 0.5);
    solver_bismuth.config.k_points.row(2).array() = 0.0;
    solver_bismuth.config.spin_major = true;
    //
    timer.tic("spin_bismuth");
    data = solver_bismuth.calc_spin_texture();
    timer.toc("spin_bismuth");
    if (solver_bismuth.is_master()) {
        timer.report_total_time();
    }
}

void test_z2()
{
    model_t model = make_graphene_soc();
    Z2<model_t> z2(model);
    z2.config.prefix = "kane_mele";
    z2.config.num_occ = 2;
    z2.config.ka_array = Eigen::VectorXd::LinSpaced(2000, -0.5, 0.5);
    z2.config.kb_array = Eigen::VectorXd::LinSpaced(1000, 0.0, 0.5);
    Timer timer;
    timer.tic("z2");
    auto data = z2.calc_phases();
    timer.toc("z2");
    if (z2.is_master()) {
        timer.report_total_time();
    }
}

void test_diag_dyn_pol()
{
    // Reference: tbplas/tests/test_units/test_lindhard.py
    double t = 3.0;
    double a = 0.142;
    model_t model = make_graphene_diamond(t);
    Lindhard<model_t> lind(model);
    lind.config.prefix = "graphene";
    lind.config.e_min = 0.0;
    lind.config.e_max = 10.0;
    lind.config.e_step = 10.0 / 2048;
    lind.config.dimension = 2;
    lind.config.k_grid_size = { 1024, 1024, 1 };
    lind.config.q_points = Eigen::Matrix3Xd(3, 2);
    lind.config.q_points.col(0) = Eigen::Vector3d(0.86602540, 0.5, 0.0) / a;
    lind.config.q_points.col(1) = Eigen::Vector3d(4.122280922013927, 2.38, 0.0);

    Timer timer;
    timer.tic("dyn_pol");
    auto data = lind.calc_dyn_pol();
    auto epsi = lind.calc_epsilon(data);
    timer.toc("dyn_pol");
    if (lind.is_master()) {
        timer.report_total_time();
    }
}

void test_diag_ac_cond()
{
    // Reference: tbplas/tests/test_units/test_lindhard.py
    double t = 3.0;
    model_t model = make_graphene_diamond(t);
    Lindhard<model_t> lind(model);
    lind.config.prefix = "graphene";
    lind.config.e_min = 0.0;
    lind.config.e_max = 10.0;
    lind.config.e_step = 10.0 / 2048;
    lind.config.dimension = 2;
    lind.config.num_spin = 2;
    lind.config.k_grid_size = { 1024, 1024, 1 };

    // AC Cond
    Timer timer;
    timer.tic("ac_cond");
    auto ac_cond = lind.calc_ac_cond();
    timer.toc("ac_cond");
    if (lind.is_master()) {
        timer.report_total_time();
    }

    // Epsilon_q0
    lind.config.dimension = 3;
    auto eps0 = lind.calc_epsilon_q0(ac_cond);
}

} // namespace my::tests
