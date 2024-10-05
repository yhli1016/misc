#include "model.h"

#include <base/lattice.h>
#include <base/utils.h>
#include <builder/advanced.h>
#include <builder/base.h>
#include <builder/kdtree.h>
#include <builder/materials.h>
#include <builder/primitive.h>

namespace my::tests {

using tbplas::base::frac2cart;
using tbplas::base::gen_lattice_vectors;
using tbplas::base::Timer;
using tbplas::builder::cell_index_t;
using tbplas::builder::extend_prim_cell;
using tbplas::builder::make_graphene_diamond;
using tbplas::builder::make_hetero_layer;
using tbplas::builder::merge_prim_cell;
using tbplas::builder::orbital_index_t;
using tbplas::builder::PCInterHopping;
using tbplas::builder::Term;

double calc_hop(const Term& term)
{
    constexpr double a0 = 0.1418;
    constexpr double a1 = 0.3349;
    constexpr double r_c = 0.6140;
    constexpr double l_c = 0.0265;
    constexpr double gamma0 = 2.7;
    constexpr double gamma1 = 0.48;
    constexpr double decay = 22.18;
    constexpr double q_pi = decay * a0;
    constexpr double q_sigma = decay * a1;
    double dr = term.distance;
    double n = term.rij(2) / dr;
    double v_pp_pi = -gamma0 * exp(q_pi * (1 - dr / a0));
    double v_pp_sigma = gamma1 * exp(q_sigma * (1 - dr / a1));
    double fc = 1 / (1 + exp((dr - r_c) / l_c));
    double hop = (n * n * v_pp_sigma + (1 - n * n) * v_pp_pi) * fc;
    return hop;
}

model_t make_tbg(int idx)
{
    auto fixed_cell = make_graphene_diamond();
    PrimitiveCell<complex_t> twisted_cell = fixed_cell;

    // Shift on-site energies by -0.78 eV to make Efermi locate at 0 eV
    for (int i = 0; i < 2; ++i) {
        fixed_cell.set_orbital_energy(i, -0.78);
        twisted_cell.set_orbital_energy(i, -0.78);
    }

    // Twist top layer
    double cos_ang = (3 * idx * idx + 3 * idx + 0.5) / (3 * idx * idx + 3 * idx + 1);
    double angle = acos(cos_ang);
    double shift = 0.3349;
    twisted_cell.rotate(angle, Eigen::Vector3d::Zero(), shift);

    // Evaluate hetero-lattice
    double i = idx * 1.0;
    Eigen::Matrix3d hetero_lattice { { i, i + 1, 0 }, { -(i + 1), 2 * i + 1, 0 }, { 0, 0, 1 } };
    hetero_lattice = frac2cart(fixed_cell.get_lattice(), hetero_lattice.transpose()).transpose();

    // Build layers
    PrimitiveCell<complex_t> layer_fixed = make_hetero_layer(fixed_cell, hetero_lattice, 1e-5);
    PrimitiveCell<complex_t> layer_twisted = make_hetero_layer(twisted_cell, hetero_lattice, 1e-5);

    // Build inter-hopping terms
    PCInterHopping<complex_t> inter_hop(0, 1);
    auto neighbors = find_neighbors(layer_fixed, layer_twisted, 1, 1, 0, 0.75);
    for (const auto& n : neighbors) {
        cell_index_t ra, rb, rc;
        orbital_index_t orb_i, orb_j;
        std::tie(ra, rb, rc) = n.rn;
        std::tie(orb_i, orb_j) = n.pair;
        inter_hop.add_hopping(ra, rb, rc, orb_i, orb_j, calc_hop(n));
    }

    // Merge cells
    std::vector<const PrimitiveCell<complex_t>*> prim_cells = { &layer_fixed, &layer_twisted };
    std::vector<const PCInterHopping<complex_t>*> inter_hops = { &inter_hop };
    PrimitiveCell<complex_t> merged_cell = merge_prim_cell(prim_cells, inter_hops);
    return merged_cell;
}

model_t make_graphene_sc(
    double t,
    double a,
    bool with_onsite,
    std::tuple<int, int, int> dim)
{
    Timer timer;
    timer.tic("model");

    double lat = a * sqrt(3);
    Eigen::Matrix3d lat_vec = gen_lattice_vectors(lat, lat, 1.0, 90.0, 90.0, 60.0);
    Eigen::Vector3d origin(0.0, 0.0, 0.0);
    PrimitiveCell<complex_t> prim_cell(2, lat_vec, origin, tbplas::base::NM);
    prim_cell.set_orbital_position(0, 0., 0., 0.0);
    prim_cell.set_orbital_position(1, 1. / 3, 1. / 3, 0.0);
    if (with_onsite) {
        prim_cell.set_orbital_energy(0, -0.1);
        prim_cell.set_orbital_energy(1, 0.1);
    }
    prim_cell.add_hopping(0, 0, 0, 0, 1, t);
    prim_cell.add_hopping(1, 0, 0, 1, 0, t);
    prim_cell.add_hopping(0, 1, 0, 1, 0, t);
    prim_cell = extend_prim_cell(prim_cell, dim);

    timer.toc("model");
    timer.report_total_time();
    return prim_cell;
}

} // namespace my::tests