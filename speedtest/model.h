#ifndef TBPLAS_TESTS_MODEL_H
#define TBPLAS_TESTS_MODEL_H

#include <tuple>

#include <base/datatypes.h>
#include <builder/primitive.h>

namespace my::tests {

using tbplas::base::complex_t;
using tbplas::builder::PrimitiveCell;
using model_t = PrimitiveCell<complex_t>;

model_t make_tbg(int idx = 5);

model_t make_graphene_sc(
    double t = -2.7,
    double a = 0.142,
    bool with_onsite = false,
    std::tuple<int, int, int> dim = { 1024, 1024, 1 });

} // namespace my::tests
#endif