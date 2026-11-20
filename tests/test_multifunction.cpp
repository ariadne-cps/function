/***************************************************************************
 *            test_multifunction.cpp
 *
 *  Copyright  2020  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Function.
 *
 *  Ariadne Function is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Function is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Function.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cassert>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

#include <fenv.h>

#include "algebra/sweeper.hpp"
#include "function/multifunction.hpp"
#include "function/taylor_multifunction.hpp"
#include "numeric/numeric.hpp"

#include "ariadne-testing.hpp"

using namespace std;
using namespace Ariadne;

// Explicit instantiation must not require operations absent from interval factories.
namespace Ariadne {
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedIntervalTaylorModelDP>, RealVector>;
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedIntervalTaylorModelMP>, RealVector>;
} // namespace Ariadne

template<class Creator>
concept CanCreateZero = requires(Creator& creator) { creator.create_zero(); };
static_assert(CanCreateZero<ScaledFunctionPatchCreator<ValidatedTaylorModelDP>>);
static_assert(!CanCreateZero<ScaledFunctionPatchCreator<ValidatedIntervalTaylorModelDP>>);
static_assert(!CanCreateZero<ScaledFunctionPatchCreator<ValidatedIntervalTaylorModelMP>>);

template<class Model>
concept CanMeasureGenericDistance = requires(const Model& f, const ValidatedScalarMultivariateFunction& g) { distance(f, g); };
static_assert(CanMeasureGenericDistance<ScaledFunctionPatch<ValidatedTaylorModelDP>>);
static_assert(!CanMeasureGenericDistance<ScaledFunctionPatch<ValidatedIntervalTaylorModelDP>>);
static_assert(!CanMeasureGenericDistance<ScaledFunctionPatch<ValidatedIntervalTaylorModelMP>>);

class TestMultifunction : public Ariadne::Testing::ContextFixture {
    using ScalarIntervalFunctionModel = ValidatedIntervalTaylorFunctionModel<FloatMP>;
    using VectorIntervalFunctionModel = ValidatedVectorIntervalTaylorFunctionModel<FloatMP>;
    using UpperInterval = Interval<FloatMPUpperBound>;
    UpperInterval create_zero(const UpperInterval& ivl);

  public:
    using ContextFixture::ContextFixture;
    void set_up() { Dyadic::set_default_writer(DecimalWriter()); }

    Void test_taylor_concept();
    Void test_taylor_evaluate();
    Void test_double_precision_evaluate();
    Void test_function_set();
    Void test_multiflow_representation();
    Void test_multiflow_dimensions();
};

Void TestMultifunction::test_taylor_concept() {
    MultiplePrecision pr(128);
    BoxDomainType dom({
        {-1, +1},
        {-1, +1}
    });
    ThresholdSweeper<FloatMP> swp(pr, 1e-10);
    UpperInterval c(pr);
    ScalarIntervalFunctionModel f(dom, swp);
    Bounds<FloatMP> x(pr);
    ValidatedNumber y;

    +f;
    -f;
    f + f;
    f - f;
    f* f;
    f / f;
    f + c;
    f - c;
    f* c;
    f / c;
    c + f;
    c - f;
    c* f;
    c / f;
    f + x;
    f - x;
    f* x;
    f / x;
    x + f;
    x - f;
    x* f;
    x / f;
    f + y;
    f - y;
    f* y;
    f / y;
    y + f;
    y - f;
    y* f;
    y / f;

    nul(f);
    pos(f);
    neg(f);
    sqr(f);
    rec(f);
    pow(f, 0u);
    pow(f, 0);
    sqrt(f);
    exp(f);
    log(f);
    sin(f);
    cos(f);
    tan(f);
    atan(f);
}

auto TestMultifunction::create_zero(const UpperInterval& ivl) -> UpperInterval {
    ARIADNE_TESTING_NOTE("create_zero(ivl)");
    ARIADNE_TESTING_NOTE("ivl=" << ivl << ", nul(ivl.lower_bound())=" << nul(ivl.lower_bound()));
    return UpperInterval(nul(ivl.lower_bound()), nul(ivl.upper_bound()));
}

Void TestMultifunction::test_taylor_evaluate() {
    MultiplePrecision pr(128);
    BoxDomainType dom({
        {0, 1},
        {1, 3}
    });
    ThresholdSweeper<FloatMP> swp(pr, 1e-10);
    ScalarIntervalFunctionModel x0 = ScalarIntervalFunctionModel::coordinate(dom, 0, swp);
    ScalarIntervalFunctionModel x1 = ScalarIntervalFunctionModel::coordinate(dom, 1, swp);
    UpperInterval one(1, 1, pr);
    UpperInterval c01(3, 5, pr);
    UpperInterval c(7, 11, pr);
    ScalarIntervalFunctionModel xc = ScalarIntervalFunctionModel::constant(dom, c, swp);

    Vector<ValidatedNumber> u({0.75_x, 2.5_x});
    // Vector<ValidatedNumber> u({0.75_x,-0.5_x});
    Vector<FloatMPBounds> v(u, pr);
    Vector<FloatDPBounds> w(u, dp);

    ARIADNE_TESTING_EXPECT_SAME(xc(v), c);
    ARIADNE_TESTING_EXPECT_SAME(x0(v), one * v[0]);
    ARIADNE_TESTING_EXPECT_SAME(x1(v), one * v[1]);
    [[maybe_unused]] ScalarIntervalFunctionModel sf(c01 * x0 * x1 + c);
    ARIADNE_TESTING_EXPECT_SAME(sf(v), c01 * v[0] * v[1] + c);
    // ARIADNE_TESTING_EXPECT_SAME(sf(w), c01*w[0]*w[1]+c);
    ARIADNE_TESTING_INSPECT("sf(w)", sf(w));

    [[maybe_unused]] VectorIntervalFunctionModel vf({c01 * x0 * x1 + c, x0 + x1});
    ARIADNE_TESTING_EXPECT_SAME(vf(v), Box<UpperInterval>({c01 * v[0] * v[1] + c, one * v[0] + v[1]}));
}

Void TestMultifunction::test_double_precision_evaluate() {
    using ScalarDP = ValidatedIntervalTaylorFunctionModel<FloatDP>;
    using VectorDP = ValidatedVectorIntervalTaylorFunctionModel<FloatDP>;
    using IntervalDP = Interval<FloatDPUpperBound>;

    BoxDomainType dom({
        {0, 1},
        {1, 3}
    });
    ThresholdSweeper<FloatDP> swp(dp, 1e-10);
    ScalarDP x0 = ScalarDP::coordinate(dom, 0, swp);
    ScalarDP x1 = ScalarDP::coordinate(dom, 1, swp);
    IntervalDP c01(3, 5, dp);
    IntervalDP c(7, 11, dp);
    IntervalDP one(1, 1, dp);

    Vector<ValidatedNumber> point({0.75_x, 2.5_x});
    Vector<FloatDPBounds> x(point, dp);
    ScalarDP sf(c01 * x0 * x1 + c);
    VectorDP vf({sf, x0 + x1});

    ARIADNE_TESTING_EXPECT_SAME(sf(x), c01 * x[0] * x[1] + c);
    ARIADNE_TESTING_EXPECT_SAME(vf(x), Box<IntervalDP>({c01 * x[0] * x[1] + c, one * x[0] + x[1]}));
}

Void TestMultifunction::test_function_set() {
    using Ivl = Interval<FloatMPUpperBound>;

    MultiplePrecision pr(128);
    BoxDomainType dom({
        {-1, +1},
        {-1, +1}
    });
    //    BoxDomainType dom({{0,1},{1,3}});
    ThresholdSweeper<FloatMP> swp(pr, 1e-10);
    ScalarIntervalFunctionModel x0 = ScalarIntervalFunctionModel::coordinate(dom, 0, swp);
    ScalarIntervalFunctionModel x1 = ScalarIntervalFunctionModel::coordinate(dom, 1, swp);

    ScalarIntervalFunctionModel sf = Ivl(3, 5, pr) * x0 * x1 + Ivl(7, 11, pr);

    ValidatedIntervalTaylorFunctionSet<FloatMP> fset(sf);
    FunctionSet<ValidatedTag, Real(RealVector), CompactSet> gset(fset);
    CompactSet<ValidatedTag, Real(RealVector)> gfset(fset);
    const Vector<ValidatedNumber> origin{0, 0};
    const CompactSet<ValidatedTag, Real(RealVector)>::Interface& scalar_set_interface = fset;
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(scalar_set_interface._call(origin).bounding_box()), ExactIntervalType(7, 11)));
    const MultifunctionInterface<ValidatedTag, Real(RealVector), CompactSet>& scalar_interface = sf;
    ScaledFunctionPatchFactory<ValidatedIntervalTaylorModelMP> interval_factory(swp);
    ARIADNE_TESTING_EXPECT_EQUAL(interval_factory.properties().precision(), pr);
    const auto scalar_image = gset(origin);
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(scalar_image.bounding_box()), ExactIntervalType(7, 11)));
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(gfset(origin).bounding_box()), ExactIntervalType(7, 11)));
    ValidatedCompactScalarMultivariateMultifunction multifunction(sf);
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(multifunction(origin).bounding_box()), ExactIntervalType(7, 11)));
    VectorIntervalFunctionModel vf{sf, x0 + x1};
    ValidatedCompactVectorMultivariateMultifunction vector_multifunction(vf);
    ARIADNE_TESTING_PROVE(same(vector_multifunction(origin).bounding_box(), vf(origin).bounding_box()));
    ValidatedVectorIntervalTaylorFunctionSet<FloatMP> vfset(vf);
    CompactSet<ValidatedTag, RealVector(RealVector)> vector_set(vfset);
    const CompactSet<ValidatedTag, RealVector(RealVector)>::Interface& vector_set_interface = vfset;
    ARIADNE_TESTING_PROVE(same(vector_set_interface._call(origin).bounding_box(), vector_set(origin).bounding_box()));
    const MultifunctionInterface<ValidatedTag, RealVector(RealVector), CompactSet>& vector_interface = vf;
    const auto vector_image = vector_set(origin).bounding_box();
    ARIADNE_TESTING_EXPECT_EQUAL(vector_image.dimension(), 2u);
    ARIADNE_TESTING_PROVE(vector_image[0].lower_bound().raw() <= FloatDP(7, dp));
    ARIADNE_TESTING_PROVE(vector_image[0].upper_bound().raw() >= FloatDP(11, dp));
    ARIADNE_TESTING_EXPECT_EQUAL(vector_image[1].lower_bound().raw(), FloatDP(0, dp));
    ARIADNE_TESTING_EXPECT_EQUAL(vector_image[1].upper_bound().raw(), FloatDP(0, dp));
    std::ostringstream text;
    text << gset << gfset << multifunction << vector_set;
    scalar_interface._write(text);
    vector_interface._write(text);
    scalar_set_interface._write(text);
    vector_set_interface._write(text);
    ARIADNE_TESTING_PROVE(text.str().find("FunctionSet(") != String::npos);
}

Void TestMultifunction::test_multiflow_representation() {
    using Ivl = Interval<FloatMPUpperBound>;
    using P = ValidatedTag;

    MultiplePrecision pr(128);
    BoxDomainType dom({
        {-1, +1},
        {-1, +1},
        { 0, +1}
    });
    //    BoxDomainType dom({{0,1},{1,3}});
    ThresholdSweeper<FloatMP> swp(pr, 1e-10);
    ScalarIntervalFunctionModel x0 = ScalarIntervalFunctionModel::coordinate(dom, 0, swp);
    ScalarIntervalFunctionModel x1 = ScalarIntervalFunctionModel::coordinate(dom, 1, swp);
    // The time coordinate follows the two initial-state coordinates.
    ScalarIntervalFunctionModel t = ScalarIntervalFunctionModel::coordinate(dom, 2, swp);

    VectorIntervalFunctionModel vf = {Ivl(3, 5, pr) * x0 * x1 * +Ivl(7, 11, pr), x1 * cos(t)};

    Vector<ValidatedNumber> w0({.375_x, 0.625_x});
    ValidatedNumber s(0.25_x);

    ValidatedIntervalTaylorMultiflow<FloatMP> phi(vf);
    ValidatedIntervalTaylorCurriedTrajectorySet<FloatMP> phi_w0 = phi(w0);
    CompactSet<P, RealVector> phi_w0_s = phi_w0(s);

    //    Function<P,CompactSet<P,RealVector(Real)>(RealVector)> const& phif=phi;
    const ValidatedCompactMultiflow::Interface& phif = phi;
    CompactSet<P, RealVector(Real)> phif_w0 = phif._call(w0);
    CompactSet<P, RealVector> phif_w0_s = phif_w0(s);

    // A multiflow built from a model must use the full (initial state, time)
    // argument, including the distinct time coordinate.
    const auto direct_image = vf(join(w0, s));
    ARIADNE_TESTING_PROVE(same(direct_image.bounding_box(), phi_w0_s.bounding_box()));
    const ValidatedNumber later_time(0.75_x);
    const auto later_image = phi_w0(later_time).bounding_box();
    ARIADNE_TESTING_PROVE(same(phi_w0_s.bounding_box()[0u], later_image[0u]));
    ARIADNE_TESTING_PROVE(phi_w0_s.bounding_box()[1u].lower_bound().raw() > later_image[1u].upper_bound().raw());

    const auto reach = phi.curry();
    const auto erased_reach = phif._curry();
    ARIADNE_TESTING_PROVE(same(reach(join(w0, s)).bounding_box(), phi_w0_s.bounding_box()));
    ARIADNE_TESTING_PROVE(same(erased_reach(join(w0, s)).bounding_box(), phi_w0_s.bounding_box()));
    std::ostringstream flow_text;
    flow_text << phi;
    ARIADNE_TESTING_PROVE(flow_text.str().find("Multiflow(") == 0u);

    ARIADNE_TESTING_INSPECT("phi_w0", phi_w0);
    ARIADNE_TESTING_INSPECT("phif_w0", phif_w0);
    ARIADNE_TESTING_INSPECT("phi_w0_s", phi_w0_s);
    ARIADNE_TESTING_INSPECT("phif_w0_s", phif_w0_s);
}

Void TestMultifunction::test_multiflow_dimensions() {
    using Ivl = Interval<FloatMPUpperBound>;
    MultiplePrecision pr(128);
    ThresholdSweeper<FloatMP> swp(pr, 1e-10);
    BoxDomainType state_domain({{-1, 1}, {-1, 1}});
    ScalarIntervalFunctionModel x0 = ScalarIntervalFunctionModel::coordinate(state_domain, 0, swp);
    ScalarIntervalFunctionModel x1 = ScalarIntervalFunctionModel::coordinate(state_domain, 1, swp);

    // A two-dimensional state cannot have a flow model with only two arguments.
    VectorIntervalFunctionModel missing_time({x0, x1});
    bool rejected_missing_time = false;
    try {
        [[maybe_unused]] ValidatedIntervalTaylorMultiflow<FloatMP> invalid(missing_time);
    } catch (const std::exception&) {
        rejected_missing_time = true;
    }
    ARIADNE_TESTING_PROVE(rejected_missing_time);

    BoxDomainType flow_domain({{-1, 1}, {-1, 1}, {0, 1}});
    ScalarIntervalFunctionModel state0 = ScalarIntervalFunctionModel::coordinate(flow_domain, 0, swp);
    ScalarIntervalFunctionModel state1 = ScalarIntervalFunctionModel::coordinate(flow_domain, 1, swp);
    VectorIntervalFunctionModel valid_model({state0 + Ivl(0, 0, pr), state1});
    ValidatedIntervalTaylorMultiflow<FloatMP> valid(valid_model);

    // A curried trajectory requires exactly one initial value per state variable.
    bool rejected_wrong_state_size = false;
    try {
        [[maybe_unused]] auto invalid = valid(Vector<ValidatedNumber>({0.25_x}));
    } catch (const std::exception&) {
        rejected_wrong_state_size = true;
    }
    ARIADNE_TESTING_PROVE(rejected_wrong_state_size);
}

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_multifunction");
    // test_taylor_concept deliberately remains compile-time only.
    ARIADNE_TESTING_ADD_TEST(TestMultifunction, test_taylor_evaluate);
    ARIADNE_TESTING_ADD_TEST(TestMultifunction, test_double_precision_evaluate);
    ARIADNE_TESTING_ADD_TEST(TestMultifunction, test_function_set);
    ARIADNE_TESTING_ADD_TEST(TestMultifunction, test_multiflow_representation);
    ARIADNE_TESTING_ADD_TEST(TestMultifunction, test_multiflow_dimensions);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
