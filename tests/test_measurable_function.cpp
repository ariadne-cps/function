/***************************************************************************
 *            test_measurable_function.cpp
 *
 *  Copyright  2020  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <fenv.h>
#include "function/measurable_function.hpp"
#include "function/formula.hpp"
#include "function/symbolic_function.hpp"
#include "function/taylor_model.hpp"
#include "algebra/algebra.hpp"

#include "numeric/numeric.hpp"


#include "macro/testing.hpp"

using namespace std;
using namespace Ariadne;

template<class X> inline OutputStream& operator<<(OutputStream& os, Sequence<X> const& seq) {
    SequenceWriter writer(4u); return os << writer(seq);
}

class TestMeasurableFunction
{
  public:
    Void test();
  private:
    Void test_concept();
    Void test_scalar_univariate_function();
};

Void TestMeasurableFunction::test()
{
    Dyadic::set_default_writer(DecimalWriter());
    ARIADNE_TEST_CALL(test_scalar_univariate_function());
}

Void TestMeasurableFunction::test_concept()
{
}

Void TestMeasurableFunction::test_scalar_univariate_function()
{
    {
        ARIADNE_TEST_CONSTRUCT(UnionOfIntervals<Dyadic>,ivls,({{1,2},{1.5_dy,2.5_dy},{0,0.5_dy}}));
        ARIADNE_TEST_CONSTRUCT(ValidatedOpenSet<Real>,vops,(wrap_open(ivls)));

        auto x=ValidatedScalarUnivariateFunction::coordinate();
        ValidatedContinuousFunction<Real(Real)> f=x*x;
        ARIADNE_TEST_PRINT(f);
        ARIADNE_TEST_PRINT(vops);
        ValidatedOpenSet<Real> prevops=preimage(f,vops,IntervalDomainType(-1,3),Accuracy(exp2(-4)));
        ARIADNE_TEST_PRINT(prevops);

        IntervalDomainType dom(-1,3);
        Error<FloatDP> err(0.093750_dy,dp);
        FanModel<Real(Real),DP> fmodel(dom,f,err);

        ARIADNE_TEST_PRINT(fmodel);
        ARIADNE_TEST_PRINT(fmodel.preimage(vops));
        ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedLowerMeasurableSet<Real>,prim,fmodel.preimage(vops));
        ARIADNE_TEST_PRINT(prim.measure());
    }
}


Void test_preimage_paths() {
    const auto identity=ValidatedScalarUnivariateFunction::coordinate();
    const IntervalDomainType domain(0,1);
    const Accuracy accuracy(exp2(-4));
    const UnionOfIntervals<Dyadic> covering({{-1,2}}),away({{2,3}}),middle({{0.25_dy,0.75_dy}});
    const ValidatedRegularSet<Real> regular_cover(covering),regular_away(away),regular_middle(middle);
    const ValidatedOpenSet<Real> open_cover(regular_cover.managed_pointer());
    const ValidatedOpenSet<Real> open_away(regular_away.managed_pointer());
    const ValidatedOpenSet<Real> open_middle(regular_middle.managed_pointer());
    // A regular set can reject a disjoint image without bisecting. A plain
    // open set uses bisection; both must enclose the same known preimages.
    const auto covered=preimage_intervals(identity,regular_cover,domain,accuracy);
    const auto disjoint=preimage_intervals(identity,regular_away,domain,accuracy);
    const auto boundary=preimage_intervals(identity,regular_middle,domain,accuracy);
    ARIADNE_TEST_EQUALS(covered.measure(),FloatDP(1,dp));
    ARIADNE_TEST_EQUALS(disjoint.measure(),FloatDP(0,dp));
    const auto unresolved=preimage_intervals(identity,regular_middle,
        IntervalDomainType(0.234375_x,0.265625_x),accuracy);
    ARIADNE_TEST_EQUALS(unresolved.measure(),FloatDP(0,dp));
    ARIADNE_TEST_ASSERT(boundary.measure()<=FloatDP(0.5_x,dp));
    ARIADNE_TEST_ASSERT(boundary.measure()>=FloatDP(0.375_x,dp));
    const auto measure_of=[&](const ValidatedOpenSet<Real>& set) {
        const auto computed_preimage=preimage(identity,set,domain,accuracy);
        const auto intervals=dynamic_pointer_cast<const OpenSetWrapper<UnionOfIntervals<FloatDP>,ValidatedTag,Real>>(computed_preimage.managed_pointer());
        ARIADNE_TEST_ASSERT(intervals!=nullptr);
        return intervals->measure();
    };
    ARIADNE_TEST_EQUALS(measure_of(open_cover),FloatDP(1,dp));
    ARIADNE_TEST_EQUALS(measure_of(open_away),FloatDP(0,dp));
    ARIADNE_TEST_ASSERT(measure_of(open_middle)>=FloatDP(0.375_x,dp));
    ARIADNE_TEST_EQUALS(measure_of(wrap_open(covering)),FloatDP(1,dp));
    ARIADNE_TEST_EQUALS(measure_of(wrap_open(away)),FloatDP(0,dp));
    ARIADNE_TEST_ASSERT(measure_of(wrap_open(middle))>=FloatDP(0.375_x,dp));
    const FanModel<Real(Real),DP> model(domain,identity,FloatDPError(0.125_x,dp));
    ARIADNE_TEST_EQUALS(model.domain(),domain);
    ARIADNE_TEST_EQUALS(model.function()(FloatDPBounds(0.5_x,dp)),FloatDPBounds(0.5_x,dp));
    ARIADNE_TEST_EQUALS(model.measure_error().raw(),FloatDP(0.125_x,dp));
    ARIADNE_TEST_EQUALS(preimage(open_cover,model).measure().get(dp).raw(),FloatDP(0.875_x,dp));
}

Int main() {
    ARIADNE_TEST_CALL(test_preimage_paths());
    TestMeasurableFunction().test();

    return ARIADNE_TEST_FAILURES;
}

