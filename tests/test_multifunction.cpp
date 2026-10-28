/***************************************************************************
 *            test_multifunction.cpp
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
#include "algebra/sweeper.hpp"
#include "function/multifunction.hpp"
#include "function/taylor_multifunction.hpp"
#include "numeric/numeric.hpp"


#include "macro/testing.hpp"

using namespace std;
using namespace Ariadne;

// Explicit instantiation must not require operations absent from interval factories.
namespace Ariadne {
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedIntervalTaylorModelDP>,RealVector>;
template class FunctionModelCreator<ScaledFunctionPatchFactory<ValidatedIntervalTaylorModelMP>,RealVector>;
}

template<class Creator> concept CanCreateZero = requires(Creator& creator) { creator.create_zero(); };
static_assert(CanCreateZero<ScaledFunctionPatchCreator<ValidatedTaylorModelDP>>);
static_assert(!CanCreateZero<ScaledFunctionPatchCreator<ValidatedIntervalTaylorModelDP>>);
static_assert(!CanCreateZero<ScaledFunctionPatchCreator<ValidatedIntervalTaylorModelMP>>);

template<class Model> concept CanMeasureGenericDistance = requires(Model const& f, ValidatedScalarMultivariateFunction const& g) { distance(f,g); };
static_assert(CanMeasureGenericDistance<ScaledFunctionPatch<ValidatedTaylorModelDP>>);
static_assert(!CanMeasureGenericDistance<ScaledFunctionPatch<ValidatedIntervalTaylorModelDP>>);
static_assert(!CanMeasureGenericDistance<ScaledFunctionPatch<ValidatedIntervalTaylorModelMP>>);

class TestMultifunction
{
    using ScalarIntervalFunctionModel = ValidatedIntervalTaylorFunctionModel<FloatMP>;
    using VectorIntervalFunctionModel = ValidatedVectorIntervalTaylorFunctionModel<FloatMP>;
    using UpperInterval = Interval<FloatMPUpperBound>;
  public:
    Void test();
  private:
    Void test_taylor_concept();
    Void test_taylor_evaluate();
    Void test_function_set();
    Void test_inclusion_solutions();
};

Void TestMultifunction::test()
{
    Dyadic::set_default_writer(DecimalWriter());
    ARIADNE_TEST_CALL(test_taylor_evaluate());
    ARIADNE_TEST_CALL(test_function_set());
    ARIADNE_TEST_CALL(test_inclusion_solutions());
}

Void TestMultifunction::test_taylor_concept()
{

    MultiplePrecision pr(128);
    BoxDomainType dom({{-1,+1},{-1,+1}});
    ThresholdSweeper<FloatMP> swp(pr,1e-10);
    UpperInterval c(pr);
    ScalarIntervalFunctionModel f(dom,swp);
    Bounds<FloatMP> x(pr);
    ValidatedNumber y;

    +f; -f;
    f+f; f-f; f*f; f/f;
    f+c; f-c; f*c; f/c;
    c+f; c-f; c*f; c/f;
    f+x; f-x; f*x; f/x;
    x+f; x-f; x*f; x/f;
    f+y; f-y; f*y; f/y;
    y+f; y-f; y*f; y/f;

    nul(f); pos(f); neg(f); sqr(f); rec(f); pow(f,0u); pow(f,0);
    sqrt(f); exp(f); log(f); sin(f); cos(f); tan(f); atan(f);
}

inline Interval<FloatMPUpperBound> create_zero(Interval<FloatMPUpperBound> const& ivl) {
    std::cerr<<"create_zero(ivl)\n";
    std::cerr<<"ivl="<<ivl<<", nul(ivl.lower_bound())="<<nul(ivl.lower_bound())<<"\n";
    return Interval<FloatMPUpperBound>(nul(ivl.lower_bound()),nul(ivl.upper_bound()));
}

Void TestMultifunction::test_taylor_evaluate()
{
    MultiplePrecision pr(128);
    BoxDomainType dom({{0,1},{1,3}});
    ThresholdSweeper<FloatMP> swp(pr,1e-10);
    ScalarIntervalFunctionModel x0=ScalarIntervalFunctionModel::coordinate(dom,0,swp);
    ScalarIntervalFunctionModel x1=ScalarIntervalFunctionModel::coordinate(dom,1,swp);
    UpperInterval one(1,1,pr);
    UpperInterval c01(3,5,pr);
    UpperInterval c(7,11,pr);
    ScalarIntervalFunctionModel xc=ScalarIntervalFunctionModel::constant(dom,c,swp);

    Vector<ValidatedNumber> u({0.75_x,2.5_x});
    //Vector<ValidatedNumber> u({0.75_x,-0.5_x});
    Vector<FloatMPBounds> v(u,pr);
    Vector<FloatDPBounds> w(u,dp);

    ARIADNE_TEST_SAME(xc(v),c);
    ARIADNE_TEST_SAME(x0(v),one*v[0]);
    ARIADNE_TEST_SAME(x1(v),one*v[1]);
    ARIADNE_TEST_CONSTRUCT(ScalarIntervalFunctionModel,sf,(c01*x0*x1+c));
    ARIADNE_TEST_SAME(sf(v),c01*v[0]*v[1]+c);
    //ARIADNE_TEST_SAME(sf(w),c01*w[0]*w[1]+c);
    ARIADNE_TEST_PRINT(sf(w));

    ARIADNE_TEST_CONSTRUCT(VectorIntervalFunctionModel,vf,({c01*x0*x1+c,x0+x1}));
    ARIADNE_TEST_SAME(vf(v),Box<UpperInterval>({c01*v[0]*v[1]+c,one*v[0]+v[1]}));
}

Void TestMultifunction::test_function_set()
{
    using Ivl = Interval<FloatMPUpperBound>;

    MultiplePrecision pr(128);
    BoxDomainType dom({{-1,+1},{-1,+1}});
//    BoxDomainType dom({{0,1},{1,3}});
    ThresholdSweeper<FloatMP> swp(pr,1e-10);
    ScalarIntervalFunctionModel x0=ScalarIntervalFunctionModel::coordinate(dom,0,swp);
    ScalarIntervalFunctionModel x1=ScalarIntervalFunctionModel::coordinate(dom,1,swp);

    ScalarIntervalFunctionModel sf=Ivl(3,5,pr)*x0*x1+Ivl(7,11,pr);

    ValidatedIntervalTaylorFunctionSet<FloatMP> fset(sf);
    FunctionSet<ValidatedTag,Real(RealVector),CompactSet> gset(fset);
    CompactSet<ValidatedTag,Real(RealVector)> gfset(fset);
    const Vector<ValidatedNumber> origin{0,0};
    const CompactSet<ValidatedTag,Real(RealVector)>::Interface& scalar_set_interface=fset;
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(scalar_set_interface._call(origin).bounding_box()),ExactIntervalType(7,11)));
    const MultifunctionInterface<ValidatedTag,Real(RealVector),CompactSet>& scalar_interface=sf;
    ScaledFunctionPatchFactory<ValidatedIntervalTaylorModelMP> interval_factory(swp);
    ARIADNE_TEST_EQUALS(interval_factory.properties().precision(),pr);
    const auto scalar_image=gset(origin);
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(scalar_image.bounding_box()),ExactIntervalType(7,11)));
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(gfset(origin).bounding_box()),ExactIntervalType(7,11)));
    Multifunction<ValidatedTag,Real(RealVector),CompactSet> multifunction(sf);
    ARIADNE_TEST_ASSERT(same(cast_exact_interval(multifunction(origin).bounding_box()),ExactIntervalType(7,11)));
    VectorIntervalFunctionModel vf{sf,x0+x1};
    ValidatedVectorIntervalTaylorFunctionSet<FloatMP> vfset(vf);
    CompactSet<ValidatedTag,RealVector(RealVector)> vector_set(vfset);
    const CompactSet<ValidatedTag,RealVector(RealVector)>::Interface& vector_set_interface=vfset;
    ARIADNE_TEST_ASSERT(same(vector_set_interface._call(origin).bounding_box(),vector_set(origin).bounding_box()));
    const MultifunctionInterface<ValidatedTag,RealVector(RealVector),CompactSet>& vector_interface=vf;
    const auto vector_image=vector_set(origin).bounding_box();
    ARIADNE_TEST_EQUALS(vector_image.dimension(),2u);
    ARIADNE_TEST_ASSERT(vector_image[0].lower_bound().raw()<=FloatDP(7,dp));
    ARIADNE_TEST_ASSERT(vector_image[0].upper_bound().raw()>=FloatDP(11,dp));
    ARIADNE_TEST_EQUALS(vector_image[1].lower_bound().raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(vector_image[1].upper_bound().raw(),FloatDP(0,dp));
    std::ostringstream text;text<<gset<<gfset<<multifunction<<vector_set;
    scalar_interface._write(text);vector_interface._write(text);
    scalar_set_interface._write(text);vector_set_interface._write(text);
    ARIADNE_TEST_ASSERT(text.str().find("FunctionSet(")!=String::npos);


}

Void TestMultifunction::test_inclusion_solutions()
{
    using Ivl = Interval<FloatMPUpperBound>;
    using P=ValidatedTag;

    MultiplePrecision pr(128);
    BoxDomainType dom({{-1,+1},{-1,+1},{0,+1}});
//    BoxDomainType dom({{0,1},{1,3}});
    ThresholdSweeper<FloatMP> swp(pr,1e-10);
    ScalarIntervalFunctionModel x0=ScalarIntervalFunctionModel::coordinate(dom,0,swp);
    ScalarIntervalFunctionModel x1=ScalarIntervalFunctionModel::coordinate(dom,1,swp);
    ScalarIntervalFunctionModel t=ScalarIntervalFunctionModel::coordinate(dom,1,swp);

    VectorIntervalFunctionModel vf={Ivl(3,5,pr)*x0*x1*+Ivl(7,11,pr),x1*cos(t)};

    Vector<ValidatedNumber> w0({.375_x,0.625_x});
    ValidatedNumber s(0.25_x);

    ValidatedIntervalTaylorMultiflow<FloatMP> phi(vf);
    ValidatedIntervalTaylorCurriedTrajectorySet<FloatMP> phi_w0=phi(w0);
    CompactSet<P,RealVector> phi_w0_s=phi_w0(s);

    //    Function<P,CompactSet<P,RealVector(Real)>(RealVector)> const& phif=phi;
    ValidatedCompactMultiflow::Interface const& phif=phi;
    CompactSet<P,RealVector(Real)> phif_w0=phif._call(w0);
    CompactSet<P,RealVector> phif_w0_s=phif_w0(s);

    const auto reach=phi.curry();
    const auto erased_reach=phif._curry();
    ARIADNE_TEST_ASSERT(same(reach(join(w0,s)).bounding_box(),phi_w0_s.bounding_box()));
    ARIADNE_TEST_ASSERT(same(erased_reach(join(w0,s)).bounding_box(),phi_w0_s.bounding_box()));
    std::ostringstream flow_text;flow_text<<phi;
    ARIADNE_TEST_ASSERT(flow_text.str().find("Multiflow(")==0u);

    ARIADNE_TEST_PRINT(phi_w0);
    ARIADNE_TEST_PRINT(phif_w0);
    ARIADNE_TEST_PRINT(phi_w0_s);
    ARIADNE_TEST_PRINT(phif_w0_s);
}

Void test_mixed_precision_comparisons() {
    const MultiplePrecision mp(128u);
    ARIADNE_TEST_ASSERT(definitely(FloatMPUpperBound(1,mp)<FloatDPLowerBound(2,dp)));
    ARIADNE_TEST_ASSERT(is_indeterminate(FloatMPUpperBound(2,mp)<FloatDPLowerBound(2,dp)));
    ARIADNE_TEST_ASSERT(definitely(FloatMPLowerBound(2,mp)>FloatDPUpperBound(1,dp)));
    ARIADNE_TEST_ASSERT(is_indeterminate(FloatMPLowerBound(2,mp)>FloatDPUpperBound(2,dp)));
}

Int main() {
    ARIADNE_TEST_CALL(test_mixed_precision_comparisons());
    TestMultifunction().test();

    return ARIADNE_TEST_FAILURES;
}

