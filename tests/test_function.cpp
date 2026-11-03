/***************************************************************************
 *            test_function.cpp
 *
 *  Copyright  2009-20  Pieter Collins
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
#include "function/function.hpp"
#include "function/formula.hpp"
#include "function/symbolic_function.hpp"
#include "function/taylor_model.hpp"
#include "algebra/algebra.hpp"

#include "numeric/numeric.hpp"

#include "ariadne-testing.hpp"

using namespace std;
using namespace Ariadne;

class TestFunction : public Ariadne::Testing::ContextFixture
{
  public:
    using ContextFixture::ContextFixture;
    Void test_concept();
    Void test_scalar_univariate_function();
    Void test_vector_univariate_function();
    Void test_scalar_function();
    Void test_vector_function();
    Void test_conversions();
    Void test_differentiation();
};



Void TestFunction::test_concept()
{

    EffectiveScalarMultivariateFunction sf1(3);
    EffectiveScalarMultivariateFunction sf2(3);
    EffectiveScalarMultivariateFunction sf3(3);
    //Vector<EffectiveScalarMultivariateFunction> ve=join(3,*e1,*e2,*e3);

    EffectiveVectorMultivariateFunction vf=join(sf1,sf2);

    EffectiveScalarMultivariateFunction g(2);
    EffectiveScalarMultivariateFunction h=compose(g,vf);

    EffectiveVectorMultivariateFunction jf=join(sf1,sf2);
    //EffectiveScalarMultivariateFunction cf=combine(sf1,sf2);

    //MultivariatePolynomial<Real> p;
    //EffectiveScalarMultivariateFunction pf(p);

    //Vector<FloatDPApproximation> b; Matrix<FloatDPApproximation> A;
    //VectorAffineFunction aff(A,b);

}

Void TestFunction::test_scalar_univariate_function()
{
    [[maybe_unused]] FloatDPApproximation p({3.0_approx});

    [[maybe_unused]] EffectiveScalarUnivariateFunction o(EffectiveScalarUnivariateFunction::constant(1));
    [[maybe_unused]] EffectiveScalarUnivariateFunction x(EffectiveScalarUnivariateFunction::coordinate());

    ARIADNE_TESTING_EXPECT_EQUAL(o(p), 1.0_approx);
    ARIADNE_TESTING_EXPECT_EQUAL(x(p), 3.0_approx);

    [[maybe_unused]] EffectiveScalarUnivariateFunction f(o+x*x);
    ARIADNE_TESTING_EXPECT_EQUAL(f(p), 10.0_approx);

    ARIADNE_TESTING_INSPECT("cos(f)", cos(f));

    EffectiveScalarUnivariateFunction df;
    FloatDPApproximation dfp(dp);
    ARIADNE_TESTING_NOTE("No ScalarUnivariateFunction<P>::derivative() method.");
    // df=f.derivative();
    df=derivative(f);
    df=f.derivative(IndexZero());
    df=derivative(f,IndexZero());
    dfp=f.slope(p);
    dfp=slope(f,p);
    ARIADNE_TESTING_INSPECT("df", df);
    ARIADNE_TESTING_EXPECT_EQUAL(df(p), 6.0_approx);

}

Void TestFunction::test_vector_univariate_function()
{
    Vector<Real> c({7,11});

    [[maybe_unused]] EffectiveVectorUnivariateFunction f(EffectiveVectorUnivariateFunction::constant(c));

    FloatDPApproximation p=2.0_approx;

    ARIADNE_TESTING_EXPECT_EQUAL(f(p)[0], c[0]);
    ARIADNE_TESTING_EXPECT_EQUAL(f[0](p), c[0]);
    ARIADNE_TESTING_EXPECT_EQUAL(f(p)[1], c[1]);
    ARIADNE_TESTING_EXPECT_EQUAL(f[1](p), c[1]);

    EffectiveVectorUnivariateFunction df;
    Vector<FloatDPApproximation> dfp;
    derivative(f);
    ARIADNE_TESTING_NOTE("Executed derivative(f)");
    ARIADNE_TESTING_EXPECT_NO_THROW(df=derivative(f));
    df=f.derivative(IndexZero());
    df=derivative(f,IndexZero());
    dfp=f.tangent(p);
    dfp=tangent(f,p);

}

Void TestFunction::test_scalar_function()
{
    [[maybe_unused]] EffectiveScalarMultivariateFunction o(EffectiveScalarMultivariateFunction::constant(3,1));
    [[maybe_unused]] EffectiveScalarMultivariateFunction x(EffectiveScalarMultivariateFunction::coordinate(3,0));
    [[maybe_unused]] EffectiveScalarMultivariateFunction y(EffectiveScalarMultivariateFunction::coordinate(3,1));


    [[maybe_unused]] EffectiveScalarMultivariateFunction f(o+x*y);
    [[maybe_unused]] Vector<FloatDPApproximation> p({2.0_approx,3.0_approx,5.0_approx});
    ARIADNE_TESTING_EXPECT_EQUAL(f(p), 7.0_approx);

    EffectiveFormula formula_x=EffectiveFormula::coordinate(0);
    EffectiveFormula scalar_formula=EffectiveFormula::scalar(
        BinaryElementaryOperator(Mul::code()),EffectiveNumber(2),formula_x);
    EffectiveScalarMultivariateFunction scalar_formula_function(EuclideanDomain(1),scalar_formula);
    Vector<FloatDPApproximation> scalar_formula_argument={3.0_approx};
    ARIADNE_TESTING_EXPECT_EQUAL(scalar_formula_function(scalar_formula_argument), 6.0_approx);

    ARIADNE_TESTING_INSPECT("cos(f)", cos(f));

    EffectiveScalarMultivariateFunction df=f.derivative(1);
    ARIADNE_TESTING_INSPECT("df", df);
    ARIADNE_TESTING_EXPECT_EQUAL(df(p), 2.0_approx);

    Covector<FloatDPApproximation> dfp;
    SizeType k=0;
    df=f.derivative(k);
    df=derivative(f,k);
    dfp=f.gradient(p);
    dfp=gradient(f,p);


}

Void TestFunction::test_vector_function()
{
    [[maybe_unused]] EffectiveVectorMultivariateFunction id(EffectiveVectorMultivariateFunction::identity(3));

    // Regression tests for element
    ARIADNE_TESTING_INSPECT("id[0]", id[0]);
    const EffectiveVectorMultivariateFunction& id_const_ref=id;
    ARIADNE_TESTING_INSPECT("id_const_ref[0]", id_const_ref[0]);
    EffectiveVectorMultivariateFunction& id_ref=id;
    ARIADNE_TESTING_INSPECT("id_ref[0]", id_ref[0]);
    const EffectiveScalarMultivariateFunction x0sqr = sqr(id[0]);
    ARIADNE_TESTING_INSPECT("x0sqr", x0sqr);

    Vector<FloatDPApproximation> v={2.0_approx,3.0_approx,5.0_approx};

    [[maybe_unused]] EffectiveVectorMultivariateFunction f(id);
    ARIADNE_TESTING_EXPECT_EQUAL(f(v), v);
    ARIADNE_TESTING_EXPECT_EQUAL(f[0](v), v[0]);
    auto f0=f[0];
    ARIADNE_TESTING_INSPECT("f0", f0);
    ARIADNE_TESTING_INSPECT("f[0]", f[0]);
    ARIADNE_TESTING_EXPECT_NO_THROW(f[0]);
    ARIADNE_TESTING_EXPECT_NO_THROW(f[0]=x0sqr);
    ARIADNE_TESTING_EXPECT_NO_THROW(f[0]=f[1]);
    ARIADNE_TESTING_EXPECT_EQUAL(f[0](v), v[1]);
    ARIADNE_TESTING_EXPECT_EQUAL(f(v)[0], v[1]);

    // Regression test for function - vector product
    EffectiveScalarMultivariateFunction x1=EffectiveScalarMultivariateFunction::coordinate(3,1);
//    EffectiveVector e0={1,0};
    Vector<EffectiveNumber> e0={1,0};
    ARIADNE_TESTING_EXPECT_EQUAL((x1*e0)(v)[0], v[1]);
    ARIADNE_TESTING_EXPECT_EQUAL((x1*e0)(v)[1], 0);

    // Regression test for validated scalar-vector join preserving every
    // component of the vector tail.
    ValidatedScalarMultivariateFunction vx=
        ValidatedScalarMultivariateFunction::coordinate(2,0);
    ValidatedScalarMultivariateFunction vy=
        ValidatedScalarMultivariateFunction::coordinate(2,1);
    ValidatedVectorMultivariateFunction tail=join(vx+vy,vx-vy);
    ValidatedVectorMultivariateFunction joined=join(vx,tail);
    Vector<FloatDPBounds> point({
        FloatDPBounds(2,dp),FloatDPBounds(3,dp)
    });
    Vector<FloatDPBounds> joined_value=joined(point);
    ARIADNE_TESTING_EXPECT_EQUAL(joined.result_size(), 3u);
    ARIADNE_TESTING_EXPECT_EQUAL(joined_value[0], FloatDPBounds(2,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(joined_value[1], FloatDPBounds(5,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(joined_value[2], FloatDPBounds(-1,dp));
}





Void TestFunction::test_conversions()
{

}

Void TestFunction::test_differentiation()
{
    EuclideanDomain dom(2);
    DegreeType deg(3);

    EffectiveScalarMultivariateFunction z=EffectiveScalarMultivariateFunction::constant(dom,0);
    EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(dom,1);
    EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(dom,0);
    EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(dom,1);

    [[maybe_unused]] EffectiveScalarMultivariateFunction f(3*x-2*y+1);
    [[maybe_unused]] EffectiveScalarMultivariateFunction df0(EffectiveScalarMultivariateFunction::constant(dom,3));
    [[maybe_unused]] EffectiveScalarMultivariateFunction df1(EffectiveScalarMultivariateFunction::constant(dom,-2));

    [[maybe_unused]] Vector<FloatDPBounds> v({2.25_exact,1.375_exact});
    [[maybe_unused]] Vector<FloatDPApproximation> va(v);

    ARIADNE_TESTING_EXPECT_NO_THROW(f.differential(v,deg));
    ARIADNE_TESTING_EXPECT_NO_THROW(f.differential(va,deg));

    ARIADNE_TESTING_EXPECT_EQUAL(f.differential(v,deg).value(), f.evaluate(v));
    ARIADNE_TESTING_EXPECT_EQUAL(f.differential(v,deg).gradient()[0], df0.evaluate(v));
    ARIADNE_TESTING_EXPECT_EQUAL(f.differential(v,deg).gradient()[1], df1.evaluate(v));
    ARIADNE_TESTING_EXPECT_EQUAL(differential(f,v,deg), f.differential(v,deg));
    ARIADNE_TESTING_EXPECT_EQUAL(gradient(f,v), f.differential(v,deg).gradient());
    ARIADNE_TESTING_EXPECT_EQUAL(f.gradient(v), gradient(f,v));

    ARIADNE_TESTING_EXPECT_EQUAL(f.derivative(1).evaluate(v), df1.evaluate(v));
    ARIADNE_TESTING_EXPECT_EQUAL(derivative(f,1).evaluate(v), df1.evaluate(v));

    ARIADNE_TESTING_EXPECT_EQUAL(x.derivative(0).evaluate(v), 1.0_exact);
    ARIADNE_TESTING_EXPECT_EQUAL(x.derivative(1).evaluate(v), 0.0_exact);

}


int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_function");
    ARIADNE_TESTING_ADD_TEST(TestFunction, test_scalar_univariate_function);
    ARIADNE_TESTING_ADD_TEST(TestFunction, test_vector_univariate_function);
    ARIADNE_TESTING_ADD_TEST(TestFunction, test_scalar_function);
    ARIADNE_TESTING_ADD_TEST(TestFunction, test_vector_function);
    ARIADNE_TESTING_ADD_TEST(TestFunction, test_conversions);
    ARIADNE_TESTING_ADD_TEST(TestFunction, test_differentiation);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
