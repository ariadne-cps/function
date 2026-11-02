/***************************************************************************
 *            test_affine.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "numeric/numeric.hpp"
#include "function/affine.hpp"
#include "ariadne-testing.hpp"
#include <sstream>
using namespace Ariadne;

class TestAffine : public Ariadne::Testing::ContextFixture {
public:
    using ContextFixture::ContextFixture;
Void test_affine_evaluation() {
    using A=Affine<FloatDPApproximation>;
    const FloatDPApproximation c(5,dp);
    A f(c,{FloatDPApproximation(2,dp),FloatDPApproximation(-3,dp)});
    Vector<FloatDPApproximation> point({4,1},dp);
    ARIADNE_TESTING_EXPECT_EQUAL(f.evaluate(point), FloatDPApproximation(10,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(A::constant(2u,c).evaluate(point), c);
    Vector<FloatDPApproximation> empty_point(0u,dp);
    ARIADNE_TESTING_EXPECT_EQUAL(A::constant(0u,c).evaluate(empty_point), c);
    ARIADNE_TESTING_EXPECT_EQUAL(f.argument_size(), 2u);
    ARIADNE_TESTING_EXPECT_EQUAL(f.value(), c);
    ARIADNE_TESTING_EXPECT_EQUAL(f.b(), c);
    ARIADNE_TESTING_EXPECT_EQUAL(f.gradient(0u), FloatDPApproximation(2,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(f.a(), f.gradient());
    ARIADNE_TESTING_EXPECT_EQUAL(f.derivative(1u), FloatDPApproximation(-3,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(derivative(f,0u), FloatDPApproximation(2,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(f[1u], FloatDPApproximation(-3,dp));
    const A& const_f=f;
    ARIADNE_TESTING_EXPECT_EQUAL(const_f[0u], FloatDPApproximation(2,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(f.coefficient_characteristics(), dp);
    ARIADNE_TESTING_EXPECT_EQUAL(std::get<0>(f.characteristics()), 2u);
    A assigned(2u,dp); assigned=c;
    ARIADNE_TESTING_EXPECT_EQUAL(assigned.evaluate(point), c);
    ARIADNE_TESTING_EXPECT_EQUAL(assigned.gradient(0u), FloatDPApproximation(0,dp));
    A constant_zero_dim(0u,dp); constant_zero_dim=c;
    ARIADNE_TESTING_EXPECT_EQUAL(constant_zero_dim.evaluate(empty_point), c);
    auto coords=A::coordinates(2u,dp),vars=A::variables(2u,dp);
    ARIADNE_TESTING_EXPECT_EQUAL(coords[0].evaluate(point), point[0]);
    ARIADNE_TESTING_EXPECT_EQUAL(vars[1].evaluate(point), point[1]);
    ARIADNE_TESTING_EXPECT_EQUAL(A::variable(2u,0u,dp).evaluate(point), point[0]);
    Affine<FloatDP> exact(FloatDP(5,dp),{FloatDP(2,dp),FloatDP(-3,dp)});
    Affine<FloatDPBounds> validated(exact);
    ARIADNE_TESTING_PROVE(definitely(validated.evaluate(Vector<FloatDPBounds>({4,1},dp))==FloatDPBounds(10,dp)));
    std::ostringstream text; text<<f;
    ARIADNE_TESTING_PROVE(text.str().find("*x0")!=String::npos);
}

Void test_affine_algebra() {
    using A=Affine<FloatDPApproximation>;
    A f(FloatDPApproximation(5,dp),{FloatDPApproximation(2,dp),FloatDPApproximation(-3,dp)});
    A g(Covector<FloatDPApproximation>({FloatDPApproximation(1,dp),FloatDPApproximation(2,dp)}),FloatDPApproximation(4,dp));
    Vector<FloatDPApproximation> point({4,1},dp);
    FloatDPApproximation c(2,dp);
    ARIADNE_TESTING_EXPECT_EQUAL((-f).evaluate(point), FloatDPApproximation(-10,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((f+g).evaluate(point), FloatDPApproximation(20,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((f-g).evaluate(point), FloatDPApproximation(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((f+c).evaluate(point), FloatDPApproximation(12,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((c+f).evaluate(point), FloatDPApproximation(12,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((f-c).evaluate(point), FloatDPApproximation(8,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((c-f).evaluate(point), FloatDPApproximation(-8,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((f*c).evaluate(point), FloatDPApproximation(20,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((c*f).evaluate(point), FloatDPApproximation(20,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((f/c).evaluate(point), FloatDPApproximation(5,dp));
    ARIADNE_TESTING_PROVE(f==f);
    ARIADNE_TESTING_PROVE(!(f==g));
    ARIADNE_TESTING_EXPECT_EQUAL((f+Integer(2)).evaluate(point), FloatDPApproximation(12,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((Integer(2)+f).evaluate(point), FloatDPApproximation(12,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((f-Integer(2)).evaluate(point), FloatDPApproximation(8,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((Integer(2)-f).evaluate(point), FloatDPApproximation(-8,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((f*Integer(2)).evaluate(point), FloatDPApproximation(20,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((Integer(2)*f).evaluate(point), FloatDPApproximation(20,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((f/Integer(2)).evaluate(point), FloatDPApproximation(5,dp));
}
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_affine");
    ARIADNE_TESTING_ADD_TEST(TestAffine, test_affine_evaluation);
    ARIADNE_TESTING_ADD_TEST(TestAffine, test_affine_algebra);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
