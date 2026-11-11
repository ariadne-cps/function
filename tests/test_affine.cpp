/***************************************************************************
 *            test_affine.cpp
 *
 *  Copyright  2026  Luca Geretti
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Function. If not, see <https://www.gnu.org/licenses/>.
 */

#include <sstream>

#include "function/affine.hpp"
#include "numeric/numeric.hpp"

#include "ariadne-testing.hpp"
using namespace Ariadne;

class TestAffine : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    Void test_affine_evaluation() {
        using A = Affine<FloatDPApproximation>;
        const FloatDPApproximation c(5, dp);
        A f(c, {FloatDPApproximation(2, dp), FloatDPApproximation(-3, dp)});
        Vector<FloatDPApproximation> point({4, 1}, dp);
        ARIADNE_TESTING_EXPECT_SAME(f.evaluate(point), FloatDPApproximation(10, dp));
        ARIADNE_TESTING_EXPECT_SAME(A::constant(2u, c).evaluate(point), c);
        Vector<FloatDPApproximation> empty_point(0u, dp);
        ARIADNE_TESTING_EXPECT_SAME(A::constant(0u, c).evaluate(empty_point), c);
        ARIADNE_TESTING_EXPECT_EQUAL(f.argument_size(), 2u);
        ARIADNE_TESTING_EXPECT_SAME(f.value(), c);
        ARIADNE_TESTING_EXPECT_SAME(f.b(), c);
        ARIADNE_TESTING_EXPECT_SAME(f.gradient(0u), FloatDPApproximation(2, dp));
        ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (f.a()); return snapshot.array(); }()), ([&] { auto snapshot = (f.gradient()); return snapshot.array(); }()), [](const auto& left, const auto& right) { return same(left, right); });
        ARIADNE_TESTING_EXPECT_SAME(f.derivative(1u), FloatDPApproximation(-3, dp));
        ARIADNE_TESTING_EXPECT_SAME(derivative(f, 0u), FloatDPApproximation(2, dp));
        ARIADNE_TESTING_EXPECT_SAME(f[1u], FloatDPApproximation(-3, dp));
        const A& const_f = f;
        ARIADNE_TESTING_EXPECT_SAME(const_f[0u], FloatDPApproximation(2, dp));
        ARIADNE_TESTING_EXPECT_EQUAL(f.coefficient_characteristics(), dp);
        ARIADNE_TESTING_EXPECT_EQUAL(std::get<0>(f.characteristics()), 2u);
        A assigned(2u, dp);
        assigned = c;
        ARIADNE_TESTING_EXPECT_SAME(assigned.evaluate(point), c);
        ARIADNE_TESTING_EXPECT_SAME(assigned.gradient(0u), FloatDPApproximation(0, dp));
        A constant_zero_dim(0u, dp);
        constant_zero_dim = c;
        ARIADNE_TESTING_EXPECT_SAME(constant_zero_dim.evaluate(empty_point), c);
        auto coords = A::coordinates(2u, dp), vars = A::variables(2u, dp);
        ARIADNE_TESTING_EXPECT_SAME(coords[0].evaluate(point), point[0]);
        ARIADNE_TESTING_EXPECT_SAME(vars[1].evaluate(point), point[1]);
        ARIADNE_TESTING_EXPECT_SAME(A::variable(2u, 0u, dp).evaluate(point), point[0]);
        Affine<FloatDP> exact(FloatDP(5, dp), {FloatDP(2, dp), FloatDP(-3, dp)});
        Affine<FloatDPBounds> validated(exact);
        ARIADNE_TESTING_PROVE(definitely(validated.evaluate(Vector<FloatDPBounds>({4, 1}, dp)) == FloatDPBounds(10, dp)));
        std::ostringstream text;
        text << f;
        ARIADNE_TESTING_PROVE(text.str().find("*x0") != String::npos);
    }

    Void test_affine_algebra() {
        using A = Affine<FloatDPApproximation>;
        A f(FloatDPApproximation(5, dp), {FloatDPApproximation(2, dp), FloatDPApproximation(-3, dp)});
        A g(Covector<FloatDPApproximation>({FloatDPApproximation(1, dp), FloatDPApproximation(2, dp)}), FloatDPApproximation(4, dp));
        Vector<FloatDPApproximation> point({4, 1}, dp);
        FloatDPApproximation c(2, dp);
        ARIADNE_TESTING_EXPECT_SAME((-f).evaluate(point), FloatDPApproximation(-10, dp));
        ARIADNE_TESTING_EXPECT_SAME((f + g).evaluate(point), FloatDPApproximation(20, dp));
        ARIADNE_TESTING_EXPECT_SAME((f - g).evaluate(point), FloatDPApproximation(0, dp));
        ARIADNE_TESTING_EXPECT_SAME((f + c).evaluate(point), FloatDPApproximation(12, dp));
        ARIADNE_TESTING_EXPECT_SAME((c + f).evaluate(point), FloatDPApproximation(12, dp));
        ARIADNE_TESTING_EXPECT_SAME((f - c).evaluate(point), FloatDPApproximation(8, dp));
        ARIADNE_TESTING_EXPECT_SAME((c - f).evaluate(point), FloatDPApproximation(-8, dp));
        ARIADNE_TESTING_EXPECT_SAME((f * c).evaluate(point), FloatDPApproximation(20, dp));
        ARIADNE_TESTING_EXPECT_SAME((c * f).evaluate(point), FloatDPApproximation(20, dp));
        ARIADNE_TESTING_EXPECT_SAME((f / c).evaluate(point), FloatDPApproximation(5, dp));
        // Approximate equality may be indeterminate even for a self-comparison.
        // Require non-false evidence and independently verify the representation.
        ARIADNE_TESTING_EXPECT_POSSIBLY(f == f);
        const A duplicate = f;
        ARIADNE_TESTING_EXPECT_SAME(f.value(), duplicate.value());
        const auto coefficients = f.gradient();
        const auto duplicate_coefficients = duplicate.gradient();
        ARIADNE_TESTING_EXPECT_RANGES_BY(coefficients.array(), duplicate_coefficients.array(),
                                         [](const auto& left, const auto& right) { return same(left, right); });
        // Approximate inequality may be undecidable even for distinct representations.
        // Keep the operator evidence, and require a concrete representational difference.
        ARIADNE_TESTING_EXPECT_POSSIBLY(!(f == g));
        ARIADNE_TESTING_PROVE(!same(f.value(), g.value()));
        ARIADNE_TESTING_EXPECT_SAME((f + Integer(2)).evaluate(point), FloatDPApproximation(12, dp));
        ARIADNE_TESTING_EXPECT_SAME((Integer(2) + f).evaluate(point), FloatDPApproximation(12, dp));
        ARIADNE_TESTING_EXPECT_SAME((f - Integer(2)).evaluate(point), FloatDPApproximation(8, dp));
        ARIADNE_TESTING_EXPECT_SAME((Integer(2) - f).evaluate(point), FloatDPApproximation(-8, dp));
        ARIADNE_TESTING_EXPECT_SAME((f * Integer(2)).evaluate(point), FloatDPApproximation(20, dp));
        ARIADNE_TESTING_EXPECT_SAME((Integer(2) * f).evaluate(point), FloatDPApproximation(20, dp));
        ARIADNE_TESTING_EXPECT_SAME((f / Integer(2)).evaluate(point), FloatDPApproximation(5, dp));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_affine");
    ARIADNE_TESTING_ADD_TEST(TestAffine, test_affine_evaluation);
    ARIADNE_TESTING_ADD_TEST(TestAffine, test_affine_algebra);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
