/***************************************************************************
 * test_symbolic_function.cpp
 * or (at your option) any later version.
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

#include "function/formula.tpl.hpp"
#include "function/function.hpp"
#include "function/symbolic_function.hpp"

#include "ariadne-testing.hpp"
using namespace Ariadne;
using S = EffectiveScalarMultivariateFunction;
using SymbolicVectorFunctionType = EffectiveVectorMultivariateFunction;
using Y = EffectiveNumber;
using F = EffectiveFormula;
using X = FloatDPApproximation;

class TestSymbolicFunction : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    Void test_formula_function_interfaces() {
        const auto x = F::coordinate(0u), y = F::coordinate(1u);
        const ScalarUnivariateFormulaFunction<Y> scalar(sqr(x) + Y(1));
        const VectorUnivariateFormulaFunction<Y> curve(Vector<F>{x, sqr(x)});
        const ScalarMultivariateFormulaFunction<Y> sum(2u, x + y);
        const VectorMultivariateFormulaFunction<Y> vector(2u, Vector<F>{x + y, x - y});
        const X parameter(3, dp);
        const Vector<X> point({3, 1}, dp);
        ARIADNE_TESTING_PROVE(identical(F(scalar), sqr(x) + Y(1)));
        ARIADNE_TESTING_EXPECT_EQUAL(scalar.argument_size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(scalar.result_size(), 1u);
        ARIADNE_TESTING_EXPECT_SAME(derivative(scalar, IndexZero())(parameter), X(6, dp));
        ARIADNE_TESTING_EXPECT_EQUAL(curve.argument_size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(curve.result_size(), 2u);
        ARIADNE_TESTING_EXPECT_SAME(curve[1u](parameter), X(9, dp));
        ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (derivative(curve,IndexZero())(parameter)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({1,6},dp)); return snapshot.array(); }()), [](const auto& left, const auto& right) { return same(left, right); });
        ARIADNE_TESTING_EXPECT_EQUAL(vector.result_size(), 2u);
        ARIADNE_TESTING_PROVE(identical(vector.formulae(), Vector<F>({x + y, x - y})));
        ARIADNE_TESTING_EXPECT_SAME(vector[1u](point), X(2, dp));
        ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (derivative(vector,1u)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({1,-1},dp)); return snapshot.array(); }()), [](const auto& left, const auto& right) { return same(left, right); });
        std::ostringstream text;
        text << scalar << curve << vector << representation(scalar) << representation(curve) << representation(sum) << representation(vector);
        ARIADNE_TESTING_PROVE(text.str().find("ScalarUnivariateFormulaFunction(") != String::npos);
        ARIADNE_TESTING_PROVE(text.str().find("VectorUnivariateFormulaFunction(") != String::npos);
        ARIADNE_TESTING_PROVE(text.str().find("ScalarMultivariateFormulaFunction(") != String::npos);
        ARIADNE_TESTING_PROVE(text.str().find("VectorMultivariateFormulaFunction(") != String::npos);
    }

    Void test_input_components() {
        const auto x = F::coordinate(0u), y = F::coordinate(1u), u = F::coordinate(2u);
        const SymbolicVectorFunctionType function = EffectiveVectorMultivariateFormulaFunction(3u, Vector<F>{sqr(x) + 2 * u, y - u});
        const Vector<X> point({3, 4, 5}, dp);
        ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (noise_independent_component(function,1u)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({9,4},dp)); return snapshot.array(); }()), [](const auto& left, const auto& right) { return same(left, right); });
        const auto sensitivities = input_derivatives(function, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(sensitivities.size(), 1u);
        ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (sensitivities[0u](point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({2,-1},dp)); return snapshot.array(); }()), [](const auto& left, const auto& right) { return same(left, right); });
        ARIADNE_TESTING_EXPECT_EQUAL(input_derivatives(function, 0u).size(), 0u);
        ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (noise_independent_component(function,0u)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (function(point)); return snapshot.array(); }()), [](const auto& left, const auto& right) { return same(left, right); });
        const SymbolicVectorFunctionType empty = EffectiveVectorMultivariateFormulaFunction(1u, Vector<F>{});
        ARIADNE_TESTING_EXPECT_EQUAL(substitute(Vector<F>{}, CoordinateFormulaPairs{}).size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(simplify(Vector<F>{}).size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(noise_independent_component(empty, 1u)(Vector<X>({0}, dp)).size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(input_derivatives(empty, 1u)[0u](Vector<X>({0}, dp)).size(), 0u);
        const SymbolicVectorFunctionType non_formula{S::coordinate(3u, 0u), S::coordinate(3u, 1u)};
        ARIADNE_TESTING_EXPECT_EXCEPTION(noise_independent_component(non_formula, 1u));
        ARIADNE_TESTING_EXPECT_EXCEPTION(input_derivatives(non_formula, 1u));
    }

    Void test_symbolic_operations_and_contracts() {
        const auto x = S::coordinate(2u, 0u), y = S::coordinate(2u, 1u);
        const Vector<X> point({3, 1}, dp);
        const ConstantFunction<Y, RealVector> constant(EuclideanDomain(2u), Y(7));
        ARIADNE_TESTING_EXPECT_EQUAL(constant.domain().dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(constant.result_size(), 1u);
        ARIADNE_TESTING_PROVE(same(constant.value(), Y(7)));
        const UnaryFunction<EffectiveTag, RealVector> unary(UnaryElementaryOperator(Neg()), x);
        const BinaryFunction<EffectiveTag, RealVector> binary(BinaryElementaryOperator(Add()), x, y);
        const GradedFunction<EffectiveTag, RealVector> graded(GradedElementaryOperator(Pow()), x, 3);
        ARIADNE_TESTING_EXPECT_EQUAL(unary.result_size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(graded.result_size(), 1u);
        ARIADNE_TESTING_EXPECT_SAME(derivative(graded, 0u)(point), X(27, dp));
        std::ostringstream text;
        text << representation(constant) << representation(unary) << representation(binary) << graded << representation(graded);
        ARIADNE_TESTING_PROVE(text.str().find("ConstantFunction(") != String::npos);
        ARIADNE_TESTING_PROVE(text.str().find("UnaryFunction(") != String::npos);
        ARIADNE_TESTING_PROVE(text.str().find("BinaryFunction(") != String::npos);
        ARIADNE_TESTING_PROVE(text.str().find("GradedFunction(") != String::npos);
        const UnaryFunction<EffectiveTag, RealVector> absolute(UnaryElementaryOperator(Abs()), x);
        const BinaryFunction<EffectiveTag, RealVector> maximum(BinaryElementaryOperator(Max()), x, y);
        const BinaryFunction<EffectiveTag, RealVector> minimum(BinaryElementaryOperator(Min()), x, y);
        ARIADNE_TESTING_EXPECT_EXCEPTION(derivative(absolute, 0u));
        ARIADNE_TESTING_EXPECT_EXCEPTION(derivative(maximum, 0u));
        ARIADNE_TESTING_EXPECT_EXCEPTION(derivative(minimum, 0u));
        const SymbolicVectorFunctionType vector{x, y};
        const FunctionElement<EffectiveTag, RealVector> element(vector, 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(element.argument_size(), 2u);
        ARIADNE_TESTING_EXPECT_SAME(element(point), X(1, dp));
        using Element = FunctionElement<EffectiveTag, RealVector>;
        ARIADNE_TESTING_EXPECT_EXCEPTION(Element(vector, 2u));
        using Embedded = EmbeddedFunction<EffectiveTag, EuclideanDomain, EuclideanDomain, EuclideanDomain, RealDomain>;
        const Embedded embedded(1u, S::coordinate(1u, 0u), 1u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(derivative(embedded, 3u));
        using Composed = ComposedFunction<EffectiveTag, Real, RealVector, RealVector>;
        ARIADNE_TESTING_EXPECT_EXCEPTION(Composed(S::coordinate(3u, 0u), vector));
        using Joined = JoinedFunction<EffectiveTag, EuclideanDomain, EuclideanDomain, EuclideanDomain>;
        const Joined joined(vector, SymbolicVectorFunctionType{x + y});
        ARIADNE_TESTING_EXPECT_EQUAL(joined.argument_size(), 2u);
        ARIADNE_TESTING_EXPECT_SAME(joined[0u](point), X(3, dp));
        ARIADNE_TESTING_EXPECT_SAME(joined[2u](point), X(4, dp));
        ARIADNE_TESTING_EXPECT_EXCEPTION(Joined(vector, SymbolicVectorFunctionType::identity(3u)));
        VectorOfScalarFunction<EffectiveTag, RealVector> mutable_vector(Vector<S>{x, y});
        ARIADNE_TESTING_EXPECT_EXCEPTION(mutable_vector[0u] = S::coordinate(3u, 0u));
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_symbolic_function");
    ARIADNE_TESTING_ADD_TEST(TestSymbolicFunction, test_formula_function_interfaces);
    ARIADNE_TESTING_ADD_TEST(TestSymbolicFunction, test_input_components);
    ARIADNE_TESTING_ADD_TEST(TestSymbolicFunction, test_symbolic_operations_and_contracts);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
