/***************************************************************************
 *            test_function_names.cpp
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

#include "function/function.hpp"
#include "function/symbolic_function.hpp"
#include "function/taylor_model.hpp"

#include "ariadne-testing.hpp"
using namespace Ariadne;

class TestFunctionNames : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    // Diagnostics identify a function's paradigm, signature and implementation family.
    template<class T>
    Void check_name(const String& expected) {
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<T>(), expected);
    }
    template<class P>
    Void check_paradigm_names(const String& prefix) {
        using N = Number<P>;
        check_name<ScalarUnivariateFunction<P>>(prefix + "ScalarUnivariateFunction");
        check_name<VectorUnivariateFunction<P>>(prefix + "VectorUnivariateFunction");
        check_name<ScalarMultivariateFunction<P>>(prefix + "ScalarMultivariateFunction");
        check_name<VectorMultivariateFunction<P>>(prefix + "VectorMultivariateFunction");
        check_name<ScalarUnivariateFormulaFunction<N>>(prefix + "ScalarUnivariateFormulaFunction");
        check_name<VectorUnivariateFormulaFunction<N>>(prefix + "VectorUnivariateFormulaFunction");
        check_name<ScalarFormulaFunction<N>>(prefix + "ScalarMultivariateFormulaFunction");
        check_name<VectorFormulaFunction<N>>(prefix + "VectorMultivariateFormulaFunction");
        check_name<ConstantFunction<N, RealScalar>>(prefix + "UnivariateConstantFunction");
        check_name<ConstantFunction<N, RealVector>>(prefix + "MultivariateConstantFunction");
        check_name<CoordinateFunction<P, RealScalar>>(prefix + "UnivariateCoordinateFunction");
        check_name<CoordinateFunction<P, RealVector>>(prefix + "MultivariateCoordinateFunction");
        check_name<UnaryFunction<P, RealScalar>>(prefix + "UnivariateUnaryFunction");
        check_name<UnaryFunction<P, RealVector>>(prefix + "MultivariateUnaryFunction");
        check_name<BinaryFunction<P, RealScalar>>(prefix + "UnivariateBinaryFunction");
        check_name<BinaryFunction<P, RealVector>>(prefix + "MultivariateBinaryFunction");
        check_name<GradedFunction<P, RealScalar>>(prefix + "UnivariateGradedFunction");
        check_name<GradedFunction<P, RealVector>>(prefix + "MultivariateGradedFunction");
        check_name<VectorOfScalarFunction<P, RealScalar>>(prefix + "VectorOfScalarUnivariateFunction");
        check_name<VectorOfScalarFunction<P, RealVector>>(prefix + "VectorOfScalarMultivariateFunction");
        check_name<ComposedFunction<P, RealScalar, RealScalar, RealScalar>>(prefix + "ScalarUnivariateComposedFunction");
        check_name<ComposedFunction<P, RealScalar, RealVector, RealScalar>>(prefix + "ScalarUnivariateComposedFunction");
        check_name<ComposedFunction<P, RealVector, RealScalar, RealScalar>>(prefix + "VectorUnivariateComposedFunction");
        check_name<ComposedFunction<P, RealVector, RealVector, RealScalar>>(prefix + "VectorUnivariateComposedFunction");
        check_name<ComposedFunction<P, RealScalar, RealScalar, RealVector>>(prefix + "ScalarMultivariateComposedFunction");
        check_name<ComposedFunction<P, RealScalar, RealVector, RealVector>>(prefix + "ScalarMultivariateComposedFunction");
        check_name<ComposedFunction<P, RealVector, RealScalar, RealVector>>(prefix + "VectorMultivariateComposedFunction");
        check_name<ComposedFunction<P, RealVector, RealVector, RealVector>>(prefix + "VectorMultivariateComposedFunction");
    }
    void test_all_names() {
        check_name<UnknownError<FloatDP>>("UnknownError<FloatDP>");
        check_name<UnknownError<FloatMP>>("UnknownError<FloatMP>");
        check_paradigm_names<ApproximateTag>("Approximate");
        check_paradigm_names<ValidatedTag>("Validated");
        check_paradigm_names<EffectiveTag>("Effective");
        check_name<JoinedFunction<ValidatedTag, EuclideanDomain, EuclideanDomain, EuclideanDomain>>("ValidatedJoinedFunction");
        check_name<JoinedFunction<EffectiveTag, EuclideanDomain, EuclideanDomain, EuclideanDomain>>("EffectiveJoinedFunction");
        check_name<EmbeddedFunction<EffectiveTag, EuclideanDomain, EuclideanDomain, EuclideanDomain, RealDomain>>("EffectiveEmbeddedFunction");
        check_name<EmbeddedFunction<EffectiveTag, EuclideanDomain, EuclideanDomain, EuclideanDomain, EuclideanDomain>>("EffectiveEmbeddedFunction");
        check_name<LieDerivativeFunction<ValidatedTag>>("ValidatedLieDerivativeFunction");
        check_name<LieDerivativeFunction<EffectiveTag>>("EffectiveLieDerivativeFunction");
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_function_names");
    ARIADNE_TESTING_ADD_TEST(TestFunctionNames, test_all_names);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
