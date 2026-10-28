/***************************************************************************
 * This file is part of Ariadne, distributed under the GNU GPL version 3
 * or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/symbolic_function.hpp"
#include "function/formula.tpl.hpp"
#include "macro/testing.hpp"
#include <sstream>
using namespace Ariadne;
using S=EffectiveScalarMultivariateFunction;
using V=EffectiveVectorMultivariateFunction;
using Y=EffectiveNumber;
using F=EffectiveFormula;
using X=FloatDPApproximation;

Void test_formula_function_interfaces() {
    const auto x=F::coordinate(0u),y=F::coordinate(1u);
    const ScalarUnivariateFormulaFunction<Y> scalar(sqr(x)+Y(1));
    const VectorUnivariateFormulaFunction<Y> curve(Vector<F>{x,sqr(x)});
    const ScalarMultivariateFormulaFunction<Y> sum(2u,x+y);
    const VectorMultivariateFormulaFunction<Y> vector(2u,Vector<F>{x+y,x-y});
    const X parameter(3,dp);const Vector<X> point({3,1},dp);
    ARIADNE_TEST_ASSERT(identical(F(scalar),sqr(x)+Y(1)));
    ARIADNE_TEST_EQUALS(scalar.argument_size(),1u);
    ARIADNE_TEST_EQUALS(scalar.result_size(),1u);
    ARIADNE_TEST_EQUALS(derivative(scalar,IndexZero())(parameter),X(6,dp));
    ARIADNE_TEST_EQUALS(curve.argument_size(),1u);
    ARIADNE_TEST_EQUALS(curve.result_size(),2u);
    ARIADNE_TEST_EQUALS(curve[1u](parameter),X(9,dp));
    ARIADNE_TEST_EQUALS(derivative(curve,IndexZero())(parameter),Vector<X>({1,6},dp));
    ARIADNE_TEST_EQUALS(vector.result_size(),2u);
    ARIADNE_TEST_ASSERT(identical(vector.formulae(),Vector<F>({x+y,x-y})));
    ARIADNE_TEST_EQUALS(vector[1u](point),X(2,dp));
    ARIADNE_TEST_EQUALS(derivative(vector,1u)(point),Vector<X>({1,-1},dp));
    std::ostringstream text;
    text<<scalar<<curve<<vector<<representation(scalar)<<representation(curve)<<representation(sum)<<representation(vector);
    ARIADNE_TEST_ASSERT(text.str().find("ScalarUnivariateFormulaFunction(")!=String::npos);
    ARIADNE_TEST_ASSERT(text.str().find("VectorUnivariateFormulaFunction(")!=String::npos);
    ARIADNE_TEST_ASSERT(text.str().find("ScalarMultivariateFormulaFunction(")!=String::npos);
    ARIADNE_TEST_ASSERT(text.str().find("VectorMultivariateFormulaFunction(")!=String::npos);
}

Void test_input_components() {
    const auto x=F::coordinate(0u),y=F::coordinate(1u),u=F::coordinate(2u);
    const V function=EffectiveVectorMultivariateFormulaFunction(3u,Vector<F>{sqr(x)+2*u,y-u});
    const Vector<X> point({3,4,5},dp);
    ARIADNE_TEST_EQUALS(noise_independent_component(function,1u)(point),Vector<X>({9,4},dp));
    const auto sensitivities=input_derivatives(function,1u);
    ARIADNE_TEST_EQUALS(sensitivities.size(),1u);
    ARIADNE_TEST_EQUALS(sensitivities[0u](point),Vector<X>({2,-1},dp));
    ARIADNE_TEST_EQUALS(input_derivatives(function,0u).size(),0u);
    ARIADNE_TEST_EQUALS(noise_independent_component(function,0u)(point),function(point));
    const V empty=EffectiveVectorMultivariateFormulaFunction(1u,Vector<F>{});
    ARIADNE_TEST_EQUALS(noise_independent_component(empty,1u)(Vector<X>({0},dp)).size(),0u);
    ARIADNE_TEST_EQUALS(input_derivatives(empty,1u)[0u](Vector<X>({0},dp)).size(),0u);
    const V non_formula{S::coordinate(3u,0u),S::coordinate(3u,1u)};
    ARIADNE_TEST_FAIL(noise_independent_component(non_formula,1u));
    ARIADNE_TEST_FAIL(input_derivatives(non_formula,1u));
}

Void test_symbolic_operations_and_contracts() {
    const auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    const Vector<X> point({3,1},dp);
    const ConstantFunction<Y,RealVector> constant(EuclideanDomain(2u),Y(7));
    ARIADNE_TEST_EQUALS(constant.domain().dimension(),2u);
    ARIADNE_TEST_EQUALS(constant.result_size(),1u);
    ARIADNE_TEST_ASSERT(same(constant.value(),Y(7)));
    const UnaryFunction<EffectiveTag,RealVector> unary(UnaryElementaryOperator(Neg()),x);
    const BinaryFunction<EffectiveTag,RealVector> binary(BinaryElementaryOperator(Add()),x,y);
    const GradedFunction<EffectiveTag,RealVector> graded(GradedElementaryOperator(Pow()),x,3);
    ARIADNE_TEST_EQUALS(unary.result_size(),1u);
    ARIADNE_TEST_EQUALS(graded.result_size(),1u);
    ARIADNE_TEST_EQUALS(derivative(graded,0u)(point),X(27,dp));
    std::ostringstream text;text<<representation(constant)<<representation(unary)<<representation(binary)<<graded<<representation(graded);
    ARIADNE_TEST_ASSERT(text.str().find("ConstantFunction(")!=String::npos);
    ARIADNE_TEST_ASSERT(text.str().find("UnaryFunction(")!=String::npos);
    ARIADNE_TEST_ASSERT(text.str().find("BinaryFunction(")!=String::npos);
    ARIADNE_TEST_ASSERT(text.str().find("GradedFunction(")!=String::npos);
    const UnaryFunction<EffectiveTag,RealVector> absolute(UnaryElementaryOperator(Abs()),x);
    const BinaryFunction<EffectiveTag,RealVector> maximum(BinaryElementaryOperator(Max()),x,y);
    const BinaryFunction<EffectiveTag,RealVector> minimum(BinaryElementaryOperator(Min()),x,y);
    ARIADNE_TEST_FAIL(derivative(absolute,0u));
    ARIADNE_TEST_FAIL(derivative(maximum,0u));
    ARIADNE_TEST_FAIL(derivative(minimum,0u));
    const V vector{x,y};
    const FunctionElement<EffectiveTag,RealVector> element(vector,1u);
    ARIADNE_TEST_EQUALS(element.argument_size(),2u);
    ARIADNE_TEST_EQUALS(element(point),X(1,dp));
    using Element=FunctionElement<EffectiveTag,RealVector>;
    ARIADNE_TEST_FAIL(Element(vector,2u));
    using Embedded=EmbeddedFunction<EffectiveTag,EuclideanDomain,EuclideanDomain,EuclideanDomain,RealDomain>;
    const Embedded embedded(1u,S::coordinate(1u,0u),1u);
    ARIADNE_TEST_FAIL(derivative(embedded,3u));
    using Composed=ComposedFunction<EffectiveTag,Real,RealVector,RealVector>;
    ARIADNE_TEST_FAIL(Composed(S::coordinate(3u,0u),vector));
    using Joined=JoinedFunction<EffectiveTag,EuclideanDomain,EuclideanDomain,EuclideanDomain>;
    const Joined joined(vector,V{x+y});
    ARIADNE_TEST_EQUALS(joined.argument_size(),2u);
    ARIADNE_TEST_EQUALS(joined[0u](point),X(3,dp));
    ARIADNE_TEST_EQUALS(joined[2u](point),X(4,dp));
    ARIADNE_TEST_FAIL(Joined(vector,V::identity(3u)));
    VectorOfScalarFunction<EffectiveTag,RealVector> mutable_vector(Vector<S>{x,y});
    ARIADNE_TEST_FAIL(mutable_vector[0u]=S::coordinate(3u,0u));
}
Int main() {
    ARIADNE_TEST_CALL(test_formula_function_interfaces());
    ARIADNE_TEST_CALL(test_input_components());
    ARIADNE_TEST_CALL(test_symbolic_operations_and_contracts());
    return ARIADNE_TEST_FAILURES;
}
