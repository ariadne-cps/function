/***************************************************************************
 *            test_function_dispatch.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "ariadne-testing.hpp"
#include <array>
using namespace Ariadne;
using X=FloatDPApproximation;

class TestFunctionDispatch : public Ariadne::Testing::ContextFixture {
public:
    using ContextFixture::ContextFixture;
template<class P> Void test_formula_constructors_and_dispatch() {
    using Y=Number<P>;using C=FunctionConstructors<P>;
    using U=ScalarUnivariateFunction<P>;using W=VectorUnivariateFunction<P>;
    using S=ScalarMultivariateFunction<P>;using V=VectorMultivariateFunction<P>;
    auto t=Formula<Y>::coordinate(0u),y=Formula<Y>::coordinate(1u);
    U scalar(RealDomain(),sqr(t)+Y(1));
    W curve(RealDomain(),Vector<Formula<Y>>{t,t+Y(2)});
    S sum(EuclideanDomain(2u),t+y);
    V identity(EuclideanDomain(2u),Vector<Formula<Y>>{t,y});
    X parameter(3,dp);Vector<X> point({3,1},dp);
    ARIADNE_TESTING_EXPECT_SAME(scalar(parameter), X(10,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (curve(parameter)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({3,5},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_SAME(sum(point), X(4,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (identity(point)); return snapshot.array(); }()), ([&] { auto snapshot = (point); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    auto coordinate=C::coordinate();
    using UA=AlgebraOperations<U,Y>;using SA=AlgebraOperations<S,Y>;
    ARIADNE_TESTING_EXPECT_SAME(UA::apply(BinaryElementaryOperator(Sub()),Y(5),coordinate)(parameter), X(2,dp));
    ARIADNE_TESTING_EXPECT_SAME(UA::apply(GradedElementaryOperator(Pow()),coordinate,3)(parameter), X(27,dp));
    std::array<S,2> first{C::coordinate(2u,0u),C::coordinate(EuclideanDomain(2u),0u)};
    std::array<S,2> second{C::coordinate(2u,1u),C::coordinate(EuclideanDomain(2u),1u)};
    for(const auto& a:first) {
        ARIADNE_TESTING_EXPECT_SAME((-a)(point), X(-3,dp));
        ARIADNE_TESTING_EXPECT_SAME(SA::apply(BinaryElementaryOperator(Sub()),Y(5),a)(point), X(2,dp));
        ARIADNE_TESTING_EXPECT_SAME(SA::apply(BinaryElementaryOperator(Sub()),a,Y(5))(point), X(-2,dp));
        V single{a};
        ARIADNE_TESTING_EXPECT_SAME(single(point)[0], X(3,dp));
        for(const auto& b:second) {
            ARIADNE_TESTING_EXPECT_SAME((a+b)(point), X(4,dp));
            ARIADNE_TESTING_EXPECT_SAME((a-b)(point), X(2,dp));
            ARIADNE_TESTING_EXPECT_SAME((a*b)(point), X(3,dp));
            ARIADNE_TESTING_EXPECT_SAME((a/b)(point), X(3,dp));
        }
    }
    ARIADNE_TESTING_EXPECT_EXCEPTION(SA::apply(BinaryElementaryOperator(Add()),C::coordinate(1u,0u),first[0]));
}

Void test_model_composition() {
    using S=ValidatedScalarMultivariateFunction;using V=ValidatedVectorMultivariateFunction;
    auto t=ValidatedScalarUnivariateFunction::coordinate();
    ValidatedVectorUnivariateFunction curve{t,ValidatedNumber(2)*t};
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    V vector{x+y,x-y};auto scalar=ValidatedNumber(2)*x-y+ValidatedNumber(1);
    BoxDomainType domain{{-1,1},{-1,1}};
    TaylorFunctionFactory factory(ThresholdSweeper<FloatDP>(dp,1e-12));
    S scalar_model=factory.create_constant(domain,ValidatedNumber(2));
    V vector_model=factory.create_identity(domain);
    ARIADNE_TESTING_EXPECT_EXCEPTION(compose(S::coordinate(3u,0u),vector_model));
    Vector<FloatDPBounds> point({0.25_x,0.5_x},dp);
    ARIADNE_TESTING_PROVE(definitely(compose(t+ValidatedNumber(1),scalar_model)(point)==FloatDPBounds(3,dp)));
    ARIADNE_TESTING_EXPECT_EQUAL(compose(curve,scalar_model)(point), Vector<FloatDPBounds>({2,4},dp));
    ARIADNE_TESTING_PROVE(definitely(compose(scalar,vector_model)(point)==FloatDPBounds(1,dp)));
    ARIADNE_TESTING_EXPECT_EQUAL(compose(vector,vector_model)(point), Vector<FloatDPBounds>({0.75_x,-0.25_x},dp));
}

Void test_interval_ranges() {
    using S=ValidatedScalarMultivariateFunction;using V=ValidatedVectorMultivariateFunction;
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    auto scalar=ValidatedNumber(2)*x-y+ValidatedNumber(1);V vector{x+y,x-y};
    Vector<UpperIntervalType> input{UpperIntervalType(-1,1),UpperIntervalType(2,3)};
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(evaluate_range(scalar,input)),ExactIntervalType(-4,1)));
    auto ranges=evaluate_range(vector,input);
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(ranges[0]),ExactIntervalType(1,4)));
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(ranges[1]),ExactIntervalType(-4,-1)));
    auto gradient=gradient_range(scalar,input);auto jacobian=jacobian_range(vector,input);
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(gradient[0]),ExactIntervalType(2,2)));
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(gradient[1]),ExactIntervalType(-1,-1)));
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(jacobian[0][0]),ExactIntervalType(1,1)));
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(jacobian[1][1]),ExactIntervalType(-1,-1)));
    using D=Differential<UpperIntervalType>;
    auto variables=D::identity(1u,input);
    auto derivatives=derivative_range(vector,variables);
    ARIADNE_TESTING_EXPECT_EQUAL(derivatives.size(), 2u);
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(derivatives[0].value()),ExactIntervalType(1,4)));
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(derivatives[1].gradient()[1]),ExactIntervalType(-1,-1)));
    ARIADNE_TESTING_EXPECT_EQUAL(derivative_range(V::zeros(0u,2u),variables).size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(evaluate_range(V::zeros(0u,2u),input).size(), 0u);
    Vector<UpperIntervalType> empty(0u,UpperIntervalType(0,0));
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(evaluate_range(S::constant(0u,5),empty)),ExactIntervalType(5,5)));
    ARIADNE_TESTING_EXPECT_EQUAL(gradient_range(S::constant(0u,5),empty).size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(jacobian_range(V::zeros(0u,0u),empty).row_size(), 0u);
    Vector<D> no_variables(0u,D(0u,1u,UpperIntervalType(0,0)));
    auto constant=derivative_range(V::constant(0u,Vector<ValidatedNumber>{2,3}),no_variables);
    ARIADNE_TESTING_EXPECT_EQUAL(constant.size(), 2u);
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(constant[0].value()),ExactIntervalType(2,2)));
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(constant[1].value()),ExactIntervalType(3,3)));
    ARIADNE_TESTING_EXPECT_EQUAL(constant[0].degree(), 1u);
    Vector<D> unused_parameters(0u,D(3u,2u,UpperIntervalType(0,0)));
    auto lifted=derivative_range(V::constant(0u,Vector<ValidatedNumber>{2,3}),unused_parameters);
    ARIADNE_TESTING_EXPECT_EQUAL(lifted[0].argument_size(), 3u);
    ARIADNE_TESTING_EXPECT_EQUAL(lifted[0].degree(), 2u);
    ARIADNE_TESTING_PROVE(same(cast_exact_interval(lifted[0].gradient()[0]),ExactIntervalType(0,0)));
    Vector<UpperIntervalType> wrong_dimension(1u,UpperIntervalType(0,0));
    Vector<D> wrong_differentials(1u,D(2u,1u,UpperIntervalType(0,0)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(evaluate_range(scalar,wrong_dimension));
    ARIADNE_TESTING_EXPECT_EXCEPTION(evaluate_range(vector,wrong_dimension));
    ARIADNE_TESTING_EXPECT_EXCEPTION(gradient_range(scalar,wrong_dimension));
    ARIADNE_TESTING_EXPECT_EXCEPTION(jacobian_range(vector,wrong_dimension));
    ARIADNE_TESTING_EXPECT_EXCEPTION(derivative_range(vector,wrong_differentials));
}
template<class P> Void test_vector_contracts() {
    using S=ScalarMultivariateFunction<P>;using V=VectorMultivariateFunction<P>;
    using U=ScalarUnivariateFunction<P>;using W=VectorUnivariateFunction<P>;
    auto s1=S::coordinate(1u,0u),s2=S::coordinate(2u,0u);
    auto vector=V::identity(2u),short_result=V::zeros(1u,2u),long_result=V::zeros(3u,2u),wrong_arguments=V::zeros(2u,3u);
    if constexpr (!Same<P,ApproximateTag>) {
        ARIADNE_TESTING_EXPECT_EXCEPTION(vector+short_result);
        ARIADNE_TESTING_EXPECT_EXCEPTION(vector+long_result);
        ARIADNE_TESTING_EXPECT_EXCEPTION(vector+wrong_arguments);
        ARIADNE_TESTING_EXPECT_EXCEPTION(vector-short_result);
        ARIADNE_TESTING_EXPECT_EXCEPTION(vector-long_result);
        ARIADNE_TESTING_EXPECT_EXCEPTION(vector-wrong_arguments);
        ARIADNE_TESTING_EXPECT_EXCEPTION(s1*vector);
        ARIADNE_TESTING_EXPECT_EXCEPTION(vector*s1);
        ARIADNE_TESTING_EXPECT_EXCEPTION(vector/s1);
        auto t=U::coordinate();W curve{t,t},short_curve{t},long_curve{t,t,t};
        ARIADNE_TESTING_EXPECT_EXCEPTION(curve+short_curve);
        ARIADNE_TESTING_EXPECT_EXCEPTION(curve+long_curve);
    }
    ARIADNE_TESTING_EXPECT_EXCEPTION(join(s1,s2));
    ARIADNE_TESTING_EXPECT_EXCEPTION(join(vector,s1));
    ARIADNE_TESTING_EXPECT_EXCEPTION(join(s1,vector));
    ARIADNE_TESTING_EXPECT_EXCEPTION(join(vector,wrong_arguments));
    ARIADNE_TESTING_EXPECT_EXCEPTION(compose(s2,long_result));
}
Void test_effective_contracts_and_derivative_failures() {
    using S=EffectiveScalarMultivariateFunction;using V=EffectiveVectorMultivariateFunction;
    auto vector=V::identity(2u);Vector<EffectiveNumber> short_constant{1};
    ARIADNE_TESTING_EXPECT_EXCEPTION(short_constant+vector);
    ARIADNE_TESTING_EXPECT_EXCEPTION(vector+short_constant);
    ARIADNE_TESTING_EXPECT_EXCEPTION(short_constant-vector);
    ARIADNE_TESTING_EXPECT_EXCEPTION(vector-short_constant);
    auto s3=S::coordinate(3u,0u),s2=S::coordinate(2u,0u);
    auto non_square=V::zeros(2u,3u),wrong_g=V::zeros(2u,3u);
    ARIADNE_TESTING_EXPECT_EXCEPTION(lie_derivative(s3,vector));
    ARIADNE_TESTING_EXPECT_EXCEPTION(lie_derivative(s2,non_square));
    ARIADNE_TESTING_EXPECT_EXCEPTION(lie_derivative(wrong_g,vector));
    ARIADNE_TESTING_EXPECT_EXCEPTION(lie_derivative(vector,non_square));
    auto x=Formula<EffectiveNumber>::coordinate(0u),y=Formula<EffectiveNumber>::coordinate(1u);
    S nonsmooth(EuclideanDomain(2u),max(x,y));
    ARIADNE_TESTING_EXPECT_EXCEPTION(lie_derivative(nonsmooth,vector));
    ARIADNE_TESTING_EXPECT_EXCEPTION(lie_derivative(V{nonsmooth,s2},vector));
}
    void test_formula_constructors_and_dispatch_approximatetag() { test_formula_constructors_and_dispatch<ApproximateTag>(); }
    void test_formula_constructors_and_dispatch_validatedtag() { test_formula_constructors_and_dispatch<ValidatedTag>(); }
    void test_formula_constructors_and_dispatch_effectivetag() { test_formula_constructors_and_dispatch<EffectiveTag>(); }
    void test_vector_contracts_approximatetag() { test_vector_contracts<ApproximateTag>(); }
    void test_vector_contracts_validatedtag() { test_vector_contracts<ValidatedTag>(); }
    void test_vector_contracts_effectivetag() { test_vector_contracts<EffectiveTag>(); }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_function_dispatch");
    ARIADNE_TESTING_ADD_TEST(TestFunctionDispatch, test_formula_constructors_and_dispatch_approximatetag);
    ARIADNE_TESTING_ADD_TEST(TestFunctionDispatch, test_formula_constructors_and_dispatch_validatedtag);
    ARIADNE_TESTING_ADD_TEST(TestFunctionDispatch, test_formula_constructors_and_dispatch_effectivetag);
    ARIADNE_TESTING_ADD_TEST(TestFunctionDispatch, test_vector_contracts_approximatetag);
    ARIADNE_TESTING_ADD_TEST(TestFunctionDispatch, test_vector_contracts_validatedtag);
    ARIADNE_TESTING_ADD_TEST(TestFunctionDispatch, test_vector_contracts_effectivetag);
    ARIADNE_TESTING_ADD_TEST(TestFunctionDispatch, test_effective_contracts_and_derivative_failures);
    ARIADNE_TESTING_ADD_TEST(TestFunctionDispatch, test_model_composition);
    ARIADNE_TESTING_ADD_TEST(TestFunctionDispatch, test_interval_ranges);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
