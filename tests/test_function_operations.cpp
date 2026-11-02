/***************************************************************************
 *            test_function_operations.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/symbolic_function.hpp"
#include "ariadne-testing.hpp"
#include <sstream>
using namespace Ariadne;
using X=FloatDPApproximation;

class TestFunctionOperations : public Ariadne::Testing::ContextFixture {
public:
    using ContextFixture::ContextFixture;
template<class F> Void check_display(const F& f) {
    std::ostringstream text; text<<f<<" "<<representation(f);
    ARIADNE_TESTING_PROVE(!text.str().empty());
}

template<class P> Void test_constructors() {
    using C=FunctionConstructors<P>;using Y=Number<P>;
    using S=ScalarMultivariateFunction<P>;using V=VectorMultivariateFunction<P>;
    using U=ScalarUnivariateFunction<P>;using W=VectorUnivariateFunction<P>;
    const X scalar(3,dp);Vector<X> point({3,2},dp);
    EuclideanDomain domain(2u);RealDomain real;
    S null_scalar; V null_vector; U null_univariate; W null_curve;
    ARIADNE_TESTING_PROVE(null_scalar.raw_pointer()==nullptr);
    ARIADNE_TESTING_PROVE(null_vector.raw_pointer()==nullptr);
    ARIADNE_TESTING_PROVE(null_univariate.raw_pointer()==nullptr);
    ARIADNE_TESTING_PROVE(null_curve.raw_pointer()==nullptr);
    ARIADNE_TESTING_EXPECT_EQUAL(C::zero()(scalar), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(C::zero(real)(scalar), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(C::constant(Y(2))(scalar), X(2,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(C::constant(real,Y(2))(scalar), X(2,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(C::coordinate()(scalar), scalar);
    ARIADNE_TESTING_EXPECT_EQUAL(C::coordinate(real)(scalar), scalar);
    ARIADNE_TESTING_EXPECT_EQUAL(C::coordinate(real,IndexZero())(scalar), scalar);
    ARIADNE_TESTING_EXPECT_EQUAL(C::identity()(scalar), scalar);
    ARIADNE_TESTING_EXPECT_EQUAL(C::identity(real)(scalar), scalar);
    ARIADNE_TESTING_EXPECT_EQUAL(C::zero(2u)(point), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(C::zero(domain)(point), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(C::constant(2u,Y(2))(point), X(2,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(C::constant(domain,Y(2))(point), X(2,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(C::coordinate(2u,1u)(point), point[1]);
    ARIADNE_TESTING_EXPECT_EQUAL(C::coordinate(domain,1u)(point), point[1]);
    ARIADNE_TESTING_EXPECT_EQUAL(C::coordinates(2u)[1](point), point[1]);
    ARIADNE_TESTING_EXPECT_EQUAL(C::coordinates(domain)[0](point), point[0]);
    ARIADNE_TESTING_PROVE(C::coordinates(0u).empty());
    ARIADNE_TESTING_PROVE(C::coordinates(EuclideanDomain(0u)).empty());
    ARIADNE_TESTING_EXPECT_EQUAL(C::identity(2u)(point), point);
    ARIADNE_TESTING_EXPECT_EQUAL(C::identity(domain)(point), point);
    Vector<X> zeroes(2u,dp);Vector<Y> constants{Y(2),Y(4)};
    Vector<X> expected({2,4},dp);
    ARIADNE_TESTING_EXPECT_EQUAL(C::zeros(2u,2u)(point), zeroes);
    ARIADNE_TESTING_EXPECT_EQUAL(C::zeros(2u,domain)(point), zeroes);
    ARIADNE_TESTING_EXPECT_EQUAL(C::zeros(2u)(scalar), zeroes);
    ARIADNE_TESTING_EXPECT_EQUAL(C::zeros(2u,real)(scalar), zeroes);
    ARIADNE_TESTING_EXPECT_EQUAL(C::constant(2u,constants)(point), expected);
    ARIADNE_TESTING_EXPECT_EQUAL(C::constant(domain,constants)(point), expected);
    ARIADNE_TESTING_EXPECT_EQUAL(C::constant(constants)(scalar), expected);
    ARIADNE_TESTING_EXPECT_EQUAL(C::constant(real,constants)(scalar), expected);
    ARIADNE_TESTING_EXPECT_EQUAL(C::zeros(0u,2u)(point).size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(C::zeros(0u,domain)(point).size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(C::zeros(0u)(scalar).size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(C::zeros(0u,real)(scalar).size(), 0u);
    Vector<X> empty_point(0u,dp);
    ARIADNE_TESTING_EXPECT_EQUAL(C::identity(0u)(empty_point).size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(C::identity(EuclideanDomain(0u))(empty_point).size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(C::constant(2u,Vector<Y>{})(point).size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(C::constant(domain,Vector<Y>{})(point).size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(C::constant(real,Vector<Y>{})(scalar).size(), 0u);
    S s(domain);U u(real);V v(2u,domain);W w(2u,real);
    ARIADNE_TESTING_EXPECT_EQUAL(s(point), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(u(scalar), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(v(point), zeroes);
    ARIADNE_TESTING_EXPECT_EQUAL(w(scalar), zeroes);
    S s2(SizeOne(),2u);U u2{SizeOne(),SizeOne()};V v2(2u,2u);W w2(2u,SizeOne());
    ARIADNE_TESTING_EXPECT_EQUAL(s2(point), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(u2(scalar), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(v2(point), zeroes);
    ARIADNE_TESTING_EXPECT_EQUAL(w2(scalar), zeroes);
    V domain_only(domain);W real_only(real);
    ARIADNE_TESTING_EXPECT_EQUAL(domain_only.result_size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(real_only.result_size(), 0u);
    auto x=C::coordinate(domain,0u);auto t=C::coordinate(real);
    ARIADNE_TESTING_EXPECT_EQUAL(x.create_zero()(point), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(x.create_constant(Y(2))(point), X(2,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(x.create_coordinate(1u)(point), point[1]);
    ARIADNE_TESTING_EXPECT_EQUAL(x.create_constant(constants)(point), expected);
    ARIADNE_TESTING_EXPECT_EQUAL(t.create_zero()(scalar), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(t.create_constant(Y(2))(scalar), X(2,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(t.create_coordinate(IndexZero())(scalar), scalar);
    ARIADNE_TESTING_EXPECT_EQUAL(t.create_constant(constants)(scalar), expected);
    V repeated(2u,x);W repeated_curve(2u,t);
    ARIADNE_TESTING_EXPECT_EQUAL(repeated(point), Vector<X>({3,3},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(repeated_curve(scalar), Vector<X>({3,3},dp));
    V from_list(List<S>{x,x});W curve_list(List<U>{t,t});
    ARIADNE_TESTING_EXPECT_EQUAL(from_list(point), repeated(point));
    ARIADNE_TESTING_EXPECT_EQUAL(curve_list(scalar), repeated_curve(scalar));
    V from_vector(Vector<S>{x,x});W curve_vector(Vector<U>{t,t});
    ARIADNE_TESTING_EXPECT_EQUAL(from_vector(point), repeated(point));
    ARIADNE_TESTING_EXPECT_EQUAL(curve_vector(scalar), repeated_curve(scalar));
    s=Y(4);u=Y(4);v=constants;w=constants;
    ARIADNE_TESTING_EXPECT_EQUAL(s(point), X(4,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(u(scalar), X(4,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(v(point), expected);
    ARIADNE_TESTING_EXPECT_EQUAL(w(scalar), expected);
    ARIADNE_TESTING_EXPECT_EQUAL(x.domain().dimension(), 2u);
    ARIADNE_TESTING_EXPECT_EQUAL(x.codomain().dimension(), 1u);
    ARIADNE_TESTING_EXPECT_EQUAL(v.codomain().dimension(), 2u);
    check_display(x);check_display(v);check_display(t);check_display(w);
}

template<class P> Void test_composition() {
    using C=FunctionConstructors<P>;using Y=Number<P>;
    using V=VectorMultivariateFunction<P>;using W=VectorUnivariateFunction<P>;
    auto x=C::coordinate(2u,0u),y=C::coordinate(2u,1u);auto t=C::coordinate();
    auto f=Y(2)*x-y+Y(1);auto h=sqr(t)+Y(1);
    V g{x+y,x-y};W q{t,t+Y(2)};
    X scalar(3,dp);Vector<X> point({3,2},dp);
    auto hh=compose(h,h),fq=compose(f,q);auto qh=compose(q,h),gq=compose(g,q);
    auto hf=compose(h,f),fg=compose(f,g);auto qf=compose(q,f),gg=compose(g,g);
    ARIADNE_TESTING_EXPECT_EQUAL(hh(scalar), X(101,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(fq(scalar), X(2,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(qh(scalar), Vector<X>({10,12},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(gq(scalar), Vector<X>({8,-2},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(hf(point), X(26,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(fg(point), X(10,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(qf(point), Vector<X>({5,7},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(gg(point), Vector<X>({6,4},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(derivative(hh)(scalar), X(120,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(derivative(fq)(scalar), X(1,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(derivative(qh)(scalar), Vector<X>({6,6},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(derivative(gq)(scalar), Vector<X>({2,0},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(derivative(hf,0u)(point), X(20,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(derivative(fg,0u)(point), X(1,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(derivative(qf,0u)(point), Vector<X>({2,2},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(derivative(gg,0u)(point), Vector<X>({2,0},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(h.slope(scalar), X(6,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(q.tangent(scalar), Vector<X>({1,1},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(f.gradient(point)[1], X(-1,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(g.jacobian(point)[1][1], X(-1,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(qh.get(1u)(scalar), X(12,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(gg.get(1u)(point), X(4,dp));
    check_display(hh);check_display(fq);check_display(qh);check_display(gq);
    check_display(hf);check_display(fg);check_display(qf);check_display(gg);
    auto both=join(g,g),head=join(x,g),tail=join(g,y),pair=join(x,y);
    ARIADNE_TESTING_EXPECT_EQUAL(both(point), Vector<X>({5,1,5,1},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(head(point), Vector<X>({3,5,1},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(tail(point), Vector<X>({5,1,2},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(pair(point), point);
    ARIADNE_TESTING_EXPECT_EQUAL(both.get(3u)(point), X(1,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(derivative(both,1u)(point), Vector<X>({1,-1,1,-1},dp));
    check_display(both);
    auto scalar_expression=pow(x,3)+Y(2)-Y(4)*x/Y(2);
    ARIADNE_TESTING_EXPECT_EQUAL(scalar_expression(point), X(23,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((Y(8)/x)(point), X(8,dp)/X(3,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((x/Y(2))(point), X(1.5_x,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((Y(2)-x)(point), X(-1,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((-x)(point), X(-3,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((+x)(point), X(3,dp));
    ARIADNE_TESTING_EXPECT_EQUAL((Y(2)+x)(point), X(5,dp));
}
Void test_effective_vector_algebra() {
    using S=EffectiveScalarMultivariateFunction;using V=EffectiveVectorMultivariateFunction;
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    auto scalar=S::constant(2u,2);V vector{x+y,x-y};
    Vector<X> point({3,1},dp),expected({4,2},dp),constants({1,2},dp);
    Vector<EffectiveNumber> c{1,2};
    ARIADNE_TESTING_EXPECT_EQUAL((+vector)(point), expected);
    ARIADNE_TESTING_EXPECT_EQUAL((-vector)(point), -expected);
    ARIADNE_TESTING_EXPECT_EQUAL((vector+vector)(point), Vector<X>({8,4},dp));
    ARIADNE_TESTING_EXPECT_EQUAL((vector-vector)(point), Vector<X>({0,0},dp));
    ARIADNE_TESTING_EXPECT_EQUAL((scalar*vector)(point), Vector<X>({8,4},dp));
    ARIADNE_TESTING_EXPECT_EQUAL((vector*scalar)(point), Vector<X>({8,4},dp));
    ARIADNE_TESTING_EXPECT_EQUAL((vector/scalar)(point), Vector<X>({2,1},dp));
    ARIADNE_TESTING_EXPECT_EQUAL((c+vector)(point), expected+constants);
    ARIADNE_TESTING_EXPECT_EQUAL((vector+c)(point), expected+constants);
    ARIADNE_TESTING_EXPECT_EQUAL((c-vector)(point), constants-expected);
    ARIADNE_TESTING_EXPECT_EQUAL((vector-c)(point), expected-constants);
    ARIADNE_TESTING_EXPECT_EQUAL((EffectiveNumber(2)*vector)(point), Vector<X>({8,4},dp));
    ARIADNE_TESTING_EXPECT_EQUAL((vector*EffectiveNumber(2))(point), Vector<X>({8,4},dp));
    ARIADNE_TESTING_EXPECT_EQUAL((vector/EffectiveNumber(2))(point), Vector<X>({2,1},dp));
    auto t=EffectiveScalarUnivariateFunction::coordinate();
    EffectiveVectorUnivariateFunction curve{t,t+2};X parameter(3,dp);
    ARIADNE_TESTING_EXPECT_EQUAL((curve+curve)(parameter), Vector<X>({6,10},dp));
    ARIADNE_TESTING_EXPECT_EQUAL((curve*t)(parameter), Vector<X>({9,15},dp));
    auto embedded_scalar=embed(1u,x+y,1u);auto embedded_vector=embed(1u,vector,1u);
    Vector<X> extended({99,3,1,88},dp);
    ARIADNE_TESTING_EXPECT_EQUAL(embedded_scalar(extended), X(4,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(embedded_vector(extended), expected);
    ARIADNE_TESTING_EXPECT_EQUAL(embedded_scalar.derivative(0u)(extended), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(embedded_scalar.derivative(1u)(extended), X(1,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(embedded_scalar.derivative(3u)(extended), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(embedded_vector.derivative(0u)(extended), Vector<X>({0,0},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(embedded_vector.derivative(1u)(extended), Vector<X>({1,1},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(embedded_vector.get(1u)(extended), X(2,dp));
    auto f=sqr(x)+y;
    ARIADNE_TESTING_EXPECT_EQUAL(derivatives(f)(point), Vector<X>({6,1},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(lie_derivative(f,vector)(point), X(26,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(lie_derivative(vector,vector)(point), Vector<X>({6,2},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(evaluate(make_formula(f),point), X(10,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(evaluate(make_formula(vector),point), expected);
    auto fx=Formula<EffectiveNumber>::coordinate(0u),fy=Formula<EffectiveNumber>::coordinate(1u);
    V formula_vector(EuclideanDomain(2u),Vector<Formula<EffectiveNumber>>{fx+fy,fx-fy});
    ARIADNE_TESTING_EXPECT_EQUAL(evaluate(make_formula(formula_vector),point), expected);
    ARIADNE_TESTING_EXPECT_EQUAL(evaluate(make_formula(compose(vector,vector)),point), Vector<X>({6,2},dp));
    auto zero_vector=V::zeros(0u,0u);
    ARIADNE_TESTING_EXPECT_EQUAL(make_formula(zero_vector).size(), 0u);
    check_display(embedded_scalar);check_display(embedded_vector);
}
Void test_validated_curve_algebra() {
    auto t=ValidatedScalarUnivariateFunction::coordinate();
    ValidatedVectorUnivariateFunction curve{t,t+ValidatedNumber(2)};X parameter(3,dp);
    ARIADNE_TESTING_EXPECT_EQUAL((curve+curve)(parameter), Vector<X>({6,10},dp));
    ARIADNE_TESTING_EXPECT_EQUAL((curve*t)(parameter), Vector<X>({9,15},dp));
}
Void test_general_scalar_vector_product() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
    Vector<EffectiveNumber> coefficients{2,-3};
    Vector<X> point({3,1},dp);
    ARIADNE_TESTING_EXPECT_EQUAL((x*coefficients)(point), Vector<X>({6,-9},dp));
    ARIADNE_TESTING_EXPECT_EQUAL((x*Vector<EffectiveNumber>{})(point).size(), 0u);
}
Void test_zero_dimensional_derivatives() {
    using S=EffectiveScalarMultivariateFunction;using V=EffectiveVectorMultivariateFunction;
    auto constant=S::constant(0u,5);auto field=V::identity(0u);
    Vector<X> point(0u,dp);
    ARIADNE_TESTING_EXPECT_EQUAL(derivatives(constant)(point).size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(lie_derivative(constant,field)(point), X(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(lie_derivative(V::constant(0u,Vector<EffectiveNumber>{2,3}),field)(point), Vector<X>({0,0},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(lie_derivative(field,field)(point).size(), 0u);
}
    void test_constructors_approximatetag() { test_constructors<ApproximateTag>(); }
    void test_constructors_validatedtag() { test_constructors<ValidatedTag>(); }
    void test_constructors_effectivetag() { test_constructors<EffectiveTag>(); }
    void test_composition_approximatetag() { test_composition<ApproximateTag>(); }
    void test_composition_validatedtag() { test_composition<ValidatedTag>(); }
    void test_composition_effectivetag() { test_composition<EffectiveTag>(); }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_function_operations");
    ARIADNE_TESTING_ADD_TEST(TestFunctionOperations, test_constructors_approximatetag);
    ARIADNE_TESTING_ADD_TEST(TestFunctionOperations, test_constructors_validatedtag);
    ARIADNE_TESTING_ADD_TEST(TestFunctionOperations, test_constructors_effectivetag);
    ARIADNE_TESTING_ADD_TEST(TestFunctionOperations, test_composition_approximatetag);
    ARIADNE_TESTING_ADD_TEST(TestFunctionOperations, test_composition_validatedtag);
    ARIADNE_TESTING_ADD_TEST(TestFunctionOperations, test_composition_effectivetag);
    ARIADNE_TESTING_ADD_TEST(TestFunctionOperations, test_general_scalar_vector_product);
    ARIADNE_TESTING_ADD_TEST(TestFunctionOperations, test_zero_dimensional_derivatives);
    ARIADNE_TESTING_ADD_TEST(TestFunctionOperations, test_effective_vector_algebra);
    ARIADNE_TESTING_ADD_TEST(TestFunctionOperations, test_validated_curve_algebra);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
