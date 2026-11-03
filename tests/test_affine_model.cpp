/***************************************************************************
 *            test_affine_model.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/affine_model.hpp"
#include "function/affine.hpp"
#include "function/taylor_model.hpp"
#include "function/function.hpp"
#include "function/affine_model.tpl.hpp"
#include "ariadne-testing.hpp"
#include <sstream>
using namespace Ariadne;

class TestAffineModel : public Ariadne::Testing::ContextFixture {
public:
    using ContextFixture::ContextFixture;
template<class F> Void test_product(typename F::PrecisionType pr) {
    using A=AffineModel<ApproximateTag,F>;using X=Approximation<F>;
    A a(X(2,pr),{X(3,pr),X(-1,pr)}),b(X(5,pr),{X(7,pr),X(2,pr)});
    auto product=a*b;
    ARIADNE_TESTING_EXPECT_SAME(product.value(), X(10,pr));
    ARIADNE_TESTING_EXPECT_SAME(product.gradient(0u), X(29,pr));
    ARIADNE_TESTING_EXPECT_SAME(product.gradient(1u), X(-1,pr));
    auto constant=A::constant(0u,X(2,pr))*A::constant(0u,X(5,pr));
    ARIADNE_TESTING_EXPECT_SAME(constant.value(), X(10,pr));
    using V=AffineModel<ValidatedTag,F>;using E=typename V::ErrorType;
    V factor(1u,F(2,pr)),uncertain(1u,F(0,pr));uncertain.set_error(E(3u,pr));
    auto enclosure=factor*uncertain;
    ARIADNE_TESTING_PROVE(decide(enclosure.error().raw()>=F(6,pr)));
    ARIADNE_TESTING_PROVE(decide(enclosure.range().lower_bound()<=F(-6,pr)));
    ARIADNE_TESTING_PROVE(decide(enclosure.range().upper_bound()>=F(6,pr)));
    ARIADNE_TESTING_PROVE(decide((uncertain*factor).error().raw()>=F(6,pr)));
}

template<class F> Void test_approximate_model(typename F::PrecisionType pr) {
    using A=AffineModel<ApproximateTag,F>;using X=Approximation<F>;
    A a(X(2,pr),{X(3,pr),X(-1,pr),X(0,pr)});
    Vector<X> point({0.25_x,0.5_x,0.0_x},pr);
    ARIADNE_TESTING_EXPECT_SAME(a.evaluate(point), X(2.25_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(a.range().lower_bound(), X(-2,pr));
    ARIADNE_TESTING_EXPECT_SAME(a.range().upper_bound(), X(6,pr));
    ARIADNE_TESTING_EXPECT_SAME(a.b(), a.value());
    ARIADNE_TESTING_EXPECT_RANGES_BY(a.a(), a.gradient(), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_EQUAL(a.precision(), pr);
    ARIADNE_TESTING_EXPECT_EQUAL(std::get<0>(a.characteristics()), 3u);
    ARIADNE_TESTING_EXPECT_SAME(a.create_zero().evaluate(point), X(0,pr));
    const A& ca=a;
    ARIADNE_TESTING_EXPECT_SAME(ca[0u], X(3,pr));
    a[0u]=X(4,pr);a.set_value(X(3,pr));
    a.set_gradient(0u,X(3,pr));a.set_gradient(1u,ApproximateNumber(-1));
    ARIADNE_TESTING_EXPECT_SAME(a.evaluate(point), X(3.25_x,pr));
    auto b=A::coordinate(3u,1u,pr);X coefficient(2,pr);
    ARIADNE_TESTING_EXPECT_SAME((-a).evaluate(point), X(-3.25_x,pr));
    ARIADNE_TESTING_EXPECT_SAME((a+b).evaluate(point), X(3.75_x,pr));
    ARIADNE_TESTING_EXPECT_SAME((coefficient+a).evaluate(point), X(5.25_x,pr));
    ARIADNE_TESTING_EXPECT_SAME((a+coefficient).evaluate(point), X(5.25_x,pr));
    ARIADNE_TESTING_EXPECT_SAME((coefficient*a).evaluate(point), X(6.5_x,pr));
    ARIADNE_TESTING_EXPECT_SAME((a*coefficient).evaluate(point), X(6.5_x,pr));
    A from_covector(X(3,pr),Covector<X>({X(3,pr),X(-1,pr),X(0,pr)}));
    ARIADNE_TESTING_EXPECT_SAME(from_covector.evaluate(point), a.evaluate(point));
    a=coefficient;
    ARIADNE_TESTING_EXPECT_SAME(a.evaluate(point), coefficient);
    A constant(0u,coefficient);constant=X(5,pr);
    Vector<X> empty(0u,pr);
    ARIADNE_TESTING_EXPECT_SAME(constant.evaluate(empty), X(5,pr));
    ARIADNE_TESTING_EXPECT_SAME(constant.range().lower_bound(), X(5,pr));
    ARIADNE_TESTING_EXPECT_EXCEPTION(a+A(2u,pr));
    ARIADNE_TESTING_EXPECT_EXCEPTION(a*A(2u,pr));
    std::ostringstream text;text<<from_covector<<constant;
    ARIADNE_TESTING_INSPECT("text.str()", text.str());
    ARIADNE_TESTING_PROVE(text.str().find("*x0")!=String::npos);
    ARIADNE_TESTING_PROVE(text.str().find("*x1")!=String::npos);
    ARIADNE_TESTING_PROVE(text.str().find("*x2")==String::npos);
    BoxDomainType domain{{2,6},{-4,2}};
    auto scalings=A::scalings(domain,pr);
    ARIADNE_TESTING_EXPECT_SAME(scalings[0u].evaluate(Vector<X>({-1,1},pr)), X(2,pr));
    ARIADNE_TESTING_EXPECT_SAME(scalings[1u].evaluate(Vector<X>({-1,1},pr)), X(2,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(A::scalings(BoxDomainType(0u),pr).size(), 0u);
}
template<class F> Void test_validated_model(typename F::PrecisionType pr) {
    using A=AffineModel<ValidatedTag,F>;using X=Bounds<F>;using E=typename A::ErrorType;
    A a(F(2,pr),{F(3,pr),F(-1,pr),F(0,pr)});a.set_error(E(F(0.5_x,pr)));
    const A& ca=a;
    // The natural-number overload must replace the error, not accumulate it.
    a.set_error(2u);
    ARIADNE_TESTING_EXPECT_SAME(a.error().raw(), F(2,pr));
    ARIADNE_TESTING_PROVE(decide(a.range().lower_bound().raw()==F(-4,pr)));
    ARIADNE_TESTING_PROVE(decide(a.range().upper_bound().raw()==F(8,pr)));
    a.set_error(E(F(0.5_x,pr)));
    Vector<X> point({0.25_x,0.5_x,0.0_x},pr);
    auto value=a.evaluate(point);
    ARIADNE_TESTING_EXPECT_SAME(value, X(1.75_x,2.75_x,pr));
    ARIADNE_TESTING_PROVE(decide(a.range().lower_bound().raw()==F(-2.5_x,pr)));
    ARIADNE_TESTING_PROVE(decide(a.range().upper_bound().raw()==F(6.5_x,pr)));
    ARIADNE_TESTING_EXPECT_SAME(ca.b(), ca.value());ARIADNE_TESTING_EXPECT_RANGES_BY(ca.a(), ca.gradient(), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_SAME(ca.e().raw(), ca.error().raw());
    ARIADNE_TESTING_EXPECT_SAME(ca[0u], F(3,pr));ARIADNE_TESTING_EXPECT_SAME(ca.gradient(1u), F(-1,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(std::get<0>(ca.characteristics()), 3u);
    ARIADNE_TESTING_EXPECT_EQUAL(ca.precision(), pr);
    ARIADNE_TESTING_PROVE(definitely(a.create_zero().evaluate(point)==X(0,pr)));
    a.set_value(F(3,pr));a.set_gradient(0u,F(3,pr));a.set_gradient(1u,Dyadic(-1));a[2u]=F(0,pr);
    A from_covector(F(3,pr),Covector<F>({F(3,pr),F(-1,pr),F(0,pr)}),E(F(0.5_x,pr)));
    ARIADNE_TESTING_EXPECT_SAME(from_covector.evaluate(point), a.evaluate(point));
    auto b=A::coordinate(3u,1u,pr);X coefficient(2,pr);
    ARIADNE_TESTING_EXPECT_SAME((-a).evaluate(point), X(-3.75_x,-2.75_x,pr));
    ARIADNE_TESTING_EXPECT_SAME((a+b).evaluate(point), X(3.25_x,4.25_x,pr));
    ARIADNE_TESTING_EXPECT_SAME((coefficient+a).evaluate(point), X(4.75_x,5.75_x,pr));
    ARIADNE_TESTING_EXPECT_SAME((a+coefficient).evaluate(point), X(4.75_x,5.75_x,pr));
    ARIADNE_TESTING_EXPECT_SAME((coefficient*a).evaluate(point), X(5.5_x,7.5_x,pr));
    ARIADNE_TESTING_EXPECT_SAME((a*coefficient).evaluate(point), X(5.5_x,7.5_x,pr));
    a=F(2,pr);ARIADNE_TESTING_EXPECT_SAME(a.error().raw(), F(0,pr));
    ARIADNE_TESTING_PROVE(definitely(a.evaluate(point)==X(2,pr)));
    A constant=A::constant(0u,F(2,pr));constant=F(5,pr);
    Vector<X> empty(0u,pr);
    ARIADNE_TESTING_PROVE(definitely(constant.evaluate(empty)==X(5,pr)));
    ARIADNE_TESTING_EXPECT_SAME((constant*constant).value(), F(25,pr));
    ARIADNE_TESTING_EXPECT_SAME((constant*constant).error().raw(), F(0,pr));
    ARIADNE_TESTING_EXPECT_EXCEPTION(a+A(2u,pr));ARIADNE_TESTING_EXPECT_EXCEPTION(a*A(2u,pr));
    std::ostringstream text;text<<from_covector<<constant;
    ARIADNE_TESTING_INSPECT("text.str()", text.str());
    ARIADNE_TESTING_PROVE(text.str().find("*x0")!=String::npos);
    ARIADNE_TESTING_PROVE(text.str().find("*x1")!=String::npos);
    ARIADNE_TESTING_PROVE(text.str().find("*x2")==String::npos);
    BoxDomainType domain{{2,6},{-4,2}};
    auto scalings=A::scalings(domain,pr);
    ARIADNE_TESTING_PROVE(definitely(scalings[0u].evaluate(Vector<X>({-1,1},pr))==X(2,pr)));
    ARIADNE_TESTING_PROVE(definitely(scalings[1u].evaluate(Vector<X>({-1,1},pr))==X(2,pr)));
    ARIADNE_TESTING_EXPECT_EQUAL(A::scalings(BoxDomainType(0u),pr).size(), 0u);
    Affine<X> uncertain(X(1,3,pr),{X(2,4,pr),X(-2,0,pr)});
    auto old_rounding=F::get_rounding_mode();
    auto converted=affine_model(uncertain);
    ARIADNE_TESTING_PROVE(F::get_rounding_mode()==old_rounding);
    ARIADNE_TESTING_EXPECT_SAME(converted.value(), F(2,pr));
    ARIADNE_TESTING_EXPECT_SAME(converted.gradient(0u), F(3,pr));
    ARIADNE_TESTING_EXPECT_SAME(converted.error().raw(), F(3,pr));
    Sweeper<F> sweeper=ThresholdSweeper<F>(pr,1e-12);
    Expansion<MultiIndex,F> expansion({{{0,0},2.0_x},{{1,0},3.0_x},{{0,1},-1.0_x},{{2,0},0.5_x}},pr);
    TaylorModel<ValidatedTag,F> taylor(expansion,F(0.25_x,pr),sweeper);
    auto linear=affine_model(taylor);
    ARIADNE_TESTING_PROVE(F::get_rounding_mode()==old_rounding);
    ARIADNE_TESTING_EXPECT_SAME(linear.value(), F(2,pr));
    ARIADNE_TESTING_EXPECT_SAME(linear.gradient(1u), F(-1,pr));
    ARIADNE_TESTING_EXPECT_SAME(linear.error().raw(), F(0.75_x,pr));
    auto models=affine_models(Vector<TaylorModel<ValidatedTag,F>>{taylor});
    ARIADNE_TESTING_EXPECT_SAME(models[0u].value(), F(2,pr));
}
template<class P,class F> Void test_composition(typename F::PrecisionType pr) {
    using A=AffineModel<P,F>;using C=typename A::CoefficientType;using X=typename A::NumericType;
    using S=ScalarMultivariateFunction<P>;using V=VectorMultivariateFunction<P>;
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    Vector<A> input{A(C(2,pr),{C(1,pr),C(0,pr)}),A(C(3,pr),{C(0,pr),C(1,pr)})};
    auto scalar=compose(x+y,input);auto vector=compose(V{x+y,x-y},input);
    Vector<X> point({0.25_x,0.5_x},pr);
    ARIADNE_TESTING_PROVE(decide(scalar.evaluate(point)==X(5.75_x,pr)));
    ARIADNE_TESTING_PROVE(decide(vector[1u].evaluate(point)==X(-1.25_x,pr)));
    auto nonlinear=compose(sqr(x)+y,input);
    ARIADNE_TESTING_EXPECT_SAME(nonlinear.value(), C(7,pr));
    ARIADNE_TESTING_EXPECT_SAME(nonlinear.gradient(1u), C(1,pr));
}
Void test_function_conversion() {
    MultiplePrecision mp(128u);
    ARIADNE_TESTING_PROVE(decide(FloatMPApproximation(3,mp)<FloatDP(4,dp)));
    ARIADNE_TESTING_PROVE(decide(FloatMP(3,mp)<FloatDP(4,dp)));
    BoxDomainType domain{{-1,1},{-1,1}};
    auto x=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
    auto y=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
    auto scalar=affine_model(domain,x+y,dp);
    auto vector=affine_models(domain,EffectiveVectorMultivariateFunction{x+y,x-y},dp);
    Vector<FloatDPBounds> point({0.25_x,0.5_x},dp);
    ARIADNE_TESTING_PROVE(definitely(scalar.evaluate(point)==FloatDPBounds(0.75_x,dp)));
    ARIADNE_TESTING_PROVE(definitely(vector[1u].evaluate(point)==FloatDPBounds(-0.25_x,dp)));
}
    void test_product_floatdp_0() { test_product<FloatDP>(dp); }
    void test_product_floatmp_1() { test_product<FloatMP>(MultiplePrecision(128u)); }
    void test_approximate_model_floatdp_2() { test_approximate_model<FloatDP>(dp); }
    void test_approximate_model_floatmp_3() { test_approximate_model<FloatMP>(MultiplePrecision(128u)); }
    void test_validated_model_floatdp_4() { test_validated_model<FloatDP>(dp); }
    void test_validated_model_floatmp_5() { test_validated_model<FloatMP>(MultiplePrecision(128u)); }
    void test_composition_approximatetag_floatdp_6() { test_composition<ApproximateTag,FloatDP>(dp); }
    void test_composition_validatedtag_floatdp_7() { test_composition<ValidatedTag,FloatDP>(dp); }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_affine_model");
    ARIADNE_TESTING_ADD_TEST(TestAffineModel, test_product_floatdp_0);
    ARIADNE_TESTING_ADD_TEST(TestAffineModel, test_product_floatmp_1);
    ARIADNE_TESTING_ADD_TEST(TestAffineModel, test_approximate_model_floatdp_2);
    ARIADNE_TESTING_ADD_TEST(TestAffineModel, test_approximate_model_floatmp_3);
    ARIADNE_TESTING_ADD_TEST(TestAffineModel, test_validated_model_floatdp_4);
    ARIADNE_TESTING_ADD_TEST(TestAffineModel, test_validated_model_floatmp_5);
    ARIADNE_TESTING_ADD_TEST(TestAffineModel, test_composition_approximatetag_floatdp_6);
    ARIADNE_TESTING_ADD_TEST(TestAffineModel, test_composition_validatedtag_floatdp_7);
    ARIADNE_TESTING_ADD_TEST(TestAffineModel, test_function_conversion);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
