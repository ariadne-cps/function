/***************************************************************************
 *            test_function_operations.cpp
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
    ARIADNE_TESTING_EXPECT_SAME(C::zero()(scalar), X(0,dp));
    ARIADNE_TESTING_EXPECT_SAME(C::zero(real)(scalar), X(0,dp));
    ARIADNE_TESTING_EXPECT_SAME(C::constant(Y(2))(scalar), X(2,dp));
    ARIADNE_TESTING_EXPECT_SAME(C::constant(real,Y(2))(scalar), X(2,dp));
    ARIADNE_TESTING_EXPECT_SAME(C::coordinate()(scalar), scalar);
    ARIADNE_TESTING_EXPECT_SAME(C::coordinate(real)(scalar), scalar);
    ARIADNE_TESTING_EXPECT_SAME(C::coordinate(real,IndexZero())(scalar), scalar);
    ARIADNE_TESTING_EXPECT_SAME(C::identity()(scalar), scalar);
    ARIADNE_TESTING_EXPECT_SAME(C::identity(real)(scalar), scalar);
    ARIADNE_TESTING_EXPECT_SAME(C::zero(2u)(point), X(0,dp));
    ARIADNE_TESTING_EXPECT_SAME(C::zero(domain)(point), X(0,dp));
    ARIADNE_TESTING_EXPECT_SAME(C::constant(2u,Y(2))(point), X(2,dp));
    ARIADNE_TESTING_EXPECT_SAME(C::constant(domain,Y(2))(point), X(2,dp));
    ARIADNE_TESTING_EXPECT_SAME(C::coordinate(2u,1u)(point), point[1]);
    ARIADNE_TESTING_EXPECT_SAME(C::coordinate(domain,1u)(point), point[1]);
    ARIADNE_TESTING_EXPECT_SAME(C::coordinates(2u)[1](point), point[1]);
    ARIADNE_TESTING_EXPECT_SAME(C::coordinates(domain)[0](point), point[0]);
    ARIADNE_TESTING_PROVE(C::coordinates(0u).empty());
    ARIADNE_TESTING_PROVE(C::coordinates(EuclideanDomain(0u)).empty());
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (C::identity(2u)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (point); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (C::identity(domain)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (point); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    Vector<X> zeroes(2u,dp);Vector<Y> constants{Y(2),Y(4)};
    Vector<X> expected({2,4},dp);
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (C::zeros(2u,2u)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (zeroes); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (C::zeros(2u,domain)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (zeroes); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (C::zeros(2u)(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (zeroes); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (C::zeros(2u,real)(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (zeroes); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (C::constant(2u,constants)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (C::constant(domain,constants)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (C::constant(constants)(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (C::constant(real,constants)(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
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
    ARIADNE_TESTING_EXPECT_SAME(s(point), X(0,dp));
    ARIADNE_TESTING_EXPECT_SAME(u(scalar), X(0,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (v(point)); return snapshot.array(); }()), ([&] { auto snapshot = (zeroes); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (w(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (zeroes); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    S s2(SizeOne(),2u);U u2{SizeOne(),SizeOne()};V v2(2u,2u);W w2(2u,SizeOne());
    ARIADNE_TESTING_EXPECT_SAME(s2(point), X(0,dp));
    ARIADNE_TESTING_EXPECT_SAME(u2(scalar), X(0,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (v2(point)); return snapshot.array(); }()), ([&] { auto snapshot = (zeroes); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (w2(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (zeroes); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    V domain_only(domain);W real_only(real);
    ARIADNE_TESTING_EXPECT_EQUAL(domain_only.result_size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(real_only.result_size(), 0u);
    auto x=C::coordinate(domain,0u);auto t=C::coordinate(real);
    ARIADNE_TESTING_EXPECT_SAME(x.create_zero()(point), X(0,dp));
    ARIADNE_TESTING_EXPECT_SAME(x.create_constant(Y(2))(point), X(2,dp));
    ARIADNE_TESTING_EXPECT_SAME(x.create_coordinate(1u)(point), point[1]);
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (x.create_constant(constants)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_SAME(t.create_zero()(scalar), X(0,dp));
    ARIADNE_TESTING_EXPECT_SAME(t.create_constant(Y(2))(scalar), X(2,dp));
    ARIADNE_TESTING_EXPECT_SAME(t.create_coordinate(IndexZero())(scalar), scalar);
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (t.create_constant(constants)(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    V repeated(2u,x);W repeated_curve(2u,t);
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (repeated(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({3,3},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (repeated_curve(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({3,3},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    V from_list(List<S>{x,x});W curve_list(List<U>{t,t});
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (from_list(point)); return snapshot.array(); }()), ([&] { auto snapshot = (repeated(point)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (curve_list(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (repeated_curve(scalar)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    V from_vector(Vector<S>{x,x});W curve_vector(Vector<U>{t,t});
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (from_vector(point)); return snapshot.array(); }()), ([&] { auto snapshot = (repeated(point)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (curve_vector(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (repeated_curve(scalar)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    s=Y(4);u=Y(4);v=constants;w=constants;
    ARIADNE_TESTING_EXPECT_SAME(s(point), X(4,dp));
    ARIADNE_TESTING_EXPECT_SAME(u(scalar), X(4,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (v(point)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (w(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
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
    ARIADNE_TESTING_EXPECT_SAME(hh(scalar), X(101,dp));
    ARIADNE_TESTING_EXPECT_SAME(fq(scalar), X(2,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (qh(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({10,12},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (gq(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({8,-2},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_SAME(hf(point), X(26,dp));
    ARIADNE_TESTING_EXPECT_SAME(fg(point), X(10,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (qf(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({5,7},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (gg(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({6,4},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_SAME(derivative(hh)(scalar), X(120,dp));
    ARIADNE_TESTING_EXPECT_SAME(derivative(fq)(scalar), X(1,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (derivative(qh)(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({6,6},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (derivative(gq)(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({2,0},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_SAME(derivative(hf,0u)(point), X(20,dp));
    ARIADNE_TESTING_EXPECT_SAME(derivative(fg,0u)(point), X(1,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (derivative(qf,0u)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({2,2},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (derivative(gg,0u)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({2,0},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_SAME(h.slope(scalar), X(6,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (q.tangent(scalar)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({1,1},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    const auto f_gradient=f.gradient(point);
    const auto g_jacobian=g.jacobian(point);
    ARIADNE_TESTING_EXPECT_SAME(f_gradient[1], X(-1,dp));
    ARIADNE_TESTING_EXPECT_SAME(g_jacobian[1][1], X(-1,dp));
    ARIADNE_TESTING_EXPECT_SAME(qh.get(1u)(scalar), X(12,dp));
    ARIADNE_TESTING_EXPECT_SAME(gg.get(1u)(point), X(4,dp));
    check_display(hh);check_display(fq);check_display(qh);check_display(gq);
    check_display(hf);check_display(fg);check_display(qf);check_display(gg);
    auto both=join(g,g),head=join(x,g),tail=join(g,y),pair=join(x,y);
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (both(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({5,1,5,1},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (head(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({3,5,1},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (tail(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({5,1,2},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (pair(point)); return snapshot.array(); }()), ([&] { auto snapshot = (point); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_SAME(both.get(3u)(point), X(1,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (derivative(both,1u)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({1,-1,1,-1},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    check_display(both);
    auto scalar_expression=pow(x,3)+Y(2)-Y(4)*x/Y(2);
    ARIADNE_TESTING_EXPECT_SAME(scalar_expression(point), X(23,dp));
    ARIADNE_TESTING_EXPECT_SAME((Y(8)/x)(point), X(8,dp)/X(3,dp));
    ARIADNE_TESTING_EXPECT_SAME((x/Y(2))(point), X(1.5_x,dp));
    ARIADNE_TESTING_EXPECT_SAME((Y(2)-x)(point), X(-1,dp));
    ARIADNE_TESTING_EXPECT_SAME((-x)(point), X(-3,dp));
    ARIADNE_TESTING_EXPECT_SAME((+x)(point), X(3,dp));
    ARIADNE_TESTING_EXPECT_SAME((Y(2)+x)(point), X(5,dp));
}
Void test_effective_vector_algebra() {
    using S=EffectiveScalarMultivariateFunction;using V=EffectiveVectorMultivariateFunction;
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    auto scalar=S::constant(2u,2);V vector{x+y,x-y};
    Vector<X> point({3,1},dp),expected({4,2},dp),constants({1,2},dp);
    Vector<EffectiveNumber> c{1,2};
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((+vector)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((-vector)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (-expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((vector+vector)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({8,4},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((vector-vector)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({0,0},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((scalar*vector)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({8,4},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((vector*scalar)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({8,4},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((vector/scalar)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({2,1},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((c+vector)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (expected+constants); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((vector+c)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (expected+constants); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((c-vector)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (constants-expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((vector-c)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (expected-constants); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((EffectiveNumber(2)*vector)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({8,4},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((vector*EffectiveNumber(2))(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({8,4},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((vector/EffectiveNumber(2))(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({2,1},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    auto t=EffectiveScalarUnivariateFunction::coordinate();
    EffectiveVectorUnivariateFunction curve{t,t+2};X parameter(3,dp);
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((curve+curve)(parameter)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({6,10},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((curve*t)(parameter)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({9,15},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    auto embedded_scalar=embed(1u,x+y,1u);auto embedded_vector=embed(1u,vector,1u);
    Vector<X> extended({99,3,1,88},dp);
    ARIADNE_TESTING_EXPECT_SAME(embedded_scalar(extended), X(4,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (embedded_vector(extended)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_SAME(embedded_scalar.derivative(0u)(extended), X(0,dp));
    ARIADNE_TESTING_EXPECT_SAME(embedded_scalar.derivative(1u)(extended), X(1,dp));
    ARIADNE_TESTING_EXPECT_SAME(embedded_scalar.derivative(3u)(extended), X(0,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (embedded_vector.derivative(0u)(extended)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({0,0},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (embedded_vector.derivative(1u)(extended)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({1,1},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_SAME(embedded_vector.get(1u)(extended), X(2,dp));
    auto f=sqr(x)+y;
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (derivatives(f)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({6,1},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_SAME(lie_derivative(f,vector)(point), X(26,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (lie_derivative(vector,vector)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({6,2},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_SAME(evaluate(make_formula(f),point), X(10,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (evaluate(make_formula(vector),point)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    auto fx=Formula<EffectiveNumber>::coordinate(0u),fy=Formula<EffectiveNumber>::coordinate(1u);
    V formula_vector(EuclideanDomain(2u),Vector<Formula<EffectiveNumber>>{fx+fy,fx-fy});
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (evaluate(make_formula(formula_vector),point)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (evaluate(make_formula(compose(vector,vector)),point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({6,2},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    auto zero_vector=V::zeros(0u,0u);
    ARIADNE_TESTING_EXPECT_EQUAL(make_formula(zero_vector).size(), 0u);
    check_display(embedded_scalar);check_display(embedded_vector);
}
Void test_validated_curve_algebra() {
    auto t=ValidatedScalarUnivariateFunction::coordinate();
    ValidatedVectorUnivariateFunction curve{t,t+ValidatedNumber(2)};X parameter(3,dp);
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((curve+curve)(parameter)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({6,10},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((curve*t)(parameter)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({9,15},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
}
Void test_general_scalar_vector_product() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
    Vector<EffectiveNumber> coefficients{2,-3};
    Vector<X> point({3,1},dp);
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = ((x*coefficients)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({6,-9},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_EQUAL((x*Vector<EffectiveNumber>{})(point).size(), 0u);
}
Void test_zero_dimensional_derivatives() {
    using S=EffectiveScalarMultivariateFunction;using V=EffectiveVectorMultivariateFunction;
    auto constant=S::constant(0u,5);auto field=V::identity(0u);
    Vector<X> point(0u,dp);
    ARIADNE_TESTING_EXPECT_EQUAL(derivatives(constant)(point).size(), 0u);
    ARIADNE_TESTING_EXPECT_SAME(lie_derivative(constant,field)(point), X(0,dp));
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (lie_derivative(V::constant(0u,Vector<EffectiveNumber>{2,3}),field)(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<X>({0,0},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
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
