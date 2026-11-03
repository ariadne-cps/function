/***************************************************************************
 *            test_function_interface.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/function_mixin.tpl.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.tpl.hpp"
#include "algebra/algebra.hpp"
#include "ariadne-testing.hpp"
#include <array>
#include <memory>
#include <sstream>
using namespace Ariadne;

class TestFunctionInterface : public Ariadne::Testing::ContextFixture {
public:
    using ContextFixture::ContextFixture;
    template<class P> void check_vector_api();
    void test_optional_getter();
    void test_function_factory();
    template<class Y,class P> void check_generic_polynomial(const ScalarMultivariateFunction<P>& function);
    template<class P> void test_generic_numeric_interfaces();
    template<class MT,class P> void check_taylor_argument(const ScalarMultivariateFunction<P>& function, typename MT::PrecisionType pr);
    template<class P> void test_numeric_jets_and_taylor_arguments();
    template<class P> void test_handle_derivative_mixin();
    void test_model_generic_getter();
    void test_mixin_domains();
    void check_vector_api_approximate() { check_vector_api<ApproximateTag>(); }
    void check_vector_api_validated() { check_vector_api<ValidatedTag>(); }
    void check_vector_api_effective() { check_vector_api<EffectiveTag>(); }
    void test_generic_numeric_interfaces_approximate() { test_generic_numeric_interfaces<ApproximateTag>(); }
    void test_generic_numeric_interfaces_validated() { test_generic_numeric_interfaces<ValidatedTag>(); }
    void test_generic_numeric_interfaces_effective() { test_generic_numeric_interfaces<EffectiveTag>(); }
    void test_numeric_jets_and_taylor_arguments_approximate() { test_numeric_jets_and_taylor_arguments<ApproximateTag>(); }
    void test_numeric_jets_and_taylor_arguments_validated() { test_numeric_jets_and_taylor_arguments<ValidatedTag>(); }
    void test_numeric_jets_and_taylor_arguments_effective() { test_numeric_jets_and_taylor_arguments<EffectiveTag>(); }
    void test_handle_derivative_mixin_approximate() { test_handle_derivative_mixin<ApproximateTag>(); }
    void test_handle_derivative_mixin_validated() { test_handle_derivative_mixin<ValidatedTag>(); }
    void test_handle_derivative_mixin_effective() { test_handle_derivative_mixin<EffectiveTag>(); }
};

template<class F> concept HasVectorIndexing = requires(F& f,const F& cf) {
    cf.get(0u);f.set(0u,cf.create_zero());f[0u];cf[0u];
};
template<class P> Void TestFunctionInterface::check_vector_api() {
    using S=ScalarMultivariateFunction<P>;using V=VectorMultivariateFunction<P>;
    using U=ScalarUnivariateFunction<P>;using W=VectorUnivariateFunction<P>;
    static_assert(!HasVectorIndexing<S> && !HasVectorIndexing<U>);
    static_assert(HasVectorIndexing<V> && HasVectorIndexing<W>);
    static_assert(!Constructible<S,List<S>> && !Constructible<S,Vector<S>>);
    static_assert(!Constructible<U,List<U>> && !Constructible<U,Vector<U>>);
    static_assert(Constructible<V,List<S>> && Constructible<W,List<U>>);
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    V vector{x,y};const V& cv=vector;
    Vector<FloatDPApproximation> point({3,4},dp);
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (cv[Range(0u,2u)](point)); return snapshot.array(); }()), ([&] { auto snapshot = (point); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    ARIADNE_TESTING_EXPECT_EQUAL(cv[Range(1u,1u)](point).size(), 0u);
    ARIADNE_TESTING_EXPECT_EXCEPTION(vector.get(2u));
    ARIADNE_TESTING_EXPECT_EXCEPTION(cv[Range(0u,3u)]);
    ARIADNE_TESTING_EXPECT_EXCEPTION(S::coordinate(2u,2u));
    ARIADNE_TESTING_EXPECT_EXCEPTION(vector.set(2u,x));
    ARIADNE_TESTING_EXPECT_EXCEPTION(vector.set(0u,S::coordinate(1u,0u)));
    vector[0u]=vector[1u];
    ARIADNE_TESTING_EXPECT_SAME(vector[0u](point), point[1]);
    ARIADNE_TESTING_EXPECT_SAME(vector[0u].evaluate(point), point[1]);
    ARIADNE_TESTING_EXPECT_EQUAL(vector[0u].domain().dimension(), 2u);
    std::ostringstream text;text<<vector[0u];
    ARIADNE_TESTING_PROVE(!text.str().empty());
    auto composed=compose(V::identity(2u),vector);
    ARIADNE_TESTING_EXPECT_EXCEPTION(composed.set(0u,x));
    auto t=U::coordinate();W curve{t,t};
    auto composed_curve=compose(curve,t);
    ARIADNE_TESTING_EXPECT_EXCEPTION(composed_curve.set(0u,t));
    ARIADNE_TESTING_EXPECT_EXCEPTION(curve.get(2u));
    ARIADNE_TESTING_EXPECT_EXCEPTION(curve.set(2u,t));
}

// A valid directly implemented vector interface can evaluate and differentiate
// without implementing the optional VectorOfFunctionInterface getter.
class EvaluatedVectorFunction
    : public FunctionMixin<EvaluatedVectorFunction,EffectiveTag,RealVector(RealVector)>
{
    SharedPointer<const EffectiveVectorMultivariateFunction> _function;
  public:
    explicit EvaluatedVectorFunction(SharedPointer<const EffectiveVectorMultivariateFunction> f) : _function(f) { }
    SizeType argument_size() const override { return _function->argument_size(); }
    SizeType result_size() const override { return _function->result_size(); }
    template<class X> Vector<X> operator()(const Vector<X>& x) const { return (*_function)(x); }
    friend EvaluatedVectorFunction derivative(const EvaluatedVectorFunction& f,SizeType j) {
        return EvaluatedVectorFunction(std::make_shared<const EffectiveVectorMultivariateFunction>(derivative(*f._function,j))); }
    friend OutputStream& operator<<(OutputStream& os,const EvaluatedVectorFunction& f) {
        return os<<"EvaluatedVectorFunction("<<*f._function<<")"; }
};
class EvaluatedCurveFunction
    : public FunctionMixin<EvaluatedCurveFunction,EffectiveTag,RealVector(RealScalar)>
{
    SharedPointer<const EffectiveVectorUnivariateFunction> _function;
  public:
    explicit EvaluatedCurveFunction(SharedPointer<const EffectiveVectorUnivariateFunction> f) : _function(f) { }
    SizeOne argument_size() const override { return SizeOne(); }
    SizeType result_size() const override { return _function->result_size(); }
    template<class X> Vector<X> operator()(const X& x) const { return (*_function)(x); }
    friend EvaluatedCurveFunction derivative(const EvaluatedCurveFunction& f,IndexZero j) {
        return EvaluatedCurveFunction(std::make_shared<const EffectiveVectorUnivariateFunction>(derivative(*f._function,j))); }
    friend OutputStream& operator<<(OutputStream& os,const EvaluatedCurveFunction& f) {
        return os<<"EvaluatedCurveFunction("<<*f._function<<")"; }
};
Void TestFunctionInterface::test_optional_getter() {
    using S=EffectiveScalarMultivariateFunction;using V=EffectiveVectorMultivariateFunction;
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    V function(new EvaluatedVectorFunction(std::make_shared<const V>(V{x+y,x-y})));
    Vector<FloatDPApproximation> point({3,1},dp);
    ARIADNE_TESTING_EXPECT_SAME(function.get(0u)(point), FloatDPApproximation(4,dp));
    ARIADNE_TESTING_EXPECT_SAME(function.get(1u)(point), FloatDPApproximation(2,dp));
    ARIADNE_TESTING_EXPECT_SAME(function.get(1u).derivative(1u)(point), FloatDPApproximation(-1,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(function.get(0u).result_size(), 1u);
    auto check_field=[&](const auto& field) {
        ARIADNE_TESTING_EXPECT_SAME(field.get(0u)(point), FloatDPApproximation(4,dp));
        ARIADNE_TESTING_EXPECT_SAME(field.get(1u)(point), FloatDPApproximation(2,dp));
    };
    check_field(ValidatedVectorMultivariateFunction(function));
    check_field(ApproximateVectorMultivariateFunction(function));
    auto t=EffectiveScalarUnivariateFunction::coordinate();
    EffectiveVectorUnivariateFunction curve(new EvaluatedCurveFunction(
        std::make_shared<const EffectiveVectorUnivariateFunction>(EffectiveVectorUnivariateFunction{t,t+2})));
    FloatDPApproximation parameter(3,dp);
    auto check_curve=[&](const auto& input_curve) {
        ARIADNE_TESTING_EXPECT_SAME(input_curve.get(0u)(parameter), FloatDPApproximation(3,dp));
        ARIADNE_TESTING_EXPECT_SAME(input_curve.get(1u)(parameter), FloatDPApproximation(5,dp));
        ARIADNE_TESTING_EXPECT_SAME(input_curve.get(1u).derivative(IndexZero())(parameter), FloatDPApproximation(1,dp));
    };
    check_curve(curve);
    check_curve(ValidatedVectorUnivariateFunction(curve));
    check_curve(ApproximateVectorUnivariateFunction(curve));
    std::ostringstream text;text<<function.get(0u);
    ARIADNE_TESTING_PROVE(!text.str().empty());
}

class TaylorBackedFactory : public FunctionFactoryInterface<ValidatedTag> {
    TaylorFunctionFactory _factory;
  public:
    TaylorBackedFactory() : _factory(Sweeper<FloatDP>()) { }
    TaylorBackedFactory* clone() const override { return new TaylorBackedFactory(*this); }
    OutputStream& _write(OutputStream& os) const override { return os<<"TaylorBackedFactory"; }
  private:
    ValidatedScalarMultivariateFunction::Interface* _create(const BoxDomainType& domain,const ValidatedScalarMultivariateFunction::Interface& f) const override {
        auto model=_factory.create(domain,f);return new decltype(model)(model); }
    ValidatedVectorMultivariateFunction::Interface* _create(const BoxDomainType& domain,const ValidatedVectorMultivariateFunction::Interface& f) const override {
        auto model=_factory.create(domain,f);return new decltype(model)(model); }
};
Void TestFunctionInterface::test_function_factory() {
    TaylorBackedFactory backend;
    std::array<ValidatedFunctionFactory,3> factories{
        ValidatedFunctionFactory(backend),
        ValidatedFunctionFactory(new TaylorBackedFactory()),
        ValidatedFunctionFactory(std::make_shared<const TaylorBackedFactory>())};
    BoxDomainType domain{{-1,1},{-1,1}};
    auto source=ValidatedScalarMultivariateFunction::coordinate(2u,0u)+ValidatedScalarMultivariateFunction::coordinate(2u,1u);
    auto identity=ValidatedVectorMultivariateFunction::identity(2u);
    Vector<FloatDPBounds> point({0.25_x,0.5_x},dp);
    for(const auto& factory:factories) {
        ARIADNE_TESTING_PROVE(definitely(factory.create(domain,source.reference())(point)==FloatDPBounds(0.75_x,dp)));
        ARIADNE_TESTING_EXPECT_EQUAL(factory.create(domain,identity.reference())(point), point);
        ARIADNE_TESTING_PROVE(definitely(factory.create_zero(domain)(point)==FloatDPBounds(0,dp)));
        ARIADNE_TESTING_EXPECT_EQUAL(factory.create_identity(domain)(point), point);
        std::ostringstream text;text<<factory;
        ARIADNE_TESTING_EXPECT_EQUAL(text.str(), String("TaylorBackedFactory"));
    }
    std::unique_ptr<FunctionFactoryInterface<ValidatedTag>> copy(backend._copy());
    ARIADNE_TESTING_EXPECT_EQUAL(copy->create_identity(domain)(point), point);
    std::ostringstream text;text<<*copy;
    ARIADNE_TESTING_EXPECT_EQUAL(text.str(), String("TaylorBackedFactory"));
}

// Each numeric paradigm promises evaluation with its own numbers and weaker
// paradigms. Formula and type-erased algebra arguments must retain that same
// polynomial, rather than dropping the constant or changing its variables.
template<class Y,class P> Void TestFunctionInterface::check_generic_polynomial(const ScalarMultivariateFunction<P>& function) {
    using X=FloatDPApproximation;using F=Formula<Y>;using A=ElementaryAlgebra<Y>;
    Vector<Y> numbers{Y(2),Y(3)};
    ARIADNE_TESTING_EXPECT_SAME(X(function(numbers),dp), X(8,dp));
    Vector<X> point({2,3},dp);
    auto formula=function(F::identity(2u));
    ARIADNE_TESTING_EXPECT_SAME(evaluate(formula,point), X(8,dp));
    Vector<A> algebra_arguments{A(F::coordinate(0u)),A(F::coordinate(1u))};
    auto algebra_result=function(algebra_arguments);
    auto extracted=algebra_result.template extract<F>();
    ARIADNE_TESTING_EXPECT_SAME(evaluate(extracted,point), X(8,dp));
    ARIADNE_TESTING_EXPECT_SAME(evaluate(derivative(extracted,0u),point), X(4,dp));
}
template<class P> Void TestFunctionInterface::test_generic_numeric_interfaces() {
    using S=ScalarMultivariateFunction<P>;
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);auto polynomial=sqr(x)+y+Number<P>(1);
    check_generic_polynomial<ApproximateNumber>(polynomial);
    if constexpr (!Same<P,ApproximateTag>) { check_generic_polynomial<ValidatedNumber>(polynomial); }
    if constexpr (Same<P,EffectiveTag>) {
        check_generic_polynomial<EffectiveNumber>(polynomial);
        check_generic_polynomial<Real>(polynomial);
    }
}


static_assert(Same<typename TaylorModel<ApproximateTag,FloatDP>::Paradigm,ApproximateTag>);
static_assert(Same<typename TaylorModel<ValidatedTag,FloatDP>::Paradigm,ValidatedTag>);

// These are separate advertised virtual overloads. Jets must retain derivatives,
// and Taylor arguments must retain the same polynomial and its enclosure when
// the coefficient representation changes.
template<class MT,class P> Void TestFunctionInterface::check_taylor_argument(const ScalarMultivariateFunction<P>& function,typename MT::PrecisionType pr) {
    using F=typename MT::RawFloatType;
    Sweeper<F> sweeper=ThresholdSweeper<F>(pr,1e-12);
    auto variables=MT::coordinates(2u,sweeper);
    auto polynomial=function(variables);
    ARIADNE_TESTING_EXPECT_EQUAL(polynomial.argument_size(), 2u);
    ARIADNE_TESTING_EXPECT_EQUAL(polynomial.precision(), pr);
    if constexpr (Same<typename MT::Paradigm,ApproximateTag>) {
        Vector<Approximation<F>> point({0.25_x,0.5_x},pr);
        ARIADNE_TESTING_EXPECT_SAME(evaluate(polynomial,point), Approximation<F>(1.5625_x,pr));
        ARIADNE_TESTING_EXPECT_SAME(evaluate(derivative(polynomial,0u),point), Approximation<F>(0.5_x,pr));
    } else {
        Vector<Bounds<F>> point({0.25_x,0.5_x},pr);
        auto as_interval=[](const auto& value) {
            if constexpr (requires { value.lower_bound(); }) { return value; }
            else { return make_interval(value); }
        };
        auto actual=as_interval(evaluate(polynomial,point));
        ARIADNE_TESTING_PROVE(decide(actual.lower_bound()<=F(1.5625_x,pr)));
        ARIADNE_TESTING_PROVE(decide(actual.upper_bound()>=F(1.5625_x,pr)));
        auto gradient=as_interval(evaluate(derivative(polynomial,0u),point));
        ARIADNE_TESTING_PROVE(decide(gradient.lower_bound()<=F(0.5_x,pr)));
        ARIADNE_TESTING_PROVE(decide(gradient.upper_bound()>=F(0.5_x,pr)));
    }
}
template<class P> Void TestFunctionInterface::test_numeric_jets_and_taylor_arguments() {
    using S=ScalarMultivariateFunction<P>;
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    auto polynomial=sqr(x)+y+Number<P>(1);
    MultiplePrecision mp(128u);
    Vector<FloatMPApproximation> point({0.25_x,0.5_x},mp);
    ARIADNE_TESTING_EXPECT_SAME(polynomial(point), FloatMPApproximation(1.5625_x,mp));
    auto jet=polynomial(Differential<FloatMPApproximation>::variables(2u,point));
    ARIADNE_TESTING_EXPECT_SAME(jet.value(), FloatMPApproximation(1.5625_x,mp));
    ARIADNE_TESTING_EXPECT_SAME(jet.gradient()[0u], FloatMPApproximation(0.5_x,mp));
    ARIADNE_TESTING_EXPECT_SAME(jet.gradient()[1u], FloatMPApproximation(1,mp));
    ARIADNE_TESTING_EXPECT_SAME(jet.hessian().get(0u,0u), FloatMPApproximation(2,mp));
    check_taylor_argument<TaylorModel<ApproximateTag,FloatDP>>(polynomial,dp);
    check_taylor_argument<TaylorModel<ApproximateTag,FloatMP>>(polynomial,mp);
    if constexpr (!Same<P,ApproximateTag>) {
        Vector<FloatMPBounds> bounded_point({0.25_x,0.5_x},mp);
        ARIADNE_TESTING_PROVE(definitely(polynomial(bounded_point)==FloatMPBounds(1.5625_x,mp)));
        auto bounded_jet=polynomial(Differential<FloatMPBounds>::variables(2u,bounded_point));
        ARIADNE_TESTING_PROVE(definitely(bounded_jet.gradient()[0u]==FloatMPBounds(0.5_x,mp)));
        ARIADNE_TESTING_PROVE(definitely(bounded_jet.hessian().get(0u,0u)==FloatMPBounds(2,mp)));
        check_taylor_argument<TaylorModel<ValidatedTag,FloatMP>>(polynomial,mp);
        check_taylor_argument<TaylorModel<ValidatedTag,FloatDPBounds>>(polynomial,dp);
        check_taylor_argument<TaylorModel<ValidatedTag,FloatMPBounds>>(polynomial,mp);
        check_taylor_argument<TaylorModel<ValidatedTag,FloatDPUpperInterval>>(polynomial,dp);
        check_taylor_argument<TaylorModel<ValidatedTag,FloatMPUpperInterval>>(polynomial,mp);
    }
}

// A mixin-backed implementation may return its derivative as a public handle.
// Ownership conversion must work for handles as well as concrete interfaces.
template<class P> class HandleDifferentiableFunction
    : public FunctionMixin<HandleDifferentiableFunction<P>,P,Real(RealVector)>
{
    SharedPointer<const ScalarMultivariateFunction<P>> _function;
  public:
    explicit HandleDifferentiableFunction(SharedPointer<const ScalarMultivariateFunction<P>> f) : _function(f) { }
    SizeType argument_size() const override { return _function->argument_size(); }
    SizeOne result_size() const override { return SizeOne(); }
    template<class X> X operator()(const Vector<X>& values) const { return (*_function)(values); }
    friend ScalarMultivariateFunction<P> derivative(const HandleDifferentiableFunction& f,SizeType j) {
        return f._function->derivative(j); }
    friend OutputStream& operator<<(OutputStream& os,const HandleDifferentiableFunction& f) { return os<<*f._function; }
};

template<class P> Void TestFunctionInterface::test_handle_derivative_mixin() {
    using S=ScalarMultivariateFunction<P>;
    auto x=S::coordinate(2u,0u),y=S::coordinate(2u,1u);
    const S function(new HandleDifferentiableFunction<P>(std::make_shared<const S>(sqr(x)+y)));
    const Vector<FloatDPApproximation> point({3,1},dp);
    ARIADNE_TESTING_EXPECT_SAME(function.derivative(0u)(point), FloatDPApproximation(6,dp));
    const FunctionInterface<ApproximateTag,Real(RealVector)>& approximate=function.reference();
    std::unique_ptr<FunctionInterface<ApproximateTag,Real(RealVector)>> copy(approximate._copy());
    ARIADNE_TESTING_EXPECT_SAME(copy->_call(point), FloatDPApproximation(10,dp));
    if constexpr (!Same<P,ApproximateTag>) {
        const FunctionInterface<ValidatedTag,Real(RealVector)>& validated=function.reference();
        const MultiplePrecision mp(128u);
        ARIADNE_TESTING_EXPECT_SAME(validated._call(Vector<FloatMP>{FloatMP(3,mp),FloatMP(1,mp)}), FloatMPBounds(10,mp));
        const auto a=ValidatedScalarMultivariateFunction::coordinate(2u,0u);
        const auto b=ValidatedScalarMultivariateFunction::coordinate(2u,1u);
        const auto composed=validated._call(Vector<ValidatedScalarMultivariateFunction>{a+b,b});
        ARIADNE_TESTING_EXPECT_SAME(composed(point), FloatDPApproximation(17,dp));
    }
}
Void TestFunctionInterface::test_model_generic_getter() {
    const BoxDomainType domain{{-1,1},{-1,1}};
    const TaylorFunctionFactory backend(ThresholdSweeper<FloatDP>(dp,1e-12));
    const ValidatedVectorMultivariateFunction function=backend.create_identity(domain);
    const Vector<FloatDPBounds> point({0.25_x,0.5_x},dp);
    ARIADNE_TESTING_EXPECT_SAME(function.get(1u)(point), point[1u]);
}
Void TestFunctionInterface::test_mixin_domains() {
    const auto line=make_domain<IntervalDomainType>(1u);
    ARIADNE_TESTING_EXPECT_EQUAL(line, IntervalDomainType(-inf,+inf));
    ARIADNE_TESTING_EXPECT_EXCEPTION(make_domain<IntervalDomainType>(0u));
    ARIADNE_TESTING_EXPECT_EXCEPTION(make_domain<IntervalDomainType>(2u));
    ARIADNE_TESTING_EXPECT_EQUAL(make_domain<BoxDomainType>(0u).dimension(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(make_domain<BoxDomainType>(2u), BoxDomainType(2u,IntervalDomainType(-inf,+inf)));
}

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_function_interface");
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, check_vector_api_approximate);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, check_vector_api_validated);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, check_vector_api_effective);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_generic_numeric_interfaces_approximate);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_generic_numeric_interfaces_validated);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_generic_numeric_interfaces_effective);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_numeric_jets_and_taylor_arguments_approximate);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_numeric_jets_and_taylor_arguments_validated);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_numeric_jets_and_taylor_arguments_effective);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_handle_derivative_mixin_approximate);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_handle_derivative_mixin_validated);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_handle_derivative_mixin_effective);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_model_generic_getter);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_mixin_domains);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_optional_getter);
    ARIADNE_TESTING_ADD_TEST(TestFunctionInterface, test_function_factory);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
