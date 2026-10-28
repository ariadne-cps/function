/***************************************************************************
 *            test_function_models.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "macro/testing.hpp"
#include <sstream>
#include <memory>
#include <array>
using namespace Ariadne;
using ScalarModel=ScalarMultivariateFunctionModel<ValidatedTag,DoublePrecision>;
using VectorModel=VectorMultivariateFunctionModel<ValidatedTag,DoublePrecision>;
using ScalarTaylor=ValidatedScalarMultivariateTaylorFunctionModelDP;
using VectorTaylor=ValidatedVectorMultivariateTaylorFunctionModelDP;

Void test_generic_vector_derivative() {
    BoxDomainType domain{{2,6},{-2,2}};
    TaylorFunctionFactory backend(ThresholdSweeper<FloatDP>(dp,1e-12));
    auto x=backend.create_coordinate(domain,0u),y=backend.create_coordinate(domain,1u);
    VectorModel model=VectorTaylor{2*x+y,x-y};
    Vector<FloatDPBounds> point({3,1},dp);
    VectorModel differentiated(model.reference()._concrete_derivative(0u));
    ARIADNE_TEST_EQUALS(differentiated(point),Vector<FloatDPBounds>({2,1},dp));
}
Void test_generic_scalar_composition() {
    BoxDomainType domain{{2,6},{-2,2}};
    TaylorFunctionFactory backend(ThresholdSweeper<FloatDP>(dp,1e-12));
    ScalarModel model=backend.create_coordinate(domain,0u);
    auto t=ValidatedScalarUnivariateFunction::coordinate();
    ValidatedVectorUnivariateFunction curve{t+1,2*t};
    Vector<FloatDPBounds> point({3,1},dp);
    ScalarModel scalar(model.reference()._compose(t+1));
    VectorModel vector(model.reference()._compose(curve));
    ARIADNE_TEST_ASSERT(definitely(scalar(point)==FloatDPBounds(4,dp)));
    ARIADNE_TEST_EQUALS(vector(point),Vector<FloatDPBounds>({4,6},dp));
    ScalarModel unchecked_scalar(model.reference()._unchecked_compose(t+1));
    VectorModel unchecked_vector(model.reference()._unchecked_compose(curve));
    ARIADNE_TEST_ASSERT(definitely(unchecked_scalar(point)==FloatDPBounds(4,dp)));
    ARIADNE_TEST_EQUALS(unchecked_vector(point),Vector<FloatDPBounds>({4,6},dp));
}
Void test_erased_factory_projection() {
    FunctionModelFactory<ValidatedTag,DoublePrecision> backend(new TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-12)));
    BoxDomainType domain{{2,6},{-2,2}};
    auto projected=backend.create_projection(domain,Range(1u,2u));
    ARIADNE_TEST_EQUALS(projected(Vector<FloatDPBounds>({3,1},dp)),Vector<FloatDPBounds>({FloatDPBounds(1,dp)}));
}

Void test_concrete_scaled_calculus() {
    BoxDomainType domain{{2,6},{-2,2}};
    Sweeper<FloatDP> sweeper=ThresholdSweeper<FloatDP>(dp,1e-12);
    auto x=ScalarTaylor::coordinate(domain,0u,sweeper),y=ScalarTaylor::coordinate(domain,1u,sweeper);
    auto f=2*x+y;VectorTaylor vector{f,x-y};const VectorTaylor& cv=vector;
    Vector<FloatDPBounds> point({3,1},dp),reduced{FloatDPBounds(3,dp)};
    Vector<ValidatedNumber> generic_point{3,1};
    ARIADNE_TEST_ASSERT(definitely(f.create_constant(FloatDPBounds(5,dp))(point)==FloatDPBounds(5,dp)));
    ARIADNE_TEST_ASSERT(definitely(f.create_constant(ValidatedNumber(5))(point)==FloatDPBounds(5,dp)));
    ARIADNE_TEST_ASSERT(definitely(partial_evaluate(f,1u,FloatDPBounds(1,dp))(reduced)==FloatDPBounds(7,dp)));
    ARIADNE_TEST_ASSERT(definitely(partial_evaluate(f,1u,ValidatedNumber(1))(reduced)==FloatDPBounds(7,dp)));
    ARIADNE_TEST_EQUALS(partial_evaluate(vector,1u,FloatDPBounds(1,dp))(reduced),Vector<FloatDPBounds>({7,2},dp));
    ARIADNE_TEST_EQUALS(partial_evaluate(vector,1u,ValidatedNumber(1))(reduced),Vector<FloatDPBounds>({7,2},dp));
    ARIADNE_TEST_EQUALS(derivative(vector,0u)(point),Vector<FloatDPBounds>({2,1},dp));
    ARIADNE_TEST_EQUALS(antiderivative(vector,0u,FloatDPBounds(2,dp))(point),Vector<FloatDPBounds>({6,1.5_x},dp));
    ARIADNE_TEST_EQUALS(antiderivative(vector,0u,ValidatedNumber(2))(point),Vector<FloatDPBounds>({6,1.5_x},dp));
    ARIADNE_TEST_EQUALS(unchecked_evaluate(vector,point),Vector<FloatDPBounds>({7,2},dp));
    ARIADNE_TEST_EQUALS(unchecked_evaluate(vector,generic_point),Vector<FloatDPBounds>({7,2},dp));
    ARIADNE_TEST_ASSERT(definitely(evaluate(f,generic_point)==FloatDPBounds(7,dp)));
    ARIADNE_TEST_EQUALS(cv.values(),Vector<FloatDP>({8,4},dp));
    ARIADNE_TEST_EQUALS(cv.expansions().size(),2u);
    ARIADNE_TEST_EQUALS(cv.errors()[0u].raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(cv.get(1u)(point),cv[1u](point));
    BoxDomainType prefix{{-1,1}},suffix{{0,1}};
    Vector<FloatDPBounds> enlarged({0,3,1,0.5_x},dp);
    ARIADNE_TEST_ASSERT(definitely(embed(prefix,f,suffix)(enlarged)==FloatDPBounds(7,dp)));
    ARIADNE_TEST_EQUALS(embed(prefix,vector,suffix)(enlarged),Vector<FloatDPBounds>({7,2},dp));
    auto projected=ScalarTaylor::projection(domain,Range(1u,2u),sweeper);
    ARIADNE_TEST_EQUALS(projected(point),Vector<FloatDPBounds>{FloatDPBounds(1,dp)});
    VectorTaylor empty(0u,domain,sweeper);
    ARIADNE_TEST_EQUALS(empty.values().size(),0u);
    ARIADNE_TEST_EQUALS(empty.expansions().size(),0u);
    ARIADNE_TEST_EQUALS(derivative(empty,0u)(point).size(),0u);
    ARIADNE_TEST_EQUALS(antiderivative(empty,0u,FloatDPBounds(2,dp))(point).size(),0u);
    std::ostringstream text;text<<representation(f)<<representation(vector)<<representation(empty);
    auto constant=ScalarTaylor::constant(BoxDomainType(0u),FloatDPBounds(3,dp),sweeper);
    text<<representation(constant);
    ARIADNE_TEST_ASSERT(text.str().find("ScaledFunctionPatch")!=String::npos);
    ARIADNE_TEST_ASSERT(text.str().find("VectorScaledFunctionPatch")!=String::npos);
}
Void test_scaled_enclosure_comparisons() {
    BoxDomainType domain{{2,6},{-2,2}},small{{3,5},{-1,1}},other{{4,6},{-1,1}};
    Sweeper<FloatDP> sweeper=ThresholdSweeper<FloatDP>(dp,1e-12);
    auto broad=ScalarTaylor::constant(domain,FloatDPBounds(1,3,dp),sweeper);
    auto exact=ScalarTaylor::constant(domain,FloatDPBounds(2,dp),sweeper);
    auto far=ScalarTaylor::constant(domain,FloatDPBounds(7,dp),sweeper);
    Vector<FloatDPBounds> point({4,0},dp);
    ARIADNE_TEST_ASSERT(refines(exact,broad));
    ARIADNE_TEST_ASSERT(!refines(broad,exact));
    ARIADNE_TEST_ASSERT(inconsistent(exact,far));
    ARIADNE_TEST_ASSERT(!inconsistent(exact,broad));
    ARIADNE_TEST_ASSERT(refines(exact,restriction(broad,small)));
    ARIADNE_TEST_ASSERT(!refines(restriction(exact,small),broad));
    ARIADNE_TEST_ASSERT(inconsistent(restriction(exact,small),restriction(far,other)));
    ARIADNE_TEST_ASSERT(definitely(refinement(exact,broad)(point)==FloatDPBounds(2,dp)));
    VectorTaylor ve{exact,exact},vb{broad,broad},vf{exact,far};
    ARIADNE_TEST_ASSERT(refines(ve,vb));ARIADNE_TEST_ASSERT(!refines(vb,ve));
    ARIADNE_TEST_ASSERT(inconsistent(ve,vf));ARIADNE_TEST_ASSERT(!inconsistent(ve,vb));
    ARIADNE_TEST_EQUALS(refinement(ve,vb)(point),Vector<FloatDPBounds>({2,2},dp));
    ValidatedNumber scalar(5);
    ARIADNE_TEST_ASSERT(definitely(max(exact,scalar)(point)==FloatDPBounds(5,dp)));
    ARIADNE_TEST_ASSERT(definitely(min(exact,scalar)(point)==FloatDPBounds(2,dp)));
    ARIADNE_TEST_ASSERT(definitely(max(scalar,exact)(point)==FloatDPBounds(5,dp)));
    ARIADNE_TEST_ASSERT(definitely(min(scalar,exact)(point)==FloatDPBounds(2,dp)));
    using Operations=AlgebraOperations<ScalarTaylor>;
    ARIADNE_TEST_ASSERT(definitely(Operations::apply(BinaryElementaryOperator(Sub()),FloatDPBounds(5,dp),exact)(point)==FloatDPBounds(3,dp)));
}
Void test_model_handles_and_factories() {
    BoxDomainType domain{{2,6},{-2,2}},small{{3,5},{-1,1}};
    Sweeper<FloatDP> sweeper=ThresholdSweeper<FloatDP>(dp,1e-12);
    TaylorFunctionFactory concrete(sweeper);
    FunctionModelFactory<ValidatedTag,DoublePrecision> backend(new TaylorFunctionFactory(sweeper));
    auto x=concrete.create_coordinate(domain,0u),y=concrete.create_coordinate(domain,1u);
    ScalarModel scalar=2*x+y;
    VectorModel vector=VectorTaylor{2*x+y,x-y};const VectorModel& cv=vector;
    Vector<FloatDPBounds> point({3,1},dp),reduced{FloatDPBounds(3,dp)};
    Vector<ValidatedNumber> generic_point{3,1};
    ARIADNE_TEST_ASSERT(definitely(unchecked_evaluate(scalar,point)==FloatDPBounds(7,dp)));
    ARIADNE_TEST_ASSERT(definitely(unchecked_evaluate(scalar,generic_point)==ValidatedNumber(7)));
    ARIADNE_TEST_EQUALS(unchecked_evaluate(vector,point),Vector<FloatDPBounds>({7,2},dp));
    ARIADNE_TEST_ASSERT(definitely(partial_evaluate(scalar,1u,FloatDPBounds(1,dp))(reduced)==FloatDPBounds(7,dp)));
    ARIADNE_TEST_ASSERT(definitely(partial_evaluate(scalar,1u,ValidatedNumber(1))(reduced)==FloatDPBounds(7,dp)));
    ARIADNE_TEST_EQUALS(vector.values(),Vector<FloatDP>({8,4},dp));
    ARIADNE_TEST_EQUALS(vector.errors()[1u].raw(),FloatDP(0,dp));
    ARIADNE_TEST_ASSERT(decide(norm(vector).raw()>=FloatDP(7,dp)));
    ARIADNE_TEST_EQUALS(cv[Range(1u,2u)](point),Vector<FloatDPBounds>{FloatDPBounds(2,dp)});
    auto projected=backend.create_projection(domain,Range(0u,2u));
    ARIADNE_TEST_EQUALS(projected(point),point);
    ARIADNE_TEST_ASSERT(definitely(backend.create_coordinate(domain,1u)(point)==FloatDPBounds(1,dp)));
    ARIADNE_TEST_ASSERT(definitely(backend.create_zero(domain)(point)==FloatDPBounds(0,dp)));
    ARIADNE_TEST_EQUALS(backend.create_zeros(0u,domain)(point).size(),0u);
    ARIADNE_TEST_EQUALS(backend.create_constants(domain,Vector<ValidatedNumber>{2,3})(point),Vector<FloatDPBounds>({2,3},dp));
    auto creator=factory(scalar);
    ARIADNE_TEST_ASSERT(definitely(creator.create_constant(ValidatedNumber(5))(point)==FloatDPBounds(5,dp)));
    ARIADNE_TEST_ASSERT(definitely(creator.create_zero()(point)==FloatDPBounds(0,dp)));
    ARIADNE_TEST_EQUALS(creator.create_zeros(0u)(point).size(),0u);
    ARIADNE_TEST_EQUALS(creator.create_identity()(point),point);
    auto restricted=restriction(vector,small);auto restricted_scalar=restriction(scalar,small);
    Vector<FloatDPBounds> centre({4,0},dp);
    ARIADNE_TEST_EQUALS(restricted(centre),Vector<FloatDPBounds>({8,4},dp));
    ARIADNE_TEST_ASSERT(definitely(restricted_scalar(centre)==FloatDPBounds(8,dp)));
    ARIADNE_TEST_EQUALS(antiderivative(vector,0u,FloatDPBounds(2,dp))(point),Vector<FloatDPBounds>({6,1.5_x},dp));
    ARIADNE_TEST_ASSERT(definitely(antiderivative(scalar,0u,FloatDPBounds(2,dp))(point)==FloatDPBounds(6,dp)));
    ScalarModel copy=scalar;copy.clobber();
    ARIADNE_TEST_EQUALS(copy.error().raw(),FloatDP(0,dp));
    VectorModel clone(vector.reference()._concrete_clone());clone.clobber();
    ARIADNE_TEST_EQUALS(clone(point),vector(point));
    std::ostringstream text;text<<backend<<scalar.reference()<<vector.reference();
    ARIADNE_TEST_ASSERT(!text.str().empty());
}

Void test_model_factory_interfaces() {
    using FactoryInterface=FunctionModelFactoryInterface<ValidatedTag,DoublePrecision>;
    const BoxDomainType domain{{-1,1},{-1,1}};
    const Vector<FloatDPBounds> point({0.25_x,0.5_x},dp);
    const auto x=ValidatedScalarMultivariateFunction::coordinate(2u,0u);
    const auto y=ValidatedScalarMultivariateFunction::coordinate(2u,1u);
    const ValidatedVectorMultivariateFunction identity{x,y};
    const Sweeper<FloatDP> sweeper=ThresholdSweeper<FloatDP>(dp,1e-12);
    std::array<std::unique_ptr<FactoryInterface>,2> factories{
        std::unique_ptr<FactoryInterface>(make_taylor_function_factory()),
        std::unique_ptr<FactoryInterface>(make_taylor_function_factory(sweeper))};
    for(const auto& ptr:factories) {
        const auto& backend=*ptr;
        ARIADNE_TEST_SAME(backend.create(ValidatedNumber(3)),FloatDPBounds(3,dp));
        ARIADNE_TEST_SAME(backend.create_number(ValidatedNumber(3)),FloatDPBounds(3,dp));
        ARIADNE_TEST_SAME(backend.create(domain,(x+y).reference())(point),FloatDPBounds(0.75_x,dp));
        ARIADNE_TEST_EQUALS(backend.create(domain,identity.reference())(point),point);
        ARIADNE_TEST_SAME(backend.create_zero(domain)(point),FloatDPBounds(0,dp));
        ARIADNE_TEST_SAME(backend.create_constant(domain,ValidatedNumber(4))(point),FloatDPBounds(4,dp));
        ARIADNE_TEST_SAME(backend.create_constant(domain,FloatDPBounds(5,dp))(point),FloatDPBounds(5,dp));
        ARIADNE_TEST_SAME(backend.create_coordinate(domain,1u)(point),point[1u]);
        ARIADNE_TEST_EQUALS(backend.create_zeros(2u,domain)(point),Vector<FloatDPBounds>({0,0},dp));
        ARIADNE_TEST_EQUALS(backend.create_constants(domain,Vector<ValidatedNumber>{2,3})(point),Vector<FloatDPBounds>({2,3},dp));
        ARIADNE_TEST_EQUALS(backend.create_projection(domain,Range(1u,2u))(point),Vector<FloatDPBounds>{point[1u]});
        ARIADNE_TEST_EQUALS(backend.create_identity(domain)(point),point);
        const Vector<FloatDPBounds> parameter{point[0u]};
        ARIADNE_TEST_SAME(backend.create_identity(IntervalDomainType(-1,1))(parameter),point[0u]);
        std::unique_ptr<FactoryInterface> copied(backend._copy());
        ARIADNE_TEST_EQUALS(copied->create_identity(domain)(point),point);
    }
    std::array<FunctionPatchFactory<ValidatedTag>,2> patch_factories{
        FunctionPatchFactory<ValidatedTag>(make_taylor_function_patch_factory()),
        FunctionPatchFactory<ValidatedTag>(make_taylor_function_patch_factory(sweeper))};
    for(const auto& backend:patch_factories) {
        ARIADNE_TEST_SAME(backend.create_coordinate(domain,1u)(point),point[1u]);
    }
    FunctionModelFactory<ValidatedTag,DoublePrecision> backend(make_taylor_function_factory(sweeper));
    auto creator=factory(backend.create_coordinate(domain,0u));
    const FloatDPBounds coefficient(7,dp);
    ARIADNE_TEST_ASSERT(&creator.create(coefficient)==&coefficient);
}

Void test_scalar_model_interface_delegation() {
    const BoxDomainType domain{{2,6},{-2,2}},prefix{{-1,1}},suffix{{0,1}};
    const Sweeper<FloatDP> sweeper=ThresholdSweeper<FloatDP>(dp,1e-12);
    const auto x=ScalarTaylor::coordinate(domain,0u,sweeper),y=ScalarTaylor::coordinate(domain,1u,sweeper);
    const ScalarModel model=2*x+y;
    const auto& interface=model.reference();
    const ElementaryAlgebraInterface<FloatDPBounds>& canonical=interface;
    const ElementaryAlgebraInterface<ValidatedNumber>& generic=interface;
    const auto from_algebra=[](auto* algebra) {
        return ScalarModel(dynamic_cast<ScalarModel::Interface*>(algebra));
    };
    const Vector<FloatDPBounds> point({3,1},dp),reduced{FloatDPBounds(3,dp)},enlarged({0,3,1,0.5_x},dp);
    ARIADNE_TEST_SAME(from_algebra(canonical._create_copy())(point),FloatDPBounds(7,dp));
    ARIADNE_TEST_SAME(from_algebra(canonical._create_zero())(point),FloatDPBounds(0,dp));
    ARIADNE_TEST_SAME(from_algebra(generic._create_constant(ValidatedNumber(5)))(point),FloatDPBounds(5,dp));
    ARIADNE_TEST_SAME(from_algebra(canonical._create_constant(FloatDPBounds(6,dp)))(point),FloatDPBounds(6,dp));
    ARIADNE_TEST_SAME(from_algebra(canonical._apply(BinaryElementaryOperator(Add()),FloatDPBounds(5,dp)))(point),FloatDPBounds(12,dp));
    ARIADNE_TEST_SAME(from_algebra(canonical._rapply(BinaryElementaryOperator(Sub()),FloatDPBounds(5,dp)))(point),FloatDPBounds(-2,dp));
    ARIADNE_TEST_SAME(from_algebra(generic._apply(BinaryElementaryOperator(Add()),ValidatedNumber(5)))(point),FloatDPBounds(12,dp));
    ARIADNE_TEST_EQUALS(interface._errors().get(dp).raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(interface._error().get(dp).raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(interface._norm().get(dp).raw(),FloatDP(14,dp));
    ARIADNE_TEST_EQUALS(interface._concrete_norm().raw(),FloatDP(14,dp));
    const ScalarModel primitive(interface._antiderivative(0u));
    ARIADNE_TEST_SAME(derivative(primitive,0u)(point),FloatDPBounds(7,dp));
    ARIADNE_TEST_SAME(ScalarModel(interface._concrete_antiderivative(0u,FloatDPBounds(2,dp)))(point),FloatDPBounds(6,dp));
    const MultiplePrecision mp(128u);
    ARIADNE_TEST_SAME(interface._unchecked_evaluate(Vector<FloatMPBounds>({3,1},mp)),FloatMPBounds(7,mp));
    ARIADNE_TEST_SAME(interface._concrete_unchecked_evaluate(point),FloatDPBounds(7,dp));
    ARIADNE_TEST_SAME(from_algebra(interface._partial_evaluate(1u,ValidatedNumber(1)))(reduced),FloatDPBounds(7,dp));
    ARIADNE_TEST_SAME(ScalarModel(interface._embed(prefix,suffix))(enlarged),FloatDPBounds(7,dp));
    const ScalarModel broad=ScalarTaylor::constant(domain,FloatDPBounds(1,3,dp),sweeper);
    const ScalarModel exact=ScalarTaylor::constant(domain,FloatDPBounds(2,dp),sweeper);
    const ScalarModel far=ScalarTaylor::constant(domain,FloatDPBounds(7,dp),sweeper);
    ARIADNE_TEST_ASSERT(exact.reference()._refines(broad.reference()));
    ARIADNE_TEST_ASSERT(decide(exact.reference()._concrete_refines(broad.reference())));
    ARIADNE_TEST_ASSERT(decide(exact.reference()._concrete_inconsistent(far.reference())));
    ARIADNE_TEST_SAME(ScalarModel(exact.reference()._concrete_refinement(broad.reference()))(point),FloatDPBounds(2,dp));
    ScalarTaylor uncertain=2*x+y;uncertain.set_error(FloatDPError(0.125_x,dp));
    ScalarModel cleared=uncertain;cleared.reference()._clobber();
    ARIADNE_TEST_EQUALS(cleared.error().raw(),FloatDP(0,dp));
    ARIADNE_TEST_SAME(cleared(point),FloatDPBounds(7,dp));
}

Void test_vector_model_interface_delegation() {
    const BoxDomainType domain{{2,6},{-2,2}},prefix{{-1,1}},suffix{{0,1}};
    const Sweeper<FloatDP> sweeper=ThresholdSweeper<FloatDP>(dp,1e-12);
    const auto x=ScalarTaylor::coordinate(domain,0u,sweeper),y=ScalarTaylor::coordinate(domain,1u,sweeper);
    const VectorModel model=VectorTaylor{2*x+y,x-y};
    const auto& interface=model.reference();
    const auto from_patch=[](auto* patch) {
        return VectorModel(dynamic_cast<VectorModel::Interface*>(patch));
    };
    const Vector<FloatDPBounds> point({3,1},dp),expected({7,2},dp),reduced{FloatDPBounds(3,dp)},enlarged({0,3,1,0.5_x},dp);
    ARIADNE_TEST_EQUALS(interface._range().dimension(),2u);
    ARIADNE_TEST_EQUALS(interface._errors()[0u].get(dp).raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(interface._error().get(dp).raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(interface._concrete_error().raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(interface._norm().get(dp).raw(),FloatDP(14,dp));
    const VectorModel primitive(interface._antiderivative(0u));
    ARIADNE_TEST_EQUALS(VectorModel(primitive.reference()._concrete_derivative(0u))(point),expected);
    ARIADNE_TEST_EQUALS(VectorModel(interface._antiderivative(0u,ValidatedNumber(2)))(point),Vector<FloatDPBounds>({6,1.5_x},dp));
    ARIADNE_TEST_EQUALS(VectorModel(interface._concrete_antiderivative(0u,FloatDPBounds(2,dp)))(point),Vector<FloatDPBounds>({6,1.5_x},dp));
    ARIADNE_TEST_EQUALS(VectorModel(interface._embed(prefix,suffix))(enlarged),expected);
    VectorModel appended(interface._concrete_clone());
    appended.reference()._concrete_adjoin(ScalarModel(x).reference());
    ARIADNE_TEST_EQUALS(appended(point),Vector<FloatDPBounds>({7,2,3},dp));
    ARIADNE_TEST_EQUALS(VectorModel(interface._concrete_join(interface))(point),Vector<FloatDPBounds>({7,2,7,2},dp));
    const Vector<FloatDPBounds> doubled_point({3,1,3,1},dp);
    ARIADNE_TEST_EQUALS(from_patch(interface._combine(interface))(doubled_point),Vector<FloatDPBounds>({7,2,7,2},dp));
    ARIADNE_TEST_EQUALS(VectorModel(interface._concrete_combine(interface))(doubled_point),Vector<FloatDPBounds>({7,2,7,2},dp));
    ARIADNE_TEST_EQUALS(Vector<FloatDPBounds>(interface._unchecked_evaluate(Vector<ValidatedNumber>{3,1}),dp),expected);
    const MultiplePrecision mp(128u);
    ARIADNE_TEST_EQUALS(interface._unchecked_evaluate(Vector<FloatMPBounds>({3,1},mp)),Vector<FloatMPBounds>({7,2},mp));
    ARIADNE_TEST_EQUALS(interface._concrete_unchecked_evaluate(point),expected);
    const auto a=ValidatedScalarMultivariateFunction::coordinate(2u,0u),b=ValidatedScalarMultivariateFunction::coordinate(2u,1u);
    ARIADNE_TEST_SAME(ScalarModel(interface._unchecked_compose(a+b))(point),FloatDPBounds(9,dp));
    ARIADNE_TEST_EQUALS(VectorModel(interface._unchecked_compose(ValidatedVectorMultivariateFunction{a,b}))(point),expected);
    ARIADNE_TEST_EQUALS(from_patch(interface._partial_evaluate(1u,ValidatedNumber(1)))(reduced),expected);
    ARIADNE_TEST_EQUALS(VectorModel(interface._concrete_partial_evaluate(1u,FloatDPBounds(1,dp)))(reduced),expected);
    const auto broad_scalar=ScalarTaylor::constant(domain,FloatDPBounds(1,3,dp),sweeper);
    const auto exact_scalar=ScalarTaylor::constant(domain,FloatDPBounds(2,dp),sweeper);
    const auto far_scalar=ScalarTaylor::constant(domain,FloatDPBounds(7,dp),sweeper);
    const VectorModel broad=VectorTaylor{broad_scalar,broad_scalar},exact=VectorTaylor{exact_scalar,exact_scalar},far=VectorTaylor{exact_scalar,far_scalar};
    ARIADNE_TEST_ASSERT(exact.reference()._refines(broad.reference()));
    ARIADNE_TEST_ASSERT(decide(exact.reference()._concrete_refines(broad.reference())));
    ARIADNE_TEST_ASSERT(decide(exact.reference()._concrete_inconsistent(far.reference())));
    ARIADNE_TEST_EQUALS(VectorModel(exact.reference()._concrete_refinement(broad.reference()))(point),Vector<FloatDPBounds>({2,2},dp));
    ScalarTaylor uncertain=2*x+y;uncertain.set_error(FloatDPError(0.125_x,dp));
    VectorModel cleared=VectorTaylor{uncertain,x-y};cleared.reference()._clobber();
    ARIADNE_TEST_EQUALS(cleared.error().raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(cleared(point),expected);
}

Int main() {
    ARIADNE_TEST_CALL(test_scalar_model_interface_delegation());
    ARIADNE_TEST_CALL(test_vector_model_interface_delegation());
    ARIADNE_TEST_CALL(test_model_factory_interfaces());
    ARIADNE_TEST_CALL(test_erased_factory_projection());
    ARIADNE_TEST_CALL(test_concrete_scaled_calculus());
    ARIADNE_TEST_CALL(test_scaled_enclosure_comparisons());
    ARIADNE_TEST_CALL(test_model_handles_and_factories());
    ARIADNE_TEST_CALL(test_generic_vector_derivative());
    ARIADNE_TEST_CALL(test_generic_scalar_composition());
    return ARIADNE_TEST_FAILURES;
}
