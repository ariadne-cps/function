/***************************************************************************
 *            test_function_patch_operations.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/function_patch.hpp"
#include "macro/testing.hpp"
#include <array>
#include <sstream>
using namespace Ariadne;

Void test_patch_constant_dispatch() {
    BoxDomainType domain{{-1,1},{-1,1}};
    TaylorFunctionFactory factory(ThresholdSweeper<FloatDP>(dp,1e-12));
    auto model=factory.create_constant(domain,ValidatedNumber(2));
    ValidatedScalarMultivariateFunction f=model;
    Vector<FloatDPBounds> point({0,0},dp);
    ARIADNE_TEST_ASSERT(definitely((f+ValidatedNumber(3))(point)==FloatDPBounds(5,dp)));
    using A=AlgebraOperations<ValidatedScalarMultivariateFunction,ValidatedNumber>;
    ARIADNE_TEST_ASSERT(definitely(A::apply(BinaryElementaryOperator(Sub()),ValidatedNumber(3),f)(point)==FloatDPBounds(1,dp)));
    ARIADNE_TEST_ASSERT(definitely((f-ValidatedNumber(3))(point)==FloatDPBounds(-1,dp)));
    ARIADNE_TEST_ASSERT(definitely((f*ValidatedNumber(3))(point)==FloatDPBounds(6,dp)));
    ARIADNE_TEST_ASSERT(definitely((f/ValidatedNumber(2))(point)==FloatDPBounds(1,dp)));
}
Void test_patch_vector_product() {
    BoxDomainType domain{{-1,1},{-1,1}};
    TaylorFunctionFactory factory(ThresholdSweeper<FloatDP>(dp,1e-12));
    ValidatedScalarMultivariateFunction scalar=factory.create_constant(domain,ValidatedNumber(2));
    ValidatedVectorMultivariateFunction vector=factory.create_constants(domain,Vector<ValidatedNumber>{4,8});
    Vector<FloatDPBounds> point({0,0},dp);
    ARIADNE_TEST_EQUALS((scalar*vector)(point),Vector<FloatDPBounds>({8,16},dp));
}
Void test_mixed_patch_operations() {
    using S=ValidatedScalarMultivariateFunction;
    using V=ValidatedVectorMultivariateFunction;
    BoxDomainType domain{{-1,1},{-1,1}};
    TaylorFunctionFactory factory(ThresholdSweeper<FloatDP>(dp,1e-12));
    std::array<S,2> scalar_a{S::constant(2u,ValidatedNumber(2)),factory.create_constant(domain,ValidatedNumber(2))};
    std::array<S,2> scalar_b{S::constant(2u,ValidatedNumber(4)),factory.create_constant(domain,ValidatedNumber(4))};
    Vector<ValidatedNumber> va{2,6},vb{4,8};
    std::array<V,2> vector_a{V::constant(2u,va),factory.create_constants(domain,va)};
    std::array<V,2> vector_b{V::constant(2u,vb),factory.create_constants(domain,vb)};
    Vector<FloatDPBounds> point({0,0},dp);
    for(SizeType i=0;i!=2u;++i) {
        const auto& a=scalar_a[i];
        ARIADNE_TEST_ASSERT(definitely((-a)(point)==FloatDPBounds(-2,dp)));
        ARIADNE_TEST_ASSERT(definitely(pow(a,2)(point)==FloatDPBounds(4,dp)));
        ARIADNE_TEST_ASSERT(definitely((ValidatedNumber(3)+a)(point)==FloatDPBounds(5,dp)));
        ARIADNE_TEST_ASSERT(definitely((ValidatedNumber(3)-a)(point)==FloatDPBounds(1,dp)));
        for(SizeType j=0;j!=2u;++j) {
            const auto& b=scalar_b[j];const auto& av=vector_a[i];const auto& bv=vector_b[j];
            ARIADNE_TEST_ASSERT(definitely((a+b)(point)==FloatDPBounds(6,dp)));
            ARIADNE_TEST_ASSERT(definitely((a-b)(point)==FloatDPBounds(-2,dp)));
            ARIADNE_TEST_ASSERT(definitely((a*b)(point)==FloatDPBounds(8,dp)));
            ARIADNE_TEST_ASSERT(definitely((a/b)(point)==FloatDPBounds(0.5_x,dp)));
            ARIADNE_TEST_EQUALS((av+bv)(point),Vector<FloatDPBounds>({6,14},dp));
            ARIADNE_TEST_EQUALS((av-bv)(point),Vector<FloatDPBounds>({-2,-2},dp));
            ARIADNE_TEST_EQUALS((a*bv)(point),Vector<FloatDPBounds>({8,16},dp));
            ARIADNE_TEST_EQUALS((bv*a)(point),Vector<FloatDPBounds>({8,16},dp));
            ARIADNE_TEST_EQUALS((bv/a)(point),Vector<FloatDPBounds>({2,4},dp));
            ARIADNE_TEST_EQUALS(join(a,b)(point),Vector<FloatDPBounds>({2,4},dp));
            ARIADNE_TEST_EQUALS(join(a,bv)(point),Vector<FloatDPBounds>({2,4,8},dp));
            ARIADNE_TEST_EQUALS(join(av,b)(point),Vector<FloatDPBounds>({2,6,4},dp));
            ARIADNE_TEST_EQUALS(join(av,bv)(point),Vector<FloatDPBounds>({2,6,4,8},dp));
        }
    }
    ValidatedScalarMultivariateFunctionPatch patch=factory.create_constant(domain,ValidatedNumber(2));
    ValidatedScalarMultivariateFunctionModelDP model=factory.create_constant(domain,ValidatedNumber(2));
    S scalar_model=model;
    V one_component{scalar_model};
    ARIADNE_TEST_EQUALS(one_component(point),Vector<FloatDPBounds>{FloatDPBounds(2,dp)});
    std::ostringstream text;text<<patch<<model<<factory;
    ARIADNE_TEST_ASSERT(!text.str().empty());
}
Void test_patch_handle_construction_and_values() {
    using SP=ValidatedScalarMultivariateFunctionPatch;
    using VP=ValidatedVectorMultivariateFunctionPatch;
    const BoxDomainType domain{{2,6},{-2,2}};
    const TaylorFunctionFactory backend(ThresholdSweeper<FloatDP>(dp,1e-12));
    const auto x=backend.create_coordinate(domain,0u),y=backend.create_coordinate(domain,1u);
    const SP original=2*x+y,second=x-y;
    const Vector<FloatDPBounds> point({3,1},dp);
    const Vector<ValidatedNumber> generic_point{3,1};
    SP null_scalar,copied_null(null_scalar),assigned(original);
    ARIADNE_TEST_ASSERT(copied_null.raw_pointer()==nullptr);
    assigned=null_scalar;ARIADNE_TEST_ASSERT(assigned.raw_pointer()==nullptr);
    assigned=original;ARIADNE_TEST_SAME(assigned(point),FloatDPBounds(7,dp));
    SP copy(original);ARIADNE_TEST_SAME(copy(point),FloatDPBounds(7,dp));
    SP from_function(static_cast<ValidatedScalarMultivariateFunction>(original));
    ARIADNE_TEST_SAME(from_function(point),FloatDPBounds(7,dp));
    ARIADNE_TEST_FAIL(SP(ValidatedScalarMultivariateFunction::coordinate(2u,0u)));
    ARIADNE_TEST_EQUALS(original.result_size(),1u);
    ARIADNE_TEST_EQUALS(original.argument_size(),2u);
    ARIADNE_TEST_ASSERT(same(original.codomain(),IntervalDomainType(2,14)));
    ARIADNE_TEST_EQUALS(original.error().get(dp).raw(),FloatDP(0,dp));
    ARIADNE_TEST_ASSERT(original.managed_pointer()!=nullptr);
    ARIADNE_TEST_SAME(derivative(original,0u)(point),FloatDPBounds(2,dp));
    assigned=ValidatedNumber(5);ARIADNE_TEST_SAME(assigned(point),FloatDPBounds(5,dp));
    assigned=ValidatedScalarMultivariateFunction::coordinate(2u,0u);
    ARIADNE_TEST_SAME(assigned(point),FloatDPBounds(3,dp));
    assigned=original.reference();ARIADNE_TEST_SAME(assigned(point),FloatDPBounds(7,dp));
    SP::Interface& mutable_interface=assigned;
    ARIADNE_TEST_EQUALS(mutable_interface.argument_size(),2u);
    ARIADNE_TEST_EQUALS(assigned.reference().argument_size(),2u);
    using Operations=AlgebraOperations<SP,ValidatedNumber>;
    ARIADNE_TEST_SAME(Operations::apply(BinaryElementaryOperator(Sub()),ValidatedNumber(5),original)(point),FloatDPBounds(-2,dp));
    SP uncertain=backend.create_constant(domain,ValidatedNumber(FloatDPBounds(1,3,dp)));
    uncertain.clobber();ARIADNE_TEST_EQUALS(uncertain.error().get(dp).raw(),FloatDP(0,dp));
    ARIADNE_TEST_SAME(uncertain(point),FloatDPBounds(2,dp));
    const VP vector(Array<SP>{original,second});
    VP listed(List<SP>{original,second});
    ARIADNE_TEST_EQUALS(vector(point),Vector<FloatDPBounds>({7,2},dp));
    ARIADNE_TEST_EQUALS(listed(point),vector(point));
    ARIADNE_TEST_EQUALS(listed.reference().result_size(),2u);
    ARIADNE_TEST_FAIL(VP(Array<SP>{}));
    ARIADNE_TEST_FAIL(VP(List<SP>{}));
    VP null_vector,copied_vector_null(null_vector),vector_copy(vector);
    ARIADNE_TEST_ASSERT(null_vector.raw_pointer()==nullptr);
    ARIADNE_TEST_ASSERT(copied_vector_null.raw_pointer()==nullptr);
    vector_copy=null_vector;ARIADNE_TEST_ASSERT(vector_copy.raw_pointer()==nullptr);
    vector_copy=vector;ARIADNE_TEST_EQUALS(vector_copy(point),vector(point));
    ARIADNE_TEST_ASSERT(vector.managed_pointer()!=nullptr);
    ARIADNE_TEST_EQUALS(vector.result_size(),2u);
    ARIADNE_TEST_EQUALS(vector.argument_size(),2u);
    ARIADNE_TEST_EQUALS(vector.codomain().dimension(),2u);
    ARIADNE_TEST_EQUALS(vector.range().dimension(),2u);
    ARIADNE_TEST_EQUALS(vector.errors()[0u].get(dp).raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(vector.error().get(dp).raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(listed[Range(0u,2u)](point),vector(point));
    ARIADNE_TEST_EQUALS(listed[Range(1u,1u)](point).size(),0u);
    listed[0u]=listed[1u];ARIADNE_TEST_EQUALS(listed(point),Vector<FloatDPBounds>({2,2},dp));
    ARIADNE_TEST_EQUALS(listed[0u].argument_size(),2u);
    listed.clobber();ARIADNE_TEST_EQUALS(listed.error().get(dp).raw(),FloatDP(0,dp));
    const auto jacobian=vector.jacobian(generic_point);
    ARIADNE_TEST_ASSERT(definitely(jacobian[0u][0u]==ValidatedNumber(2)));
    ARIADNE_TEST_ASSERT(definitely(jacobian[0u][1u]==ValidatedNumber(1)));
    ARIADNE_TEST_ASSERT(definitely(jacobian[1u][1u]==ValidatedNumber(-1)));
    FunctionPatchFactory<ValidatedTag> patch_factory(make_taylor_function_patch_factory());
    const FunctionPatchFactory<ValidatedTag> shared(patch_factory.managed_pointer());
    ARIADNE_TEST_ASSERT(&static_cast<const FunctionPatchFactoryInterface<ValidatedTag>&>(shared)==patch_factory.managed_pointer().get());
    ARIADNE_TEST_EQUALS(shared.create_constants(domain,Vector<ValidatedNumber>{2,3})(point),Vector<FloatDPBounds>({2,3},dp));
    ARIADNE_TEST_EQUALS(shared.create_identity(domain)(point),point);
    ARIADNE_TEST_EQUALS(shared.create_projection(domain,Range(1u,2u))(point),Vector<FloatDPBounds>{point[1u]});
    auto creator=factory(original);
    ARIADNE_TEST_SAME(creator.create_zero()(point),FloatDPBounds(0,dp));
    ARIADNE_TEST_EQUALS(creator.create_zeros(0u)(point).size(),0u);
    ARIADNE_TEST_SAME(creator.create_constant(ValidatedNumber(4))(point),FloatDPBounds(4,dp));
    ARIADNE_TEST_EQUALS(creator.create_identity()(point),point);
    const auto no_results=shared.create_zeros(0u,domain);
    ARIADNE_TEST_EQUALS(no_results.jacobian(generic_point).row_size(),0u);
    const SP no_arguments=shared.create_constant(BoxDomainType(0u),ValidatedNumber(5));
    const VP constant_vector(1u,no_arguments);
    const auto zero_columns=constant_vector.jacobian(Vector<ValidatedNumber>(0u,ValidatedNumber(0)));
    ARIADNE_TEST_EQUALS(zero_columns.row_size(),1u);
    ARIADNE_TEST_EQUALS(zero_columns.column_size(),0u);
    std::ostringstream text;text<<original.reference()<<vector.reference();
    ARIADNE_TEST_ASSERT(!text.str().empty());
}

Int main() {
    ARIADNE_TEST_CALL(test_patch_handle_construction_and_values());
    ARIADNE_TEST_CALL(test_patch_constant_dispatch());
    ARIADNE_TEST_CALL(test_patch_vector_product());
    ARIADNE_TEST_CALL(test_mixed_patch_operations());
    return ARIADNE_TEST_FAILURES;
}
