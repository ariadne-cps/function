/***************************************************************************
 *            test_function_patch_operations.cpp
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

#include <array>
#include <sstream>

#include "function/function.hpp"
#include "function/function_patch.hpp"
#include "function/taylor_function.hpp"

#include "ariadne-testing.hpp"
using namespace Ariadne;

class TestFunctionPatchOperations : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    Void test_patch_constant_dispatch() {
        BoxDomainType domain{
            {-1, 1},
            {-1, 1}
        };
        TaylorFunctionFactory factory(ThresholdSweeper<FloatDP>(dp, 1e-12));
        auto model = factory.create_constant(domain, ValidatedNumber(2));
        ValidatedScalarMultivariateFunction f = model;
        Vector<FloatDPBounds> point({0, 0}, dp);
        ARIADNE_TESTING_PROVE(definitely((f + ValidatedNumber(3))(point) == FloatDPBounds(5, dp)));
        using A = AlgebraOperations<ValidatedScalarMultivariateFunction, ValidatedNumber>;
        ARIADNE_TESTING_PROVE(definitely(A::apply(BinaryElementaryOperator(Sub()), ValidatedNumber(3), f)(point) == FloatDPBounds(1, dp)));
        ARIADNE_TESTING_PROVE(definitely((f - ValidatedNumber(3))(point) == FloatDPBounds(-1, dp)));
        ARIADNE_TESTING_PROVE(definitely((f * ValidatedNumber(3))(point) == FloatDPBounds(6, dp)));
        ARIADNE_TESTING_PROVE(definitely((f / ValidatedNumber(2))(point) == FloatDPBounds(1, dp)));
    }
    Void test_patch_vector_product() {
        BoxDomainType domain{
            {-1, 1},
            {-1, 1}
        };
        TaylorFunctionFactory factory(ThresholdSweeper<FloatDP>(dp, 1e-12));
        ValidatedScalarMultivariateFunction scalar = factory.create_constant(domain, ValidatedNumber(2));
        ValidatedVectorMultivariateFunction vector = factory.create_constants(domain, Vector<ValidatedNumber>{4, 8});
        Vector<FloatDPBounds> point({0, 0}, dp);
        ARIADNE_TESTING_EXPECT_EQUAL((scalar * vector)(point), Vector<FloatDPBounds>({8, 16}, dp));
    }
    Void test_mixed_patch_operations() {
        using S = ValidatedScalarMultivariateFunction;
        using V = ValidatedVectorMultivariateFunction;
        BoxDomainType domain{
            {-1, 1},
            {-1, 1}
        };
        TaylorFunctionFactory factory(ThresholdSweeper<FloatDP>(dp, 1e-12));
        std::array<S, 2> scalar_a{S::constant(2u, ValidatedNumber(2)), factory.create_constant(domain, ValidatedNumber(2))};
        std::array<S, 2> scalar_b{S::constant(2u, ValidatedNumber(4)), factory.create_constant(domain, ValidatedNumber(4))};
        Vector<ValidatedNumber> va{2, 6}, vb{4, 8};
        std::array<V, 2> vector_a{V::constant(2u, va), factory.create_constants(domain, va)};
        std::array<V, 2> vector_b{V::constant(2u, vb), factory.create_constants(domain, vb)};
        Vector<FloatDPBounds> point({0, 0}, dp);
        for (SizeType i = 0; i != 2u; ++i) {
            const auto& a = scalar_a[i];
            ARIADNE_TESTING_PROVE(definitely((-a)(point) == FloatDPBounds(-2, dp)));
            ARIADNE_TESTING_PROVE(definitely(pow(a, 2)(point) == FloatDPBounds(4, dp)));
            ARIADNE_TESTING_PROVE(definitely((ValidatedNumber(3) + a)(point) == FloatDPBounds(5, dp)));
            ARIADNE_TESTING_PROVE(definitely((ValidatedNumber(3) - a)(point) == FloatDPBounds(1, dp)));
            for (SizeType j = 0; j != 2u; ++j) {
                const auto& b = scalar_b[j];
                const auto& av = vector_a[i];
                const auto& bv = vector_b[j];
                ARIADNE_TESTING_PROVE(definitely((a + b)(point) == FloatDPBounds(6, dp)));
                ARIADNE_TESTING_PROVE(definitely((a - b)(point) == FloatDPBounds(-2, dp)));
                ARIADNE_TESTING_PROVE(definitely((a * b)(point) == FloatDPBounds(8, dp)));
                ARIADNE_TESTING_PROVE(definitely((a / b)(point) == FloatDPBounds(0.5_x, dp)));
                ARIADNE_TESTING_EXPECT_EQUAL((av + bv)(point), Vector<FloatDPBounds>({6, 14}, dp));
                ARIADNE_TESTING_EXPECT_EQUAL((av - bv)(point), Vector<FloatDPBounds>({-2, -2}, dp));
                ARIADNE_TESTING_EXPECT_EQUAL((a * bv)(point), Vector<FloatDPBounds>({8, 16}, dp));
                ARIADNE_TESTING_EXPECT_EQUAL((bv * a)(point), Vector<FloatDPBounds>({8, 16}, dp));
                ARIADNE_TESTING_EXPECT_EQUAL((bv / a)(point), Vector<FloatDPBounds>({2, 4}, dp));
                ARIADNE_TESTING_EXPECT_EQUAL(join(a, b)(point), Vector<FloatDPBounds>({2, 4}, dp));
                ARIADNE_TESTING_EXPECT_EQUAL(join(a, bv)(point), Vector<FloatDPBounds>({2, 4, 8}, dp));
                ARIADNE_TESTING_EXPECT_EQUAL(join(av, b)(point), Vector<FloatDPBounds>({2, 6, 4}, dp));
                ARIADNE_TESTING_EXPECT_EQUAL(join(av, bv)(point), Vector<FloatDPBounds>({2, 6, 4, 8}, dp));
            }
        }
        ValidatedScalarMultivariateFunctionPatch patch = factory.create_constant(domain, ValidatedNumber(2));
        ValidatedScalarMultivariateFunctionModelDP model = factory.create_constant(domain, ValidatedNumber(2));
        S scalar_model = model;
        V one_component{scalar_model};
        ARIADNE_TESTING_EXPECT_EQUAL(one_component(point), Vector<FloatDPBounds>{FloatDPBounds(2, dp)});
        std::ostringstream text;
        text << patch << model << factory;
        ARIADNE_TESTING_PROVE(!text.str().empty());
    }
    Void test_patch_handle_construction_and_values() {
        using SP = ValidatedScalarMultivariateFunctionPatch;
        using VP = ValidatedVectorMultivariateFunctionPatch;
        const BoxDomainType domain{
            { 2, 6},
            {-2, 2}
        };
        const TaylorFunctionFactory backend(ThresholdSweeper<FloatDP>(dp, 1e-12));
        const auto x = backend.create_coordinate(domain, 0u), y = backend.create_coordinate(domain, 1u);
        const SP original = 2 * x + y, second = x - y;
        const Vector<FloatDPBounds> point({3, 1}, dp);
        const Vector<ValidatedNumber> generic_point{3, 1};
        SP null_scalar, copied_null(null_scalar), assigned(original);
        ARIADNE_TESTING_PROVE(copied_null.raw_pointer() == nullptr);
        assigned = null_scalar;
        ARIADNE_TESTING_PROVE(assigned.raw_pointer() == nullptr);
        assigned = original;
        ARIADNE_TESTING_EXPECT_SAME(assigned(point), FloatDPBounds(7, dp));
        SP copy(original);
        ARIADNE_TESTING_EXPECT_SAME(copy(point), FloatDPBounds(7, dp));
        SP from_function(static_cast<ValidatedScalarMultivariateFunction>(original));
        ARIADNE_TESTING_EXPECT_SAME(from_function(point), FloatDPBounds(7, dp));
        ARIADNE_TESTING_EXPECT_EXCEPTION(SP(ValidatedScalarMultivariateFunction::coordinate(2u, 0u)));
        ARIADNE_TESTING_EXPECT_EQUAL(original.result_size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(original.argument_size(), 2u);
        ARIADNE_TESTING_PROVE(same(original.codomain(), IntervalDomainType(2, 14)));
        ARIADNE_TESTING_EXPECT_EQUAL(original.error().get(dp).raw(), FloatDP(0, dp));
        ARIADNE_TESTING_PROVE(original.managed_pointer() != nullptr);
        ARIADNE_TESTING_EXPECT_SAME(derivative(original, 0u)(point), FloatDPBounds(2, dp));
        assigned = ValidatedNumber(5);
        ARIADNE_TESTING_EXPECT_SAME(assigned(point), FloatDPBounds(5, dp));
        assigned = ValidatedScalarMultivariateFunction::coordinate(2u, 0u);
        ARIADNE_TESTING_EXPECT_SAME(assigned(point), FloatDPBounds(3, dp));
        assigned = original.reference();
        ARIADNE_TESTING_EXPECT_SAME(assigned(point), FloatDPBounds(7, dp));
        SP::Interface& mutable_interface = assigned;
        ARIADNE_TESTING_EXPECT_EQUAL(mutable_interface.argument_size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(assigned.reference().argument_size(), 2u);
        using Operations = AlgebraOperations<SP, ValidatedNumber>;
        ARIADNE_TESTING_EXPECT_SAME(Operations::apply(BinaryElementaryOperator(Sub()), ValidatedNumber(5), original)(point), FloatDPBounds(-2, dp));
        SP uncertain = backend.create_constant(domain, ValidatedNumber(FloatDPBounds(1, 3, dp)));
        uncertain.clobber();
        ARIADNE_TESTING_EXPECT_EQUAL(uncertain.error().get(dp).raw(), FloatDP(0, dp));
        ARIADNE_TESTING_EXPECT_SAME(uncertain(point), FloatDPBounds(2, dp));
        const VP vector(Array<SP>{original, second});
        VP listed(List<SP>{original, second});
        ARIADNE_TESTING_EXPECT_EQUAL(vector(point), Vector<FloatDPBounds>({7, 2}, dp));
        ARIADNE_TESTING_EXPECT_EQUAL(listed(point), vector(point));
        ARIADNE_TESTING_EXPECT_EQUAL(listed.reference().result_size(), 2u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(VP(Array<SP>{}));
        ARIADNE_TESTING_EXPECT_EXCEPTION(VP(List<SP>{}));
        VP null_vector, copied_vector_null(null_vector), vector_copy(vector);
        ARIADNE_TESTING_PROVE(null_vector.raw_pointer() == nullptr);
        ARIADNE_TESTING_PROVE(copied_vector_null.raw_pointer() == nullptr);
        vector_copy = null_vector;
        ARIADNE_TESTING_PROVE(vector_copy.raw_pointer() == nullptr);
        vector_copy = vector;
        ARIADNE_TESTING_EXPECT_EQUAL(vector_copy(point), vector(point));
        ARIADNE_TESTING_PROVE(vector.managed_pointer() != nullptr);
        ARIADNE_TESTING_EXPECT_EQUAL(vector.result_size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(vector.argument_size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(vector.codomain().dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(vector.range().dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(vector.errors()[0u].get(dp).raw(), FloatDP(0, dp));
        ARIADNE_TESTING_EXPECT_EQUAL(vector.error().get(dp).raw(), FloatDP(0, dp));
        ARIADNE_TESTING_EXPECT_EQUAL(listed[Range(0u, 2u)](point), vector(point));
        ARIADNE_TESTING_EXPECT_EQUAL(listed[Range(1u, 1u)](point).size(), 0u);
        listed[0u] = listed[1u];
        ARIADNE_TESTING_EXPECT_EQUAL(listed(point), Vector<FloatDPBounds>({2, 2}, dp));
        ARIADNE_TESTING_EXPECT_EQUAL(listed[0u].argument_size(), 2u);
        listed.clobber();
        ARIADNE_TESTING_EXPECT_EQUAL(listed.error().get(dp).raw(), FloatDP(0, dp));
        const auto jacobian = vector.jacobian(generic_point);
        ARIADNE_TESTING_PROVE(definitely(jacobian[0u][0u] == ValidatedNumber(2)));
        ARIADNE_TESTING_PROVE(definitely(jacobian[0u][1u] == ValidatedNumber(1)));
        ARIADNE_TESTING_PROVE(definitely(jacobian[1u][1u] == ValidatedNumber(-1)));
        FunctionPatchFactory<ValidatedTag> patch_factory(make_taylor_function_patch_factory());
        const FunctionPatchFactory<ValidatedTag> shared(patch_factory.managed_pointer());
        ARIADNE_TESTING_PROVE(&static_cast<const FunctionPatchFactoryInterface<ValidatedTag>&>(shared) == patch_factory.managed_pointer().get());
        ARIADNE_TESTING_EXPECT_EQUAL(shared.create_constants(domain, Vector<ValidatedNumber>{2, 3})(point), Vector<FloatDPBounds>({2, 3}, dp));
        ARIADNE_TESTING_EXPECT_EQUAL(shared.create_identity(domain)(point), point);
        ARIADNE_TESTING_EXPECT_EQUAL(shared.create_projection(domain, Range(1u, 2u))(point), Vector<FloatDPBounds>{point[1u]});
        auto creator = factory(original);
        ARIADNE_TESTING_EXPECT_SAME(creator.create_zero()(point), FloatDPBounds(0, dp));
        ARIADNE_TESTING_EXPECT_EQUAL(creator.create_zeros(0u)(point).size(), 0u);
        ARIADNE_TESTING_EXPECT_SAME(creator.create_constant(ValidatedNumber(4))(point), FloatDPBounds(4, dp));
        ARIADNE_TESTING_EXPECT_EQUAL(creator.create_identity()(point), point);
        const auto no_results = shared.create_zeros(0u, domain);
        ARIADNE_TESTING_EXPECT_EQUAL(no_results.jacobian(generic_point).row_size(), 0u);
        const SP no_arguments = shared.create_constant(BoxDomainType(0u), ValidatedNumber(5));
        const VP constant_vector(1u, no_arguments);
        const auto zero_columns = constant_vector.jacobian(Vector<ValidatedNumber>(0u, ValidatedNumber(0)));
        ARIADNE_TESTING_EXPECT_EQUAL(zero_columns.row_size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(zero_columns.column_size(), 0u);
        std::ostringstream text;
        text << original.reference() << vector.reference();
        ARIADNE_TESTING_PROVE(!text.str().empty());
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_function_patch_operations");
    ARIADNE_TESTING_ADD_TEST(TestFunctionPatchOperations, test_patch_handle_construction_and_values);
    ARIADNE_TESTING_ADD_TEST(TestFunctionPatchOperations, test_patch_constant_dispatch);
    ARIADNE_TESTING_ADD_TEST(TestFunctionPatchOperations, test_patch_vector_product);
    ARIADNE_TESTING_ADD_TEST(TestFunctionPatchOperations, test_mixed_patch_operations);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
