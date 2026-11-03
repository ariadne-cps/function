/***************************************************************************
 * Tests of concrete scaled patches, domain contracts and numeric evaluation.
 * Distributed under the GNU GPL version 3 or later.
 ****************************************************************************/
#include "function/taylor_function.hpp"
#include "function/scaled_function_patch.tpl.hpp"
#include "ariadne-testing.hpp"
#include <sstream>
using namespace Ariadne;
using M=ValidatedTaylorModelDP;
using S=ScaledFunctionPatch<M>;
using V=VectorScaledFunctionPatch<M>;

class TestScaledFunctionPatch : public Ariadne::Testing::ContextFixture {
public:
    using ContextFixture::ContextFixture;
Void test_scaled_patch_values_and_domains() {
    const Sweeper<FloatDP> swp=ThresholdSweeper<FloatDP>(dp,1e-12);
    const BoxDomainType domain{{-1,1},{-1,1}},other{{0,2},{0,2}};
    const auto x=S::coordinate(domain,0u,swp),y=S::coordinate(domain,1u,swp);
    const V vector(domain,Vector<M>{x.model(),y.model()});
    const Vector<FloatDPBounds> point({0.5_x,-0.25_x},dp),outside({2,2},dp);
    const Vector<FloatDPApproximation> approximate({0.5_x,-0.25_x},dp),approximate_outside({2,2},dp);
    ARIADNE_TESTING_EXPECT_SAME(x(approximate), FloatDPApproximation(0.5_x,dp));
    ARIADNE_TESTING_EXPECT_THROWS(DomainException, x(approximate_outside));
    ARIADNE_TESTING_EXPECT_THROWS(DomainException, x(outside));
    ARIADNE_TESTING_EXPECT_THROWS(DomainException, vector(approximate_outside));
    ARIADNE_TESTING_EXPECT_THROWS(DomainException, vector(outside));
    ARIADNE_TESTING_EXPECT_EQUAL(vector(Vector<ValidatedNumber>{0.5_x,-0.25_x}), Vector<ValidatedNumber>({0.5_x,-0.25_x}));
    ARIADNE_TESTING_EXPECT_EQUAL(evaluate(vector,Vector<ValidatedNumber>{0.5_x,-0.25_x}), Vector<ValidatedNumber>({0.5_x,-0.25_x}));
    ARIADNE_TESTING_EXPECT_EQUAL(unchecked_evaluate(vector,Vector<ValidatedNumber>{0.5_x,-0.25_x}), Vector<ValidatedNumber>({0.5_x,-0.25_x}));
    ARIADNE_TESTING_EXPECT_EQUAL(x.characteristics(), S::CharacteristicsType(domain,swp));
    auto simplified=x;simplified.simplify(swp);simplified.simplify();
    ARIADNE_TESTING_EXPECT_EQUAL(simplified.expansion().argument_size(), 2u);
    ARIADNE_TESTING_EXPECT_SAME(simplified(point), point[0u]);
    ScaledFunctionPatchFactory<M> factory(swp);
    ARIADNE_TESTING_EXPECT_SAME(factory.create_constant(domain,FloatDPBounds(2,dp))(point), FloatDPBounds(2,dp));
    ARIADNE_TESTING_EXPECT_SAME(factory.create_identity(IntervalDomainType(-1,1))(Vector<FloatDPBounds>{point[0u]}), point[0u]);
    ARIADNE_TESTING_EXPECT_SAME(factory.create_number(ValidatedNumber(2)), FloatDPBounds(2,dp));
    ARIADNE_TESTING_EXPECT_THROWS(DomainException, restriction(x,other));
    ARIADNE_TESTING_EXPECT_EXCEPTION(antiderivative(x,2u,FloatDPBounds(0,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(antiderivative(x,0u,FloatDPBounds(2,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(partial_evaluate(x,0u,FloatDPBounds(2,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(partial_evaluate(x,0u,ValidatedNumber(2)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(partial_evaluate(vector,0u,FloatDPBounds(2,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(derivative(vector,2u));
    ARIADNE_TESTING_EXPECT_EXCEPTION(antiderivative(vector,2u));
    ARIADNE_TESTING_EXPECT_EXCEPTION(antiderivative(vector,2u,FloatDPBounds(0,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(restriction(vector,other));
    ARIADNE_TESTING_PROVE(same(restriction(vector,domain),vector));
    const auto shifted=S::coordinate(other,0u,swp);
    const V shifted_vector(other,Vector<M>{shifted.model(),shifted.model()});
    ARIADNE_TESTING_EXPECT_EXCEPTION(join(vector,shifted_vector));
    ARIADNE_TESTING_EXPECT_EXCEPTION(refinement(x,shifted));
    using O=AlgebraOperations<S>;
    ARIADNE_TESTING_EXPECT_EXCEPTION(O::apply(BinaryElementaryOperator(Add()),x,shifted));
    ARIADNE_TESTING_PROVE(!same(vector,shifted_vector));
    const V changed(domain,Vector<M>{x.model()+1,y.model()});
    ARIADNE_TESTING_PROVE(!same(vector,changed));
    const V single(domain,Vector<M>{x.model()});
    ARIADNE_TESTING_PROVE(!same(vector,single));
    ARIADNE_TESTING_PROVE(!same(single,vector));
    auto difference=vector;difference-=vector;
    ARIADNE_TESTING_EXPECT_EQUAL(difference(point), Vector<FloatDPBounds>({0,0},dp));
    ARIADNE_TESTING_EXPECT_EXCEPTION(difference-=single);
    ARIADNE_TESTING_EXPECT_EXCEPTION(difference-=shifted_vector);
    ARIADNE_TESTING_EXPECT_EXCEPTION(refines(vector,single));
    ARIADNE_TESTING_EXPECT_EXCEPTION(inconsistent(vector,single));
    ARIADNE_TESTING_EXPECT_EXCEPTION(refinement(vector,single));
    // Composition through the public model interface must preserve values.
    ValidatedVectorMultivariateFunctionModelDP abstract_vector=vector;
    ARIADNE_TESTING_EXPECT_EQUAL(compose(abstract_vector,vector)(point), point);
    // Compile-time-dispatched model evaluation, distinct from numeric evaluation.
    const Vector<ApproximateNumber> unsupported{ApproximateNumber(0),ApproximateNumber(0)};
    ARIADNE_TESTING_EXPECT_EXCEPTION(x(unsupported));
    ARIADNE_TESTING_EXPECT_EXCEPTION(vector(unsupported));
    const Vector<M> model_arguments=M::coordinates(2u,swp);
    const auto model_results=vector.template operator()<M>(model_arguments);
    ARIADNE_TESTING_EXPECT_SAME(evaluate(model_results[1u],point), point[1u]);
}

Void test_scaled_patch_construction_contracts() {
    const Sweeper<FloatDP> swp=ThresholdSweeper<FloatDP>(dp,1e-12);
    const BoxDomainType domain{{-1,1},{-1,1}},other{{0,2},{0,2}};
    const auto x=S::coordinate(domain,0u,swp);
    ARIADNE_TESTING_EXPECT_EXCEPTION(S::coordinate(domain,2u,swp));
    ARIADNE_TESTING_EXPECT_EXCEPTION(V(List<S>{}));
    const auto scalar=ValidatedScalarMultivariateFunction::coordinate(1u,0u);
    ARIADNE_TESTING_EXPECT_EXCEPTION(S(domain,scalar,swp));
    const auto model=M::coordinate(1u,0u,swp);
    ARIADNE_TESTING_EXPECT_EXCEPTION(V(domain,Vector<M>{model}));
    // Avoid UniformArray's iterator-range constructor (-Wshadow with GCC 13).
    Vector<Expansion<MultiIndex,FloatDP>> expansions(1u,x.expansion());
    Vector<M::ErrorType> errors{FloatDPError(0u,dp)};
    ARIADNE_TESTING_EXPECT_EXCEPTION(V(domain,expansions,Vector<M::ErrorType>{},swp));
    ARIADNE_TESTING_EXPECT_EXCEPTION(V(BoxDomainType{{-1,1}},expansions,errors,swp));
    const V raw_vector(domain,expansions,Vector<FloatDP>{FloatDP(0,dp)},swp);
    // Keep the evaluated vector alive while Testing compares its element.
    const auto evaluated_raw_vector=raw_vector(Vector<FloatDPBounds>({0.5_x,0},dp));
    ARIADNE_TESTING_EXPECT_SAME(evaluated_raw_vector[0u], FloatDPBounds(0.5_x,dp));
    Vector<FloatDP> raw_errors(0u,FloatDP(0,dp));
    // The API must reject unequal numbers of expansions and error bounds.
    ARIADNE_TESTING_EXPECT_EXCEPTION(V(domain,expansions,raw_errors,swp));
    const auto vector=ValidatedVectorMultivariateFunction::identity(1u);
    ARIADNE_TESTING_EXPECT_EXCEPTION(V(domain,vector,swp));
    const auto shifted=S::coordinate(other,0u,swp);
    ARIADNE_TESTING_EXPECT_EXCEPTION(V(List<S>{x,shifted}));
    ARIADNE_TESTING_EXPECT_EXCEPTION(V(Vector<S>{x,shifted}));
    V mutable_vector(domain,Vector<M>{x.model()});
    ARIADNE_TESTING_EXPECT_EXCEPTION(mutable_vector.set(1u,x));
    ARIADNE_TESTING_EXPECT_EXCEPTION(mutable_vector.set(0u,shifted));
    ARIADNE_TESTING_EXPECT_EXCEPTION(mutable_vector.adjoin(shifted));
}

Void test_scaled_patch_numeric_helpers() {
    const MultiplePrecision mp(128_bits);
    const FloatDPError error(0.25_x,dp);
    ARIADNE_TESTING_EXPECT_SAME(convert_error_to_bounds(error), FloatDPBounds(-0.25_x,0.25_x,dp));
    ARIADNE_TESTING_EXPECT_SAME(convert_error_to_bounds(PositiveFloatDPUpperBound(0.25_x,dp)), FloatDPBounds(-0.25_x,0.25_x,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(convert_error_to_bounds(PositiveFloatDPApproximation(FloatDPApproximation(0.25_x,dp))).raw(), FloatDP(0,dp));
    ARIADNE_TESTING_EXPECT_SAME(convert_error_to_bounds(FloatMPError(0.25_x,mp)), FloatMPBounds(-0.25_x,0.25_x,mp));
    ARIADNE_TESTING_EXPECT_SAME(convert_error_to_bounds(PositiveFloatMPUpperBound(0.25_x,mp)), FloatMPBounds(-0.25_x,0.25_x,mp));
    ARIADNE_TESTING_EXPECT_EQUAL(convert_error_to_bounds(PositiveFloatMPApproximation(FloatMPApproximation(0.25_x,mp))).raw(), FloatMP(0,mp));
    ARIADNE_TESTING_PROVE(contains(IntervalDomainType(-1,1),FloatMPBounds(0,mp)));
    ARIADNE_TESTING_PROVE(contains(IntervalDomainType(-1,1),FloatMPApproximation(0,mp)));
    ARIADNE_TESTING_PROVE(contains(BoxDomainType{{-1,1}},Vector<FloatMPApproximation>{FloatMPApproximation(0,mp)}));
    std::ostringstream text;text<<representation(FloatMP(2,mp));
    ARIADNE_TESTING_PROVE(!text.str().empty());
    ARIADNE_TESTING_PROVE(FloatDPApproximation(2,dp)==Int(2));
    ARIADNE_TESTING_PROVE(FloatDPBounds(2,dp)!=Int(3));
    ARIADNE_TESTING_PROVE(FloatDPBounds(3,dp)>Int(2));
}

};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_scaled_function_patch");
    ARIADNE_TESTING_ADD_TEST(TestScaledFunctionPatch, test_scaled_patch_values_and_domains);
    ARIADNE_TESTING_ADD_TEST(TestScaledFunctionPatch, test_scaled_patch_construction_contracts);
    ARIADNE_TESTING_ADD_TEST(TestScaledFunctionPatch, test_scaled_patch_numeric_helpers);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
