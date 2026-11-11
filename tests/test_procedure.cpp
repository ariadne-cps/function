/***************************************************************************
 *            test_procedure.cpp
 *
 *  Copyright  2010-20  Pieter Collins
 *
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Function.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <type_traits>
#include "function/procedure.hpp"
#include "numeric/numeric.hpp"
#include "interval/interval.hpp"

namespace Ariadne {

// Test-side overloads used to instantiate the generic backward contractors.
// The UpperIntervalType overload mirrors the production implementation in
// procedure.cpp; the FloatDPBounds overload is used only by the isolated Leq
// test, since Leq is not a ProcedureInstruction operator.
inline Void restrict(UpperIntervalType& r, UpperIntervalType const& x) {
    r.set_lower_bound(max(r.lower_bound(),x.lower_bound()));
    r.set_upper_bound(min(r.upper_bound(),x.upper_bound()));
}
inline Void restrict(FloatDPBounds& r, FloatDPBounds const& x) {
    r=refinement(r,x);
}

} // namespace Ariadne

#include "function/procedure.tpl.hpp"
#include "algebra/graded.hpp"
#include "algebra/vector.hpp"
#include "algebra/expansion.hpp"
#include "algebra/evaluate.hpp"
#include "function/formula.hpp"

#include "ariadne-testing.hpp"

#include "algebra/covector.hpp"
#include "function/function.hpp"
#include "algebra/expansion.inl.hpp"

using namespace std;
using namespace Ariadne;

namespace Ariadne {
template<class X> decltype(auto) mag(Covector<X> const& u) { return norm(transpose(u)); }
}

class TestProcedure : public Ariadne::Testing::ContextFixture
{
    DoublePrecision pr;
  public:
    explicit TestProcedure(Ariadne::Testing::Context& context);
  public:
    Void test_construct_from_formula();
    Void test_construct_from_expansion();
    Void test_evaluate();
    Void test_propagate();
    Void test_backward_contractor_soundness();
    Void test_backward_contractor_witness_preservation();
    Void test_backward_contractor_generated_witnesses();
    Void test_leq_backpropagate_soundness();
    Void test_inverse_trigonometric_empty_backpropagation();
    Void test_abs_backpropagation();
    Void test_nul_hlf_backpropagation();
    Void test_elementary_empty_backpropagation();
    Void test_derivative();
    Void test_scalar_formula_conversion();
    Void test_unsupported_graded_operations();
};

TestProcedure::TestProcedure(Ariadne::Testing::Context& context)
    : ContextFixture(context)
{
}




Void TestProcedure::test_unsupported_graded_operations()
{
    // Test the typed dispatch branches: visiting a variant operator would
    // instantiate all arithmetic alternatives, including combinations not
    // supported by this standalone Graded coefficient type.
    Graded<RoundedFloatDP> value(RoundedFloatDP(2,dp));
    RoundedFloatDP scalar(1,dp);
    ARIADNE_TESTING_EXPECT_EXCEPTION(ProcedureDetail::propagate(value,Abs(),value));
    ARIADNE_TESTING_EXPECT_EXCEPTION(ProcedureDetail::propagate(value,Max(),value,value));
    ARIADNE_TESTING_EXPECT_EXCEPTION(ProcedureDetail::propagate(value,Min(),value,value));
    ARIADNE_TESTING_EXPECT_EXCEPTION(ProcedureDetail::propagate(value,Max(),scalar,value));
    ARIADNE_TESTING_EXPECT_EXCEPTION(ProcedureDetail::propagate(value,Min(),scalar,value));
    ARIADNE_TESTING_EXPECT_EXCEPTION(ProcedureDetail::propagate(value,Sub(),scalar,value));
    ARIADNE_TESTING_EXPECT_EXCEPTION(ProcedureDetail::propagate(value,Div(),scalar,value));
}

Void TestProcedure::test_construct_from_formula()
{
    ApproximateNumber c(2.0_x);
    ApproximateFormula o(ApproximateFormula::constant(1.0_x));
    ApproximateFormula x(ApproximateFormula::coordinate(0));
    ApproximateFormula y(ApproximateFormula::coordinate(1));

    ApproximateFormula xs=pow(x,2);
    ApproximateFormula ys=pow(y,2);

    ARIADNE_TESTING_INSPECT("y*c", y*c);
    ARIADNE_TESTING_INSPECT("c*y", c*y);
    Vector<ApproximateFormula> f={sqrt(xs+c*ys), atan(y/x), xs*c-ys};
    ARIADNE_TESTING_INSPECT("f", f);

    ApproximateProcedure p0(argument_size=2u,f[0]);
    ARIADNE_TESTING_EXPECT_EQUAL(p0.argument_size(), 2u);
    ARIADNE_TESTING_INSPECT("p0", p0);
    Vector<ApproximateProcedure> p(argument_size=2u,f);
    ARIADNE_TESTING_INSPECT("p", p);
    ARIADNE_TESTING_EXPECT_EQUAL(p.result_size(), 3u);
    ARIADNE_TESTING_EXPECT_EQUAL(p.argument_size(), 2u);
}


Void TestProcedure::test_construct_from_expansion()
{
    {
        Expansion<MultiIndex,FloatDPApproximation> e({ {{0,0},1.0_x}, {{1,0},2.0_x}, {{0,2},3.0_x}, {{1,4},4.0_x} },pr);
        ARIADNE_TESTING_INSPECT("e", e);
        e.reverse_lexicographic_sort();
        ARIADNE_TESTING_INSPECT("e", e);
        Procedure<ApproximateNumber> p(e);
        ARIADNE_TESTING_INSPECT("p", p);
        Vector<FloatDPApproximation> x({2.0_x,3.0_x},pr);
        ARIADNE_TESTING_INSPECT("simple_evaluate(e,x)", simple_evaluate(e,x));
        ARIADNE_TESTING_INSPECT("evaluate(p,x)", evaluate(p,x));
        ARIADNE_TESTING_EXPECT_SAME(evaluate(p,x), simple_evaluate(e,x));
    }

    {
        Expansion<MultiIndex,FloatDPApproximation> e({ {{0,0},1.0_x}, {{1,0},2.0_x}, {{0,1},3.0_x}, {{2,0},4.0_x}, {{1,1},5.0_x}, {{0,2},6.0_x} },pr);
        e.reverse_lexicographic_sort();
        Procedure<ApproximateNumber> p(e);
        ARIADNE_TESTING_INSPECT("p", p);
        Vector<FloatDPApproximation> x({2.0_x,3.0_x},pr);
        ARIADNE_TESTING_EXPECT_SAME(evaluate(p,x), simple_evaluate(e,x));
    }
}


Void TestProcedure::test_evaluate()
{
    ApproximateProcedure p(2);
    p.new_instruction(Var(),0ul);
    p.new_instruction(Var(),1ul);
    p.new_instruction(Pow(),0ul,2);
    p.new_instruction(Sqr(),1ul);
    p.new_instruction(Add(),2ul,3ul);
    p.new_constant(9.0_x);
    p.new_instruction_scalar(Mul(),0ul,4ul);
    p.new_instruction(Sqrt(),5ul);
    ARIADNE_TESTING_INSPECT("p", p);

    Vector<FloatDPApproximation> x({3.0_x,4.0_x},pr);
    ARIADNE_TESTING_INSPECT("x", x);

    ARIADNE_TESTING_EXPECT_SAME(evaluate(p,x), FloatDPApproximation(15,pr));

    ValidatedProcedure tanh_procedure(1);
    tanh_procedure.new_instruction(Var(),0u);
    tanh_procedure.new_instruction(Tanh(),0u);
    Vector<FloatDPBounds> tanh_argument({FloatDPBounds(-1,1,dp)});
    ARIADNE_TESTING_EXPECT_SAME(evaluate(tanh_procedure,tanh_argument), tanh(tanh_argument[0]));

    // Regression test: vector procedure evaluation must support validated
    // numeric types that require an explicit precision at construction.
    ValidatedFormula px(ValidatedFormula::coordinate(0));
    ValidatedFormula py(ValidatedFormula::coordinate(1));
    Vector<ValidatedFormula> formulas={px+py,px-py};
    Vector<ValidatedProcedure> procedures(argument_size=2u,formulas);
    Vector<FloatDPBounds> arguments({
        FloatDPBounds(2,dp),FloatDPBounds(3,dp)
    });
    Vector<FloatDPBounds> values=evaluate(procedures,arguments);
    ARIADNE_TESTING_EXPECT_EQUAL(values[0], FloatDPBounds(5,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(values[1], FloatDPBounds(-1,dp));
}

Void TestProcedure::test_propagate()
{
    {
        ValidatedProcedure p(2);
        p.new_instruction(Var(),0u);
        p.new_instruction(Var(),1u);
        p.new_instruction(Sqr(),0u);
        p.new_instruction(Sqr(),1u);
        p.new_instruction(Add(),2u,3u);
        p.new_constant(1.125_decimal);
        p.new_instruction_scalar(Mul(),0ul,4ul);
        p.new_instruction(Sqrt(),5u);
        ARIADNE_TESTING_INSPECT("p", p);

        UpperBoxType x=ExactBoxType{ {0.25_x,2.0_x}, {0.5_x,3.0_x} };
        ARIADNE_TESTING_INSPECT("x", x);

        ProcedureHullReductionStatistics statistics;
        simple_hull_reduce(x,p,ExactIntervalType(1,1),statistics);
        ARIADNE_TESTING_INSPECT("x", x);
        ARIADNE_TESTING_PROVE(statistics.temporary_allocation_seconds>=0.0);
        ARIADNE_TESTING_PROVE(statistics.forward_execution_seconds>=0.0);
        ARIADNE_TESTING_PROVE(statistics.backward_propagation_seconds>=0.0);
        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TESTING_INSPECT("x", x);
    }

    ValidatedFormula x(ValidatedFormula::coordinate(0));
    ValidatedFormula y(ValidatedFormula::coordinate(1));

    Vector<ValidatedFormula> ff={sqrt(sqr(x)+sqr(y)),2*x-y};
    ARIADNE_TESTING_INSPECT("ff", ff);
    Vector<ValidatedProcedure> pp(argument_size=2,ff);
    ARIADNE_TESTING_INSPECT("pp", pp);
    UpperBoxType xx=ExactBoxType{ {0.125_x,2.0_x}, {0.25_x,3.0_x} };
    ExactBoxType cc=ExactBoxType{ {1.0_x,1.0_x}, {1.0_x,1.0_x} };
    ARIADNE_TESTING_INSPECT("xx", xx);
    ARIADNE_TESTING_INSPECT("evaluate(pp,cast_vector(xx))", evaluate(pp,cast_vector(xx)));
    simple_hull_reduce(xx,pp,cc);
    ARIADNE_TESTING_INSPECT("xx", xx);
}

Void TestProcedure::test_backward_contractor_soundness()
{
    auto check_periodic_witness = [this](auto op, ExactIntervalType domain) {
        ARIADNE_TESTING_NOTE("[contractor-soundness] periodic op=" << op
                  << " domain=" << domain
                  << " codomain={0:0}");
        ValidatedProcedure p(1);
        p.new_instruction(Var(),0u);
        p.new_instruction(op,0u);
        UpperBoxType x=ExactBoxType({domain});
        simple_hull_reduce(x,p,ExactIntervalType(0,0));
        ARIADNE_TESTING_NOTE("[contractor-soundness] periodic op=" << op
                  << " contracted=" << x[0]);
        ARIADNE_TESTING_PROVE(!definitely(x[0].is_empty()));
    };

    // Each interval contains a non-principal zero of the corresponding
    // periodic function. Backward propagation must not discard that branch.
    check_periodic_witness(Sin(),ExactIntervalType(3.125_x,3.25_x));
    check_periodic_witness(Cos(),ExactIntervalType(4.5_x,5.0_x));
    check_periodic_witness(Tan(),ExactIntervalType(3.125_x,3.25_x));

    {
        ARIADNE_TESTING_NOTE("[contractor-soundness] unary op=sqr domain={-1:-1} codomain={1:1}");
        ValidatedProcedure p(1);
        p.new_instruction(Var(),0u);
        p.new_instruction(Sqr(),0u);
        UpperBoxType x=ExactBoxType({ExactIntervalType(-1,-1)});
        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TESTING_PROVE(!definitely(x[0].is_empty()));
    }

    {
        ARIADNE_TESTING_NOTE("[contractor-soundness] power exponent=2 domain={-1:-1} codomain={1:1}");
        ValidatedProcedure p(1);
        p.new_instruction(Var(),0u);
        p.new_instruction(Pow(),0u,2);
        UpperBoxType x=ExactBoxType({ExactIntervalType(-1,-1)});
        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TESTING_PROVE(!definitely(x[0].is_empty()));
    }

    {
        ARIADNE_TESTING_NOTE("[contractor-soundness] binary op=max lhs={0:0} rhs={1:1} codomain={1:1}");
        ValidatedProcedure p(2);
        p.new_instruction(Var(),0u);
        p.new_instruction(Var(),1u);
        p.new_instruction(Max(),0u,1u);
        UpperBoxType x=ExactBoxType({ExactIntervalType(0,0),ExactIntervalType(1,1)});
        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TESTING_PROVE(!definitely(x[0].is_empty()));
        ARIADNE_TESTING_PROVE(!definitely(x[1].is_empty()));
    }

    {
        ARIADNE_TESTING_NOTE("[contractor-soundness] binary op=min lhs={1:1} rhs={0:0} codomain={0:0}");
        ValidatedProcedure p(2);
        p.new_instruction(Var(),0u);
        p.new_instruction(Var(),1u);
        p.new_instruction(Min(),0u,1u);
        UpperBoxType x=ExactBoxType({ExactIntervalType(1,1),ExactIntervalType(0,0)});
        simple_hull_reduce(x,p,ExactIntervalType(0,0));
        ARIADNE_TESTING_PROVE(!definitely(x[0].is_empty()));
        ARIADNE_TESTING_PROVE(!definitely(x[1].is_empty()));
    }

    {
        ValidatedProcedure p(1);
        p.new_instruction(Var(),0u);
        p.new_constant(1.0_x);
        p.new_instruction_scalar(Max(),0u,0u);
        UpperBoxType x=ExactBoxType({ExactIntervalType(0,0)});
        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TESTING_PROVE(!definitely(x[0].is_empty()));
    }

    {
        ValidatedProcedure p(1);
        p.new_instruction(Var(),0u);
        p.new_constant(0.0_x);
        p.new_instruction_scalar(Min(),0u,0u);
        UpperBoxType x=ExactBoxType({ExactIntervalType(1,1)});
        simple_hull_reduce(x,p,ExactIntervalType(0,0));
        ARIADNE_TESTING_PROVE(!definitely(x[0].is_empty()));
    }

    // Conservative periodic backward propagation must still propagate an
    // empty forward image so that an impossible codomain empties the domain.
    {
        ValidatedProcedure p(1);
        p.new_instruction(Var(),0u);
        p.new_instruction(Sin(),0u);
        UpperBoxType x=ExactBoxType({ExactIntervalType(0,0)});
        simple_hull_reduce(x,p,ExactIntervalType(2,2));
        ARIADNE_TESTING_PROVE(x[0].is_empty());
    }
}


Void TestProcedure::test_backward_contractor_witness_preservation()
{
    auto U = [](ExactIntervalType const& x) { return UpperIntervalType(x); };
    auto E = [](auto l, auto u) { return ExactIntervalType(l,u); };
    auto S = [](auto x) { return ExactIntervalType(x,x); };

    auto preserves = [](UpperIntervalType const& contracted, ExactIntervalType const& witness) {
        return !definitely(disjoint(contracted,witness));
    };

    auto check_unary = [&](auto op, ExactIntervalType domain, ExactIntervalType output, ExactIntervalType witness) {
        UpperIntervalType contracted=U(domain);
        ARIADNE_TESTING_NOTE("[contractor-witness] unary op=" << op
                  << " domain=" << domain
                  << " output=" << output
                  << " witness=" << witness);
        backpropagate(U(output),op,contracted);
        ARIADNE_TESTING_NOTE("[contractor-witness] unary op=" << op
                  << " contracted=" << contracted);
        ARIADNE_TESTING_PROVE(preserves(contracted,witness));
    };

    auto check_binary = [&](auto op,
                            ExactIntervalType lhs_domain, ExactIntervalType rhs_domain,
                            ExactIntervalType output,
                            ExactIntervalType lhs_witness, ExactIntervalType rhs_witness) {
        UpperIntervalType lhs=U(lhs_domain);
        UpperIntervalType rhs=U(rhs_domain);
        ARIADNE_TESTING_NOTE("[contractor-witness] binary op=" << op
                  << " lhs_domain=" << lhs_domain
                  << " rhs_domain=" << rhs_domain
                  << " output=" << output
                  << " lhs_witness=" << lhs_witness
                  << " rhs_witness=" << rhs_witness);
        backpropagate(U(output),op,lhs,rhs);
        ARIADNE_TESTING_NOTE("[contractor-witness] binary op=" << op
                  << " lhs_contracted=" << lhs
                  << " rhs_contracted=" << rhs);
        ARIADNE_TESTING_PROVE(preserves(lhs,lhs_witness));
        ARIADNE_TESTING_PROVE(preserves(rhs,rhs_witness));
    };

    auto check_power = [&](ExactIntervalType domain, Int exponent,
                           ExactIntervalType output, ExactIntervalType witness) {
        UpperIntervalType contracted=U(domain);
        ARIADNE_TESTING_NOTE("[contractor-witness] power exponent=" << exponent
                  << " domain=" << domain
                  << " output=" << output
                  << " witness=" << witness);
        backpropagate(U(output),Pow(),contracted,exponent);
        ARIADNE_TESTING_NOTE("[contractor-witness] power exponent=" << exponent
                  << " contracted=" << contracted);
        ARIADNE_TESTING_PROVE(preserves(contracted,witness));
    };

    // Unary inverse contractors.
    check_unary(Pos(), S(2), S(2), S(2));
    check_unary(Neg(), S(2), S(-2), S(2));
    check_unary(Rec(), S(2), S(0.5_x), S(2));

    check_unary(Sqr(), S(-2), S(4), S(-2));
    check_unary(Sqr(), S( 2), S(4), S(2));
    check_unary(Sqr(), S( 0), S(0), S(0));

    check_unary(Sqrt(), S(4), S(2), S(4));
    check_unary(Sqrt(), S(0), S(0), S(0));
    check_unary(Exp(), S(0), S(1), S(0));
    check_unary(Log(), S(1), S(0), S(1));

    // Principal and non-principal periodic witnesses.
    check_unary(Sin(), S(0), S(0), S(0));
    check_unary(Cos(), S(0), S(1), S(0));
    check_unary(Tan(), S(0), S(0), S(0));
    check_unary(Tanh(), S(0), S(0), S(0));
    check_unary(Sin(), E(3.125_x,3.25_x), S(0), E(3.125_x,3.25_x));
    check_unary(Cos(), E(4.5_x,5.0_x), S(0), E(4.5_x,5.0_x));
    check_unary(Tan(), E(3.125_x,3.25_x), S(0), E(3.125_x,3.25_x));

    // Inverse trigonometric functions on their natural domains.
    check_unary(Asin(), S(0), S(0), S(0));
    check_unary(Acos(), S(1), S(0), S(1));
    check_unary(Atan(), S(0), S(0), S(0));

    // Binary contractors, including signs and zero products.
    check_binary(Add(), S(2), S(3), S(5), S(2), S(3));
    check_binary(Sub(), S(2), S(3), S(-1), S(2), S(3));
    check_binary(Mul(), S(2), S(3), S(6), S(2), S(3));
    check_binary(Mul(), S(-2), S(3), S(-6), S(-2), S(3));
    check_binary(Mul(), S(0), S(3), S(0), S(0), S(3));
    check_binary(Div(), S(6), S(3), S(2), S(6), S(3));
    check_binary(Div(), S(-6), S(3), S(-2), S(-6), S(3));

    check_binary(Max(), S(0), S(1), S(1), S(0), S(1));
    check_binary(Max(), S(1), S(0), S(1), S(1), S(0));
    check_binary(Max(), S(1), S(1), S(1), S(1), S(1));
    check_binary(Min(), S(0), S(1), S(0), S(0), S(1));
    check_binary(Min(), S(1), S(0), S(0), S(1), S(0));
    check_binary(Min(), S(0), S(0), S(0), S(0), S(0));

    // Integer powers: both branches for even powers, signed odd powers,
    // zero exponent, and negative exponents.
    check_power(S(-2),  2, S(4), S(-2));
    check_power(S( 2),  2, S(4), S(2));
    check_power(S(-2),  3, S(-8), S(-2));
    check_power(S( 2),  3, S(8), S(2));
    check_power(S( 2),  0, S(1), S(2));
    check_power(S( 2), -1, S(0.5_x), S(2));
    check_power(S(-2), -2, S(0.25_x), S(-2));

    // Non-degenerate domains exercise endpoint propagation.
    check_unary(Pos(), E(1,3), E(1,3), S(2));
    check_unary(Neg(), E(1,3), E(-3,-1), S(2));
    check_unary(Rec(), E(1,4), E(0.25_x,1.0_x), S(2));
    check_unary(Sqr(), E(-3,3), E(1,4), S(-2));
    check_unary(Sqr(), E(-3,3), E(1,4), S(2));
    check_unary(Sqrt(), E(0,9), E(1,2), S(4));
    check_unary(Exp(), E(-1,2), E(0.5_x,2.0_x), S(0));
    check_unary(Log(), E(0.5_x,3.0_x), E(-0.5_x,1.0_x), S(1));
    check_unary(Sin(), E(3,4), E(-0.5_x,0.5_x), E(3.125_x,3.25_x));
    check_unary(Cos(), E(4,5), E(-0.5_x,0.5_x), E(4.5_x,5.0_x));
    check_unary(Tan(), E(3,4), E(-0.5_x,0.5_x), E(3.125_x,3.25_x));
    check_unary(Asin(), E(-1,1), E(-0.5_x,0.5_x), S(0));
    check_unary(Acos(), E(-1,1), E(0,2), S(1));
    check_unary(Atan(), E(-2,2), E(-1,1), S(0));

    check_binary(Add(), E(0,4), E(1,5), E(4,6), S(2), S(3));
    check_binary(Sub(), E(0,4), E(1,5), E(-2,0), S(2), S(3));
    check_binary(Mul(), E(-3,-1), E(2,4), E(-8,-4), S(-2), S(3));
    check_binary(Div(), E(4,8), E(2,4), E(1,3), S(6), S(3));
    check_binary(Max(), E(-2,2), E(0,4), E(1,3), S(1), S(2));
    check_binary(Min(), E(-2,2), E(0,4), E(-1,1), S(0), S(2));

    check_power(E(-3,3), 2, E(1,4), S(-2));
    check_power(E(-3,3), 2, E(1,4), S(2));
    check_power(E(-3,3), 3, E(-9,-1), S(-2));
    check_power(E(1,4), -1, E(0.25_x,1.0_x), S(2));

    // Scalar-left overloads are distinct implementations.
    {
        UpperIntervalType a=U(S(3));
        backpropagate(U(S(5)),Add(),ValidatedNumber(2),a);
        ARIADNE_TESTING_PROVE(preserves(a,S(3)));
    }
    {
        UpperIntervalType a=U(S(3));
        backpropagate(U(S(-1)),Sub(),ValidatedNumber(2),a);
        ARIADNE_TESTING_PROVE(preserves(a,S(3)));
    }
    {
        UpperIntervalType a=U(S(3));
        backpropagate(U(S(6)),Mul(),ValidatedNumber(2),a);
        ARIADNE_TESTING_PROVE(preserves(a,S(3)));
    }
    {
        UpperIntervalType a=U(S(2));
        backpropagate(U(S(3)),Div(),ValidatedNumber(6),a);
        ARIADNE_TESTING_PROVE(preserves(a,S(2)));
    }
    {
        UpperIntervalType a=U(S(0));
        backpropagate(U(S(1)),Max(),ValidatedNumber(1),a);
        ARIADNE_TESTING_PROVE(preserves(a,S(0)));
    }
    {
        UpperIntervalType a=U(S(1));
        backpropagate(U(S(0)),Min(),ValidatedNumber(0),a);
        ARIADNE_TESTING_PROVE(preserves(a,S(1)));
    }

    // Scalar-right overloads are not used by ScalarProcedureInstruction today,
    // but remain part of the generic backpropagation API.
    {
        UpperIntervalType a=U(S(3));
        backpropagate(U(S(5)),Add(),a,ValidatedNumber(2));
        ARIADNE_TESTING_PROVE(preserves(a,S(3)));
    }
    {
        UpperIntervalType a=U(S(3));
        backpropagate(U(S(1)),Sub(),a,ValidatedNumber(2));
        ARIADNE_TESTING_PROVE(preserves(a,S(3)));
    }
    {
        UpperIntervalType a=U(S(3));
        backpropagate(U(S(6)),Mul(),a,ValidatedNumber(2));
        ARIADNE_TESTING_PROVE(preserves(a,S(3)));
    }
    {
        UpperIntervalType a=U(S(6));
        backpropagate(U(S(3)),Div(),a,ValidatedNumber(2));
        ARIADNE_TESTING_PROVE(preserves(a,S(6)));
    }
    {
        UpperIntervalType a=U(S(0));
        backpropagate(U(S(1)),Max(),a,ValidatedNumber(1));
        ARIADNE_TESTING_PROVE(preserves(a,S(0)));
    }
    {
        UpperIntervalType a=U(S(1));
        backpropagate(U(S(0)),Min(),a,ValidatedNumber(0));
        ARIADNE_TESTING_PROVE(preserves(a,S(1)));
    }
}

Void TestProcedure::test_backward_contractor_generated_witnesses()
{
    auto U = [](ExactIntervalType const& x) { return UpperIntervalType(x); };
    auto E = [](auto l, auto u) { return ExactIntervalType(l,u); };
    auto S = [](auto x) { return ExactIntervalType(x,x); };

    auto preserves = [](UpperIntervalType const& contracted, ExactIntervalType const& witness) {
        return !definitely(disjoint(contracted,witness));
    };

    auto check_unary = [&](auto op, ExactIntervalType domain, ExactIntervalType witness) {
        UpperIntervalType contracted=U(domain);
        UpperIntervalType witness_interval=U(witness);
        ARIADNE_TESTING_NOTE("[contractor] unary op=" << op
                  << " domain=" << domain
                  << " witness=" << witness);
        UpperIntervalType output=op(witness_interval);
        ARIADNE_TESTING_NOTE("[contractor] unary op=" << op
                  << " output=" << output);
        backpropagate(output,op,contracted);
        ARIADNE_TESTING_NOTE("[contractor] unary op=" << op
                  << " contracted=" << contracted);
        ARIADNE_TESTING_PROVE(preserves(contracted,witness));
    };

    auto check_binary = [&](auto op,
                            ExactIntervalType lhs_domain, ExactIntervalType rhs_domain,
                            ExactIntervalType lhs_witness, ExactIntervalType rhs_witness) {
        UpperIntervalType lhs=U(lhs_domain);
        UpperIntervalType rhs=U(rhs_domain);
        ARIADNE_TESTING_NOTE("[contractor] binary op=" << op
                  << " lhs_domain=" << lhs_domain
                  << " rhs_domain=" << rhs_domain
                  << " lhs_witness=" << lhs_witness
                  << " rhs_witness=" << rhs_witness);
        UpperIntervalType output=op(U(lhs_witness),U(rhs_witness));
        ARIADNE_TESTING_NOTE("[contractor] binary op=" << op
                  << " output=" << output);
        backpropagate(output,op,lhs,rhs);
        ARIADNE_TESTING_NOTE("[contractor] binary op=" << op
                  << " lhs_contracted=" << lhs
                  << " rhs_contracted=" << rhs);
        ARIADNE_TESTING_PROVE(preserves(lhs,lhs_witness));
        ARIADNE_TESTING_PROVE(preserves(rhs,rhs_witness));
    };

    auto check_power = [&](ExactIntervalType domain, ExactIntervalType witness, Int exponent) {
        UpperIntervalType contracted=U(domain);
        ARIADNE_TESTING_NOTE("[contractor] power exponent=" << exponent
                  << " domain=" << domain
                  << " witness=" << witness);
        UpperIntervalType output=pow(U(witness),exponent);
        ARIADNE_TESTING_NOTE("[contractor] power exponent=" << exponent
                  << " output=" << output);
        backpropagate(output,Pow(),contracted,exponent);
        ARIADNE_TESTING_NOTE("[contractor] power exponent=" << exponent
                  << " contracted=" << contracted);
        ARIADNE_TESTING_PROVE(preserves(contracted,witness));
    };

    // Generic unary grid over exactly representable witnesses.
    for (auto witness : {S(-2),S(-1),S(-0.5_x),S(0),S(0.5_x),S(1),S(2)}) {
        check_unary(Pos(),E(-3,3),witness);
        check_unary(Neg(),E(-3,3),witness);
        check_unary(Sqr(),E(-3,3),witness);
        check_unary(Exp(),E(-3,3),witness);
        check_unary(Sin(),E(-3,3),witness);
        check_unary(Cos(),E(-3,3),witness);
        check_unary(Tan(),E(-3,3),witness);
        check_unary(Atan(),E(-3,3),witness);
    }

    for (auto witness : {S(-2),S(-1),S(-0.5_x),S(0.5_x),S(1),S(2)}) {
        check_unary(Rec(),E(-3,3),witness);
        for (Int exponent=-3; exponent<=4; ++exponent) {
            check_power(E(-3,3),witness,exponent);
        }
    }

    for (auto witness : {S(0),S(0.5_x),S(1),S(2),S(4)}) {
        check_unary(Sqrt(),E(0,4),witness);
    }

    for (auto witness : {S(0.5_x),S(1),S(2)}) {
        check_unary(Log(),E(0.5_x,2),witness);
    }

    for (auto witness : {S(-1),S(-0.5_x),S(0),S(0.5_x),S(1)}) {
        check_unary(Asin(),E(-1,1),witness);
        check_unary(Acos(),E(-1,1),witness);
    }

    // Exhaustive small binary grid. Every output is evaluated from the
    // singleton witnesses, so each tested tuple is valid by construction.
    for (auto lhs_witness : {S(-2),S(-1),S(0),S(1),S(2)}) {
        for (auto rhs_witness : {S(-2),S(-1),S(0),S(1),S(2)}) {
            check_binary(Add(),E(-3,3),E(-3,3),lhs_witness,rhs_witness);
            check_binary(Sub(),E(-3,3),E(-3,3),lhs_witness,rhs_witness);
            check_binary(Mul(),E(-3,3),E(-3,3),lhs_witness,rhs_witness);
            check_binary(Max(),E(-3,3),E(-3,3),lhs_witness,rhs_witness);
            check_binary(Min(),E(-3,3),E(-3,3),lhs_witness,rhs_witness);
        }
    }

    for (auto lhs_witness : {S(-2),S(-1),S(0),S(1),S(2)}) {
        for (auto rhs_witness : {S(-2),S(-1),S(1),S(2)}) {
            check_binary(Div(),E(-3,3),E(-3,3),lhs_witness,rhs_witness);
        }
    }
}


Void TestProcedure::test_leq_backpropagate_soundness()
{
    // Leq is currently not representable by ProcedureInstruction, so this
    // overload is not on the simple_hull_reduce execution path.
    ARIADNE_TESTING_COMPILE_CONCEPT((not std::is_constructible_v<BinaryElementaryOperator,Leq>));
    ARIADNE_TESTING_COMPILE_CONCEPT((not std::is_constructible_v<ProcedureInstruction,Leq,SizeType,SizeType>));

    // Independently exercise the overload itself. The valid witness
    // a1=0, a2=1 satisfies a1<=a2 and must not be removed by propagation.
    FloatDPBounds r(0,dp);
    FloatDPBounds a1(0,0,dp);
    FloatDPBounds a2(1,1,dp);

    backpropagate(r,Leq(),a1,a2);

    ARIADNE_TESTING_PROVE(not inconsistent(a1,FloatDPBounds(0,dp)));
    ARIADNE_TESTING_PROVE(not inconsistent(a2,FloatDPBounds(1,dp)));
}


Void TestProcedure::test_inverse_trigonometric_empty_backpropagation()
{
    UpperIntervalType empty=ExactIntervalType::empty_interval();

    {
        UpperIntervalType argument=ExactIntervalType(0,0);
        backpropagate(empty,Asin(),argument);
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
    {
        UpperIntervalType argument=ExactIntervalType(0,0);
        backpropagate(empty,Acos(),argument);
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
    {
        UpperIntervalType argument=ExactIntervalType(0,0);
        backpropagate(empty,Atan(),argument);
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
}


Void TestProcedure::test_abs_backpropagation()
{
    {
        UpperIntervalType argument=ExactIntervalType(-2,2);
        UpperIntervalType result=ExactIntervalType(1,1);
        backpropagate(result,Abs(),argument);
        UpperIntervalType expected=ExactIntervalType(-1,1);
        ARIADNE_TESTING_EXPECT_EQUAL(argument.lower_bound().raw(), expected.lower_bound().raw());
        ARIADNE_TESTING_EXPECT_EQUAL(argument.upper_bound().raw(), expected.upper_bound().raw());
    }
    {
        UpperIntervalType argument=ExactIntervalType(-2,2);
        UpperIntervalType target=ExactIntervalType::empty_interval();
        backpropagate(target,Abs(),argument);
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
    {
        UpperIntervalType argument=ExactIntervalType(-2,2);
        UpperIntervalType target=ExactIntervalType(-2,-1);
        backpropagate(target,Abs(),argument);
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
}


Void TestProcedure::test_nul_hlf_backpropagation()
{
    {
        UpperIntervalType argument=ExactIntervalType(-4,4);
        UpperIntervalType result=ExactIntervalType(1,1);
        backpropagate(result,Hlf(),argument);
        UpperIntervalType expected=ExactIntervalType(2,2);
        ARIADNE_TESTING_EXPECT_EQUAL(argument.lower_bound().raw(), expected.lower_bound().raw());
        ARIADNE_TESTING_EXPECT_EQUAL(argument.upper_bound().raw(), expected.upper_bound().raw());
    }
    {
        UpperIntervalType argument=ExactIntervalType(-4,4);
        UpperIntervalType target=ExactIntervalType::empty_interval();
        backpropagate(target,Hlf(),argument);
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
    {
        UpperIntervalType argument=ExactIntervalType(-4,4);
        UpperIntervalType target=ExactIntervalType(-1,1);
        backpropagate(target,Nul(),argument);
        ARIADNE_TESTING_PROVE(not definitely(argument.is_empty()));
    }
    {
        UpperIntervalType argument=ExactIntervalType(-4,4);
        UpperIntervalType target=ExactIntervalType(1,2);
        backpropagate(target,Nul(),argument);
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
}


Void TestProcedure::test_elementary_empty_backpropagation()
{
    UpperIntervalType empty=ExactIntervalType::empty_interval();
    {
        UpperIntervalType argument=ExactIntervalType(-1,1);
        backpropagate(empty,Max(),argument,ValidatedNumber(0));
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
    {
        UpperIntervalType argument=ExactIntervalType(-1,1);
        backpropagate(empty,Min(),argument,ValidatedNumber(0));
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }

    {
        UpperIntervalType argument=ExactIntervalType(-1,1);
        backpropagate(empty,Sqrt(),argument);
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
    {
        UpperIntervalType argument=ExactIntervalType(-1,1);
        backpropagate(empty,Exp(),argument);
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
    {
        UpperIntervalType argument=ExactIntervalType(-1,1);
        backpropagate(empty,Log(),argument);
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
    {
        UpperIntervalType argument=ExactIntervalType(-1,1);
        backpropagate(empty,Tanh(),argument);
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
}


Void TestProcedure::test_derivative()
{
    typedef ApproximateTag P;
    typedef Number<P> Y;
    typedef FloatDPApproximation X;

    Y c(2);
    Formula<Y> o(Formula<Y>::constant(1));
    Formula<Y> x(Formula<Y>::coordinate(0));
    Formula<Y> y(Formula<Y>::coordinate(1));

    auto xs=sqr(x); auto ys=sqr(y);
    Formula<Y> e=sqrt(xs+c*ys)+sin(y/x)*xs*c-ys;
    ARIADNE_TESTING_INSPECT("e", e);
    ScalarMultivariateFunction<P> f(EuclideanDomain(2),e);
    ARIADNE_TESTING_INSPECT("f", f);
    Procedure<Y> p(f);
    ARIADNE_TESTING_INSPECT("p", p);

    X zero(0,dp);

    Vector<X> q({2,1},dp);
    Vector<X> s({-1,3},dp);
    Vector<Differential<X>> ds=Differential<X>::variable(1,2,zero,0)*s+q;
    // Compare each gradient component against the original absolute tolerance.
    const auto procedure_gradient=gradient(p,q);
    const auto function_gradient=f.gradient(q);
    ARIADNE_TESTING_EXPECT_RANGES_BY(procedure_gradient.array(), function_gradient.array(),
        [](auto const& actual, auto const& expected) {
            return abs(actual-expected).raw() <= FloatDP(ExactDouble(8e-16),dp);
        });
    // The two evaluation paths may differ by roundoff; compare their
    // approximations using the absolute Hessian tolerance already used below.
    const auto procedure_hessian=hessian(p,q,s);
    const auto differential_result=f(ds);
    const auto function_hessian=differential_result.hessian().get(0,0);
    const auto derivative_hessian_error=abs(procedure_hessian-function_hessian);
    ARIADNE_TESTING_INSPECT("procedure_hessian", procedure_hessian);
    ARIADNE_TESTING_INSPECT("function_hessian", function_hessian);
    ARIADNE_TESTING_EXPECT_POSSIBLY(procedure_hessian==function_hessian);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(derivative_hessian_error.raw(), FloatDP(ExactDouble(1e-14),dp));
}

Void TestProcedure::test_scalar_formula_conversion() {
    using F=ApproximateFormula;
    auto x=F::coordinate(0u),y=F::coordinate(1u);
    auto formula=F::scalar(Mul(),ApproximateNumber(3),x+y);
    ApproximateProcedure procedure(2u,formula);
    Vector<FloatDPApproximation> point({2,4},dp);
    ARIADNE_TESTING_EXPECT_SAME(evaluate(procedure,point), FloatDPApproximation(18,dp));
    Vector<ApproximateProcedure> vector(procedure);
    ARIADNE_TESTING_EXPECT_EQUAL(vector.argument_size(), 2u);
    const auto converted_vector_value=evaluate(vector,point);
    ARIADNE_TESTING_EXPECT_SAME(converted_vector_value[0u], FloatDPApproximation(18,dp));
}


class TestProcedureExtra : public Ariadne::Testing::ContextFixture {
public:
    using ContextFixture::ContextFixture;
    void test_function_procedure_construction();
    void test_procedure_interval_contracts_and_empty_targets();
    void test_procedure_instruction_accessors();
    void test_procedure_builders_and_contracts();
    void test_procedure_scalar_derivatives();
};

Void TestProcedureExtra::test_procedure_instruction_accessors() {
    ProcedureInstruction constant(Cnst(),2u),variable(Var(),1u),unary(Neg(),3u);
    ProcedureInstruction binary(Add(),4u,5u),power(Pow(),6u,3);
    ProcedureInstruction scalar(ProcedureInstructionVariant(ScalarProcedureInstruction(Mul(),0u,1u)));
    ARIADNE_TESTING_EXPECT_EQUAL(constant.val(), 2u);ARIADNE_TESTING_EXPECT_EQUAL(variable.ind(), 1u);
    ARIADNE_TESTING_EXPECT_EQUAL(unary.arg(), 3u);ARIADNE_TESTING_EXPECT_EQUAL(power.arg(), 6u);
    ARIADNE_TESTING_EXPECT_EQUAL(binary.arg1(), 4u);ARIADNE_TESTING_EXPECT_EQUAL(binary.arg2(), 5u);
    ARIADNE_TESTING_EXPECT_EQUAL(power.num(), 3);
    ARIADNE_TESTING_PROVE(constant.op().code()==OperatorCode::CNST);
    ARIADNE_TESTING_PROVE(variable.op().code()==OperatorCode::VAR);
    ARIADNE_TESTING_PROVE(unary.op().code()==OperatorCode::NEG);
    ARIADNE_TESTING_PROVE(binary.op().code()==OperatorCode::ADD);
    ARIADNE_TESTING_PROVE(power.op().code()==OperatorCode::POW);
    ARIADNE_TESTING_PROVE(scalar.op().code()==OperatorCode::MUL);
}
Void TestProcedureExtra::test_procedure_builders_and_contracts() {
    using X=FloatDPApproximation;using F=ApproximateFormula;
    ApproximateProcedure constant(0u);constant.new_constant(ApproximateNumber(7));constant.new_instruction(Cnst(),0u);
    ARIADNE_TESTING_EXPECT_SAME(evaluate(constant,Vector<X>(0u,dp)), X(7,dp));
    auto x=F::coordinate(0u),y=F::coordinate(1u);
    ApproximateProcedure sum(2u,x+y);
    Vector<ApproximateProcedure> vector(sum);
    ARIADNE_TESTING_EXPECT_EQUAL(vector.argument_size(), 2u);
    auto previous=vector.temporaries_size()-1u;
    vector.new_instruction(Neg(),previous);
    vector.new_instruction(Pow(),vector.temporaries_size()-1u,2);
    vector.new_instruction(Add(),previous,vector.temporaries_size()-1u);
    vector.set_return(0u,vector.temporaries_size()-1u);
    Vector<X> point({2,4},dp);
    const auto constructed_vector_value=evaluate(vector,point);
    ARIADNE_TESTING_EXPECT_SAME(constructed_vector_value[0u], X(42,dp));
    List<X> temporaries(vector.temporaries_size(),X(0,dp));
    execute(temporaries,vector,point);
    ARIADNE_TESTING_EXPECT_SAME(temporaries.back(), X(42,dp));
    ARIADNE_TESTING_EXPECT_EXCEPTION(evaluate(sum,Vector<X>(1u,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(evaluate(vector,Vector<X>(1u,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(execute(temporaries,vector,Vector<X>(1u,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(execute(temporaries,sum,point));
    ARIADNE_TESTING_EXPECT_EXCEPTION(gradient(sum,Vector<X>(1u,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(hessian(sum,point,Vector<X>(1u,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(hessian(sum,Vector<X>(1u,dp),point));
    ApproximateProcedure empty(2u);
    ARIADNE_TESTING_EXPECT_EXCEPTION(evaluate(empty,point));
    ARIADNE_TESTING_EXPECT_EXCEPTION(gradient(empty,point));
    using I=UpperIntervalType;
    ValidatedProcedure identity(1u,ValidatedFormula::coordinate(0u));
    Vector<I> interval{I(-2,2)};List<I> values(1u,I(1,1));
    backpropagate(values,identity,interval);
    ARIADNE_TESTING_EXPECT_SAME(interval[0u], I(1,1));
    Vector<ValidatedProcedure> vector_identity(identity);
    interval[0u]=I(-2,2);backpropagate(values,vector_identity,interval);
    ARIADNE_TESTING_EXPECT_SAME(interval[0u], I(1,1));
    List<I> wrong_values;Vector<I> wrong_point(0u,I(0,0));
    ARIADNE_TESTING_EXPECT_EXCEPTION(backpropagate(wrong_values,identity,interval));
    ARIADNE_TESTING_EXPECT_EXCEPTION(backpropagate(values,identity,wrong_point));
    ARIADNE_TESTING_EXPECT_EXCEPTION(backpropagate(values,vector_identity,wrong_point));
    ARIADNE_TESTING_EXPECT_EXCEPTION(execute(values,identity,wrong_point));
}
Void TestProcedureExtra::test_procedure_scalar_derivatives() {
    using X=FloatDPApproximation;using F=ApproximateFormula;
    Vector<X> point{X(2,dp)},direction{X(1,dp)};
    struct Case { BinaryElementaryOperator op;Int value;ExactDouble first,second; };
    List<Case> cases{{Add(),5,1.0_x,0.0_x},{Sub(),1,-1.0_x,0.0_x},{Mul(),6,3.0_x,0.0_x},{Div(),0,-0.75_x,0.75_x}};
    for(const auto& example:cases) {
        auto expression=F::scalar(example.op,ApproximateNumber(3),F::coordinate(0u));
        ApproximateProcedure p(1u,expression);
        if(example.op.code()==OperatorCode::DIV) { ARIADNE_TESTING_EXPECT_SAME(evaluate(p,point), X(1.5_x,dp)); }
        else { ARIADNE_TESTING_EXPECT_SAME(evaluate(p,point), X(example.value,dp)); }
        const auto procedure_gradient=gradient(p,point);
        ARIADNE_TESTING_EXPECT_SAME(procedure_gradient[0u], X(example.first,dp));
        ARIADNE_TESTING_EXPECT_SAME(hessian(p,point,direction), X(example.second,dp));
    }
    ApproximateProcedure cube(1u,F::graded(Pow(),F::coordinate(0u),3));
    const auto cube_gradient=gradient(cube,point);
    ARIADNE_TESTING_EXPECT_SAME(cube_gradient[0u], X(12,dp));
    ApproximateProcedure polynomial(1u,pow(F::coordinate(0u),3)+ApproximateNumber(7));
    const auto polynomial_gradient=gradient(polynomial,point);
    ARIADNE_TESTING_EXPECT_SAME(polynomial_gradient[0u], X(12,dp));
    for(const auto& op:List<UnaryElementaryOperator>{Asin(),Acos()}) {
        ApproximateProcedure p(1u,F::unary(op,F::coordinate(0u)));
        Vector<X> argument{X(0.5_x,dp)},velocity{X(2,dp)};
        X expected=X(2,dp)/(X(0.75_x,dp)*sqrt(X(0.75_x,dp)));
        if(op.code()==OperatorCode::ACOS) { expected=-expected; }
        const auto computed_hessian=hessian(p,argument,velocity);
        const auto hessian_error=abs(computed_hessian-expected);
        ARIADNE_TESTING_INSPECT("computed_hessian", computed_hessian);
        ARIADNE_TESTING_INSPECT("expected_hessian", expected);
        ARIADNE_TESTING_EXPECT_LESS_EQUAL(hessian_error.raw(), FloatDP(ExactDouble(1e-14),dp));
    }
    for(const auto& op:List<BinaryElementaryOperator>{Max(),Min()}) {
        ApproximateProcedure p(1u,F::scalar(op,ApproximateNumber(3),F::coordinate(0u)));
        ARIADNE_TESTING_EXPECT_EXCEPTION(gradient(p,point));ARIADNE_TESTING_EXPECT_EXCEPTION(hessian(p,point,direction));
        ApproximateProcedure binary(1u,F::binary(op,F::coordinate(0u),F(ApproximateNumber(3))));
        ARIADNE_TESTING_EXPECT_EXCEPTION(gradient(binary,point));
    }
    ApproximateProcedure absolute(1u,abs(F::coordinate(0u)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(gradient(absolute,point));
    ValidatedProcedure validated(1u,ValidatedFormula::graded(Pow(),ValidatedFormula::coordinate(0u),3));
    Vector<FloatDPBounds> bounds{FloatDPBounds(2,dp)};
    ARIADNE_TESTING_PROVE(definitely(evaluate(validated,bounds)==FloatDPBounds(8,dp)));
    ValidatedProcedure scalar_validated(1u,ValidatedFormula::scalar(Mul(),ValidatedNumber(3),ValidatedFormula::coordinate(0u)));
    ARIADNE_TESTING_PROVE(definitely(evaluate(scalar_validated,bounds)==FloatDPBounds(6,dp)));
}


Void TestProcedureExtra::test_procedure_interval_contracts_and_empty_targets() {
    using B=FloatDPBounds;using D=UnivariateSecondDifferential<FloatDPApproximation>;
    ValidatedProcedure identity(1u,ValidatedFormula::coordinate(0u)),empty(1u);
    Vector<B> point{B(2,dp)};List<B> values(1u,B(0,dp));
    execute(values,identity,point);ARIADNE_TESTING_PROVE(definitely(values[0u]==B(2,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(execute(values,identity,Vector<B>(0u,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(evaluate(identity,Vector<B>(0u,dp)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(evaluate(empty,point));
    Vector<ValidatedProcedure> vector(identity);UpperBoxType domain{{-1,1}};
    ARIADNE_TESTING_EXPECT_EXCEPTION(simple_hull_reduce(domain,vector,ExactBoxType(0u)));
    ApproximateProcedure p(1u,ApproximateFormula::coordinate(0u)),no_instructions(1u);
    Vector<D> differentials(1u,D(FloatDPApproximation(2,dp)));
    const auto differential_value=evaluate(p,differentials);
    ARIADNE_TESTING_EXPECT_SAME(differential_value.value(), FloatDPApproximation(2,dp));
    ARIADNE_TESTING_EXPECT_EXCEPTION(evaluate(p,Vector<D>(0u,D(FloatDPApproximation(0,dp)))));
    ARIADNE_TESTING_EXPECT_EXCEPTION(evaluate(no_instructions,differentials));
    ARIADNE_TESTING_PROVE(definitely(gradient(identity,point)[0u]==B(1,dp)));
    ARIADNE_TESTING_PROVE(definitely(hessian(identity,point,Vector<B>{B(1,dp)})==B(0,dp)));
    using I=UpperIntervalType;
    I target=I::empty_interval();
    for(const auto& op:List<BinaryElementaryOperator>{Max(),Min()}) {
        I a(-2,2),b(-3,3);backpropagate(target,op,a,b);
        ARIADNE_TESTING_PROVE(definitely(a.is_empty()));ARIADNE_TESTING_PROVE(definitely(b.is_empty()));
        a=I(-2,2);const ValidatedNumber coefficient(1);backpropagate(target,op,coefficient,a);
        ARIADNE_TESTING_PROVE(definitely(a.is_empty()));
    }
    for(const auto& op:List<UnaryElementaryOperator>{Nul(),Cos(),Tan()}) {
        I argument(-2,2);backpropagate(target,op,argument);
        ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    }
    I argument(-2,2);backpropagate(I(-2,-1),Nul(),argument);
    ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
    argument=I(-2,2);backpropagate(target,Pow(),argument,3);
    ARIADNE_TESTING_PROVE(definitely(argument.is_empty()));
}
Void TestProcedureExtra::test_function_procedure_construction() {
    const auto x=ValidatedScalarMultivariateFunction::coordinate(2u,0u);
    const auto y=ValidatedScalarMultivariateFunction::coordinate(2u,1u);
    const auto scalar=make_procedure(sqr(x)+y);
    const Vector<ValidatedProcedure> vector(ValidatedVectorMultivariateFunction{sqr(x)+y,x-y});
    const Vector<FloatDPBounds> point({2,3},dp);
    ARIADNE_TESTING_EXPECT_EQUAL(scalar.argument_size(), 2u);
    ARIADNE_TESTING_EXPECT_EQUAL(evaluate(scalar,point), FloatDPBounds(7,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(vector.result_size(), 2u);
    ARIADNE_TESTING_EXPECT_EQUAL(evaluate(vector,point), Vector<FloatDPBounds>({7,-1},dp));
}

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_procedure");
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_construct_from_formula);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_construct_from_expansion);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_evaluate);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_propagate);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_backward_contractor_soundness);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_backward_contractor_witness_preservation);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_backward_contractor_generated_witnesses);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_leq_backpropagate_soundness);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_inverse_trigonometric_empty_backpropagation);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_abs_backpropagation);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_nul_hlf_backpropagation);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_elementary_empty_backpropagation);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_derivative);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_scalar_formula_conversion);
    ARIADNE_TESTING_ADD_TEST(TestProcedure, test_unsupported_graded_operations);
    ARIADNE_TESTING_ADD_TEST(TestProcedureExtra, test_function_procedure_construction);
    ARIADNE_TESTING_ADD_TEST(TestProcedureExtra, test_procedure_interval_contracts_and_empty_targets);
    ARIADNE_TESTING_ADD_TEST(TestProcedureExtra, test_procedure_instruction_accessors);
    ARIADNE_TESTING_ADD_TEST(TestProcedureExtra, test_procedure_builders_and_contracts);
    ARIADNE_TESTING_ADD_TEST(TestProcedureExtra, test_procedure_scalar_derivatives);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
