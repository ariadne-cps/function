/***************************************************************************
 *            function/procedure.hpp
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

/*! @file function/procedure.hpp
 *  \brief Procedure to compute a real function
 */

#ifndef ARIADNE_PROCEDURE_HPP
#define ARIADNE_PROCEDURE_HPP

#include <iostream>

#include "algebra/vector.hpp"
#include "function/formula.decl.hpp"
#include "function/function.decl.hpp"
#include "numeric/numeric.hpp"
#include "numeric/operators.hpp"
#include "utility/attribute.hpp"
#include "utility/container.hpp"

namespace Ariadne {

class MultiIndex;
template<class I, class X>
class Expansion;
template<class X>
class Graded;

template<class Y>
class Procedure;
//!@{
//! \ingroup FunctionModule
//! \relates Procedure
//! \name Type synonyms
using ApproximateProcedure = Procedure<ApproximateNumber>; //!< <p/>
using ValidatedProcedure = Procedure<ValidatedNumber>;     //!< <p/>
using EffectiveProcedure = Procedure<EffectiveNumber>;     //!< <p/>
//!@}

/** @brief Counters recording the work performed by hull reduction. */
struct ProcedureHullReductionStatistics {
    double temporary_allocation_seconds = 0.0;
    double forward_execution_seconds = 0.0;
    double backward_propagation_seconds = 0.0;
};

/** @brief Contract the argument hull using the procedure constraints. */
Void simple_hull_reduce(UpperBoxType& dom, const ValidatedProcedure& f, IntervalDomainType codom);
/** @brief Contract the argument hull using the procedure constraints. */
Void simple_hull_reduce(
    UpperBoxType& dom, const ValidatedProcedure& f, IntervalDomainType codom,
    ProcedureHullReductionStatistics& statistics);
/** @brief Contract the argument hull using the procedure constraints. */
Void simple_hull_reduce(UpperBoxType& dom, const Vector<ValidatedProcedure>& f, BoxDomainType codom);

/** @brief Procedure instruction loading a constant. */
struct ConstantProcedureInstruction {
    Cnst _op;
    SizeType _val;
    /** @brief Construct a \c ConstantProcedureInstruction from op, val. */
    ConstantProcedureInstruction(Cnst op, SizeType val) : _op(op), _val(val) {}
};
/** @brief Procedure instruction loading an argument coordinate. */
struct IndexProcedureInstruction {
    Var _op;
    SizeType _ind;
    /** @brief Construct a \c IndexProcedureInstruction from op, ind. */
    IndexProcedureInstruction(Var op, SizeType ind) : _op(op), _ind(ind) {}
};
/** @brief Procedure instruction applying a unary operation. */
struct UnaryProcedureInstruction {
    UnaryElementaryOperator _op;
    SizeType _arg;
    /** @brief Construct a \c UnaryProcedureInstruction from op, arg. */
    UnaryProcedureInstruction(UnaryElementaryOperator op, SizeType arg) : _op(op), _arg(arg) {}
};
/** @brief Procedure instruction applying a binary operation. */
struct BinaryProcedureInstruction {
    BinaryElementaryOperator _op;
    SizeType _arg1, _arg2;
    /** @brief Construct a \c BinaryProcedureInstruction from op, arg1, arg2. */
    BinaryProcedureInstruction(BinaryElementaryOperator op, SizeType arg1, SizeType arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) {}
};
/** @brief Procedure instruction applying an integer-graded operation. */
struct GradedProcedureInstruction {
    GradedElementaryOperator _op;
    SizeType _arg;
    Int _num;
    /** @brief Construct a \c GradedProcedureInstruction from op, arg, num. */
    GradedProcedureInstruction(GradedElementaryOperator op, SizeType arg, Int num)
        : _op(op), _arg(arg), _num(num) {}
};
/** @brief Procedure instruction combining a temporary with a scalar constant. */
struct ScalarProcedureInstruction {
    BinaryElementaryOperator _op;
    SizeType _arg1, _arg2;
    /** @brief Construct a \c ScalarProcedureInstruction from op, arg1, arg2. */
    ScalarProcedureInstruction(BinaryElementaryOperator op, SizeType arg1, SizeType arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) {}
};

typedef Variant<ConstantProcedureInstruction, IndexProcedureInstruction,
                UnaryProcedureInstruction, BinaryProcedureInstruction, GradedProcedureInstruction, ScalarProcedureInstruction>
    ProcedureInstructionVariant;

/** @brief Tagged instruction in an executable arithmetic procedure. */
struct ProcedureInstruction : public ProcedureInstructionVariant {
    /** @brief Construct a \c ProcedureInstruction from o, c. */
    explicit ProcedureInstruction(Cnst o, SizeType c)
        : ProcedureInstructionVariant(ConstantProcedureInstruction(o, c)) {}
    /** @brief Construct a \c ProcedureInstruction from o, i. */
    explicit ProcedureInstruction(Var o, SizeType i)
        : ProcedureInstructionVariant(IndexProcedureInstruction(o, i)) {}
    /** @brief Construct a \c ProcedureInstruction from o, a. */
    explicit ProcedureInstruction(UnaryElementaryOperator o, SizeType a)
        : ProcedureInstructionVariant(UnaryProcedureInstruction(o, a)) {}
    /** @brief Construct a \c ProcedureInstruction from o, a1, a2. */
    explicit ProcedureInstruction(BinaryElementaryOperator o, SizeType a1, SizeType a2)
        : ProcedureInstructionVariant(BinaryProcedureInstruction(o, a1, a2)) {}
    /** @brief Construct a \c ProcedureInstruction from o, a, n. */
    explicit ProcedureInstruction(GradedElementaryOperator o, SizeType a, Int n)
        : ProcedureInstructionVariant(GradedProcedureInstruction(o, a, n)) {}

  public:
    /** @brief Construct a \c ProcedureInstruction from var. */
    ProcedureInstruction(ProcedureInstructionVariant var) : ProcedureInstructionVariant(var) {}
    /** @brief Return the underlying base object. */
    const ProcedureInstructionVariant& base() const { return *this; }
    /** @brief Dispatch the stored alternative to the visitor. */
    template<class VIS>
    decltype(auto) accept(VIS&& vis) const {
        return std::visit(std::forward<VIS>(vis), static_cast<const ProcedureInstructionVariant&>(*this));
    }
    /** @brief Return the stored operation. */
    Operator op() const {
        return this->accept([](auto s) { return Operator(s._op); });
    }
    /** @brief Return the stored value. */
    const SizeType& val() const {
        return std::get<ConstantProcedureInstruction>(this->base())._val;
    }
    /** @brief Return the coordinate index. */
    const SizeType& ind() const {
        return std::get<IndexProcedureInstruction>(this->base())._ind;
    }
    /** @brief Return the operand formula. */
    const SizeType& arg() const {
        if (auto uptr = std::get_if<UnaryProcedureInstruction>(&this->base())) {
            return uptr->_arg;
        } else {
            return std::get<GradedProcedureInstruction>(this->base())._arg;
        }
    }
    /** @brief Return the first operand. */
    const SizeType& arg1() const { return std::get<BinaryProcedureInstruction>(this->base())._arg1; }
    /** @brief Return the second operand. */
    const SizeType& arg2() const { return std::get<BinaryProcedureInstruction>(this->base())._arg2; }
    /** @brief Return the stored numeric operand. */
    const Int& num() const { return std::get<GradedProcedureInstruction>(this->base())._num; }

  public:
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ProcedureInstruction& pri);
};

//! \ingroup FunctionModule
//! \brief An algorithmic procedure for computing a function.
//!
//! A Procedure is more efficient to compute than a Formula, since common
//! subexpressions have already been removed. However, it is also more
//! difficult to manipulate, so it should usually only be used when no further
//! manipulations of the function are possible.
template<class Y>
class Procedure {
    typedef Paradigm<Y> P;

  public:
    /** @brief Construct a \c Procedure from as. */
    explicit Procedure(SizeType as);
    /** @brief Construct a \c Procedure from f. */
    explicit Procedure(const ScalarMultivariateFunction<P>& f);
    /** @brief Construct a \c Procedure from as, f. */
    explicit Procedure(SizeType as, const Formula<Y>& f);
    /** @brief Construct a \c Procedure from e. */
    template<ConvertibleTo<Y> X>
    explicit Procedure(const Expansion<MultiIndex, X>& e);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Procedure<Y>& p) { return p._write(os); }

  public:
    /** @brief Return the number of scalar arguments. */
    SizeType argument_size() const;
    /** @brief Evaluate at the supplied argument. */
    template<class X, class YY>
    friend X evaluate(const Procedure<YY>& p, const Vector<X>& x);

  private:

  public:
    SizeType _argument_size;
    List<Y> _constants;
    List<ProcedureInstruction> _instructions;

  public:
    /** @brief Insert a constant into the procedure and return its index. */
    Void new_constant(const Y& c) { _constants.append(c); }
    /** @brief Append a procedure instruction and return its result index. */
    Void new_instruction(Cnst o, SizeType c) { _instructions.append(ProcedureInstruction(o, c)); }
    /** @brief Append a procedure instruction and return its result index. */
    Void new_instruction(Var o, SizeType i) { _instructions.append(ProcedureInstruction(o, i)); }
    /** @brief Append a procedure instruction and return its result index. */
    Void new_instruction(UnaryElementaryOperator o, SizeType a) { _instructions.append(ProcedureInstruction(o, a)); }
    /** @brief Append a procedure instruction and return its result index. */
    Void new_instruction(BinaryElementaryOperator o, SizeType a1, SizeType a2) { _instructions.append(ProcedureInstruction(o, a1, a2)); }
    /** @brief Append an instruction combining a temporary value with a scalar. */
    Void new_instruction_scalar(BinaryElementaryOperator o, SizeType c1, SizeType a2) {
        _instructions.append(ProcedureInstruction(ProcedureInstructionVariant(ScalarProcedureInstruction(o, c1, a2))));
    }
    /** @brief Append a procedure instruction and return its result index. */
    Void new_instruction(GradedElementaryOperator o, SizeType a, Int n) { _instructions.append(ProcedureInstruction(o, a, n)); }

  private:
    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream& os) const;
};

// These exported specializations are emitted by procedure.cpp.
extern template class Procedure<ApproximateNumber>;
extern template class Procedure<ValidatedNumber>;

//! \brief An algorithmic procedure for computing a function.
//!
//! A Procedure is more efficient to compute than a Formula, since common
//! subexpressions have already been removed. However, it is also more
//! difficult to manipulate, so it should usually only be used when no further
//! manipulations of the function are possible.
template<class Y>
class Vector<Procedure<Y>> {
    typedef Paradigm<Y> P;

  public:
    /** @brief Construct a \c Vector from as, f. */
    explicit Vector(SizeType as, const Vector<Formula<Y>>& f);
    /** @brief Copy-construct the object from the supplied instance. */
    explicit Vector(const VectorMultivariateFunction<P>& f);
    /** @brief Construct a \c Vector from p. */
    explicit Vector(const Procedure<Y>& p);
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Vector<Procedure<Y>>& p) { return p._write(os); }

  public:
    /** @brief Return the number of result components. */
    SizeType result_size() const { return _results.size(); }
    /** @brief Return the number of temporary slots required for procedure execution. */
    SizeType temporaries_size() const { return _instructions.size(); }
    /** @brief Return the number of scalar arguments. */
    SizeType argument_size() const { return this->_argument_size; }
    /** @brief Append a procedure instruction and return its result index. */
    Void new_instruction(UnaryElementaryOperator o, SizeType a) { _instructions.append(ProcedureInstruction(o, a)); }
    /** @brief Append a procedure instruction and return its result index. */
    Void new_instruction(GradedElementaryOperator o, SizeType a, Int n) { _instructions.append(ProcedureInstruction(o, a, n)); }
    /** @brief Append a procedure instruction and return its result index. */
    Void new_instruction(BinaryElementaryOperator o, SizeType a1, SizeType a2) { _instructions.append(ProcedureInstruction(o, a1, a2)); }
    /** @brief Set the return. */
    Void set_return(SizeType i, SizeType a) { _results[i] = a; }

  public:
    SizeType _argument_size;
    List<Y> _constants;
    List<ProcedureInstruction> _instructions;
    Vector<SizeType> _results;

  private:
    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream& os) const;
};

template<class X, class Y>
Void _execute(List<X>& v, const List<ProcedureInstruction>& p, const List<Y>& c, const Vector<X>& x);
template<class X, class Y>
Void _backpropagate(Vector<X>& x, List<X>& v, const List<ProcedureInstruction>& p, const List<Y>& c);

// \related Procedure \brief Evaluate a function \a p defined by an algorithmic procedure.
/** @brief Evaluate at the supplied argument. */
template<class X, class Y>
inline X evaluate(const Procedure<Y>& p, const Vector<X>& x) {
    ARIADNE_PRECONDITION(p.argument_size() == x.size());
    List<X> t(p._instructions.size(), x.zero_element());
    _execute(t, p._instructions, p._constants, x);
    ARIADNE_PRECONDITION(!t.empty());
    return std::move(t.back());
}

// \related Procedure \brief Evaluate a function \a p defined by an algorithmic procedure.
/** @brief Evaluate at the supplied argument. */
template<class X, class Y>
inline Vector<X> evaluate(const Vector<Procedure<Y>>& p, const Vector<X>& x) {
    ARIADNE_PRECONDITION(p.argument_size() == x.size());
    List<X> t(p._instructions.size(), x.zero_element());
    _execute(t, p._instructions, p._constants, x);
    Vector<X> r(p.result_size(), x.zero_element());
    for (SizeType i = 0; i != r.size(); ++i) {
        r[i] = std::move(t[p._results[i]]);
    }
    return r;
}

// \related Procedure \brief Execute the procedure \a p storing the results in \a t.
/** @brief Execute the procedure using the supplied arguments and workspace. */
template<class X, class Y>
Void execute(List<X>& t, const Procedure<Y>& p, const Vector<X>& x) {
    ARIADNE_PRECONDITION(p.argument_size() == x.size());
    _execute(t, p._instructions, p._constants, x);
}

/** @brief Execute the procedure using the supplied arguments and workspace. */
template<class X, class Y>
Void execute(List<X>& t, const Vector<Procedure<Y>>& p, const Vector<X>& x) {
    ARIADNE_PRECONDITION(p.argument_size() == x.size());
    _execute(t, p._instructions, p._constants, x);
}

// \related Procedure \brief Backpropagate the results of a validated procedure to the inputs \a x.
/** @brief Propagate output constraints backward through the operation. */
template<class X, class Y>
Void backpropagate(List<X>& t, const Procedure<Y>& p, Vector<X>& x) {
    ARIADNE_PRECONDITION(p.argument_size() == x.size());
    _backpropagate(x, t, p._instructions, p._constants);
}

/** @brief Propagate output constraints backward through the operation. */
template<class X, class Y>
Void backpropagate(List<X>& t, const Vector<Procedure<Y>>& p, Vector<X>& x) {
    ARIADNE_PRECONDITION(p.argument_size() == x.size());
    _backpropagate(x, t, p._instructions, p._constants);
}

// \related Procedure \brief Compute the gradient of a Procedure by backwards automatic differentiation.
template<class X, class Y>
Covector<X> gradient(const Procedure<Y>& f, const Vector<X>& x);
// \related Procedure \brief Compute the second derivative of a Procedure at \a x in direction \a s.
template<class X, class Y>
X hessian(const Procedure<Y>& f, const Vector<X>& x, const Vector<X>& s);

// \related Convert a function into a procedure.
template<class P>
Procedure<Number<P>> make_procedure(const ScalarMultivariateFunction<P>& f);

} // namespace Ariadne

#endif // ARIADNE_PROCEDURE_HPP
