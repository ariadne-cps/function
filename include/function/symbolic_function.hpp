/***************************************************************************
 *            function/symbolic_function.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! @file function/symbolic_function.hpp
 *  \brief Symbolic functions
 */

#ifndef ARIADNE_SYMBOLIC_FUNCTION_HPP
#define ARIADNE_SYMBOLIC_FUNCTION_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "algebra/vector.hpp"
#include "function/formula.hpp"
#include "function/formula.tpl.hpp"
#include "function/function_interface.hpp"
#include "function/function.hpp"
#include "function/function_mixin.hpp"
#include "function/function_wrapper.hpp"
#include "function/projection.hpp"
#include "numeric/numeric.hpp"
#include "numeric/operators.tpl.hpp"
#include "utility/container.hpp"
#include "utility/macros.hpp"
#include "utility/metaprogramming.hpp"
#include "utility/pointer.hpp"
#include "utility/stlio.hpp"

namespace Ariadne {

//------------------------ Formula functions  -----------------------------------//

//! A function defined by a formula
template<class Y>
class ScalarUnivariateFormulaFunction
//     : public ScalarUnivariateFunctionMixin<ScalarUnivariateFormulaFunction<Y>,InformationTag<Y>>
{
    typedef InformationTag<Y> P;
    Formula<Y> _formula;

  public:
    /** @brief Construct a \c ScalarUnivariateFormulaFunction from f. */
    ScalarUnivariateFormulaFunction(const Formula<Y>& f) : _formula(f) {}
    /** @brief Return the scalar formula. */
    Formula<Y> formula() const { return _formula; }
    /** @brief Convert to formula< y >. */
    operator Formula<Y>() const { return _formula; }

    /** @brief Return the domain. */
    RealDomain domain() const { return RealDomain(); }
    /** @brief Return the number of scalar arguments. */
    SizeOne argument_size() const { return SizeOne(); }
    /** @brief Return the number of result components. */
    SizeOne result_size() const { return SizeOne(); }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    X operator()(const X& x) const { return Ariadne::evaluate(_formula, Vector<X>({x})); }
    /** @brief Return the derivative with respect to the selected argument. */
    friend ScalarUnivariateFormulaFunction<Y> derivative(const ScalarUnivariateFormulaFunction<Y>& f, IndexZero) {
        return ScalarUnivariateFormulaFunction<Y>(Ariadne::derivative(f._formula, 0));
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ScalarUnivariateFormulaFunction<Y>& f) {
        return Formula<Y>::_write_univariate(os, f);
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Representation<ScalarUnivariateFormulaFunction<Y>>& fr) {
        return os << "ScalarUnivariateFormulaFunction(" << fr.reference() << ")";
    }
};

typedef ScalarUnivariateFormulaFunction<EffectiveNumber> EffectiveScalarUnivariateFormulaFunction;

//! A function defined by a formula
template<class Y>
class VectorUnivariateFormulaFunction {
    typedef InformationTag<Y> P;
    Vector<Formula<Y>> _formulae;

  public:
    /** @brief Construct a \c VectorUnivariateFormulaFunction from f. */
    VectorUnivariateFormulaFunction(const List<Formula<Y>>& f) : _formulae(f) {}
    /** @brief Construct a \c VectorUnivariateFormulaFunction from f. */
    VectorUnivariateFormulaFunction(const Vector<Formula<Y>>& f) : _formulae(f) {}

    /** @brief Return the component formulae. */
    Vector<Formula<Y>> formulae() const { return _formulae; }
    /** @brief Return the selected element or component. */
    ScalarUnivariateFormulaFunction<Y> operator[](SizeType i) const { return ScalarUnivariateFormulaFunction(this->_formulae[i]); }

    /** @brief Return the domain. */
    RealDomain domain() const { return RealDomain(); }
    /** @brief Return the number of result components. */
    SizeType result_size() const { return this->_formulae.size(); }
    /** @brief Return the number of scalar arguments. */
    SizeOne argument_size() const { return SizeOne(); }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    Vector<X> operator()(const X& x) const { return Ariadne::evaluate(this->_formulae, Vector<X>({x})); }
    /** @brief Return the derivative with respect to the selected argument. */
    friend VectorUnivariateFormulaFunction<Y> derivative(const VectorUnivariateFormulaFunction<Y>& f, IndexZero j) {
        return VectorUnivariateFormulaFunction<Y>(Vector<Formula<Y>>(f._formulae.size(), [&](SizeType i) { return derivative(f._formulae[i], j); }));
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const VectorUnivariateFormulaFunction<Y>& f) {
        return Formula<Y>::_write_univariate(os, f._formulae);
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Representation<VectorUnivariateFormulaFunction<Y>>& fr) {
        return os << "VectorUnivariateFormulaFunction(" << fr.reference() << ")";
    }
};

typedef VectorUnivariateFormulaFunction<EffectiveNumber> EffectiveVectorUnivariateFormulaFunction;

//! A function defined by a formula
template<class Y>
class ScalarMultivariateFormulaFunction {
    typedef InformationTag<Y> P;
    SizeType _argument_size;
    Formula<Y> _formula;

  public:
    /** @brief Construct a \c ScalarMultivariateFormulaFunction from as, f. */
    ScalarMultivariateFormulaFunction(SizeType as, const Formula<Y>& f) : _argument_size(as), _formula(f) {}
    /** @brief Return the scalar formula. */
    Formula<Y> formula() const { return _formula; }
    /** @brief Convert to formula< y >. */
    operator Formula<Y>() const { return _formula; }

    /** @brief Return the domain. */
    EuclideanDomain domain() const { return EuclideanDomain(this->_argument_size); }
    /** @brief Return the number of scalar arguments. */
    SizeType argument_size() const { return _argument_size; }
    /** @brief Return the number of result components. */
    SizeOne result_size() const { return SizeOne(); }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    X operator()(const Vector<X>& x) const { return Ariadne::evaluate(_formula, x); }
    /** @brief Return the derivative with respect to the selected argument. */
    friend ScalarMultivariateFormulaFunction<Y> derivative(const ScalarMultivariateFormulaFunction<Y>& f, SizeType j) {
        return ScalarMultivariateFormulaFunction<Y>(f._argument_size, Ariadne::derivative(f._formula, j));
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ScalarMultivariateFormulaFunction<Y>& f) {
        return os << f._formula;
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Representation<ScalarMultivariateFormulaFunction<Y>>& fr) {
        const ScalarMultivariateFormulaFunction<Y>& f = fr.reference();
        return os << "ScalarMultivariateFormulaFunction(" << f._argument_size << "," << f._formula << ")";
    }
};
typedef ScalarMultivariateFormulaFunction<EffectiveNumber> EffectiveScalarMultivariateFormulaFunction;

template<class Y>
using ScalarFormulaFunction = ScalarMultivariateFormulaFunction<Y>;
typedef ScalarFormulaFunction<EffectiveNumber> EffectiveScalarFormulaFunction;

//! A vector function defined by formulae
template<class Y>
class VectorMultivariateFormulaFunction {
    SizeType _argument_size;
    Vector<Formula<Y>> _formulae;

  public:
    /** @brief Construct a \c VectorMultivariateFormulaFunction from as, f. */
    VectorMultivariateFormulaFunction(SizeType as, const List<Formula<Y>>& f) : _argument_size(as), _formulae(f) {}
    /** @brief Construct a \c VectorMultivariateFormulaFunction from as, f. */
    VectorMultivariateFormulaFunction(SizeType as, const Vector<Formula<Y>>& f) : _argument_size(as), _formulae(f) {}
    /** @brief Return the component formulae. */
    Vector<Formula<Y>> formulae() const { return _formulae; }

    /** @brief Return the selected element or component. */
    ScalarMultivariateFormulaFunction<Y> operator[](SizeType i) const { return ScalarMultivariateFormulaFunction(this->_argument_size, this->_formulae[i]); }
    /** @brief Return the domain. */
    EuclideanDomain domain() const { return EuclideanDomain(this->_argument_size); }
    /** @brief Return the number of result components. */
    SizeType result_size() const { return this->_formulae.size(); }
    /** @brief Return the number of scalar arguments. */
    SizeType argument_size() const { return this->_argument_size; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    Vector<X> operator()(const Vector<X>& x) const { return Ariadne::evaluate(this->_formulae, x); }
    /** @brief Return the derivative with respect to the selected argument. */
    friend VectorMultivariateFormulaFunction<Y> derivative(const VectorMultivariateFormulaFunction<Y>& f, SizeType j) {
        return VectorMultivariateFormulaFunction<Y>(f._argument_size, Vector<Formula<Y>>(f._formulae.size(), [&](SizeType i) { return derivative(f._formulae[i], j); }));
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const VectorMultivariateFormulaFunction<Y>& f) {
        return os << f._formulae;
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Representation<VectorMultivariateFormulaFunction<Y>>& fr) {
        const VectorMultivariateFormulaFunction<Y>& f = fr.reference();
        return os << "VectorMultivariateFormulaFunction(" << f._argument_size << "," << f._formulae << ")";
    }
};
typedef VectorMultivariateFormulaFunction<EffectiveNumber> EffectiveVectorMultivariateFormulaFunction;

template<class Y>
using VectorFormulaFunction = VectorMultivariateFormulaFunction<Y>;
typedef VectorFormulaFunction<EffectiveNumber> EffectiveVectorFormulaFunction;

typedef Pair<Nat, EffectiveFormula> CoordinateFormulaPair;
typedef List<CoordinateFormulaPair> CoordinateFormulaPairs;

/** @brief Exception reporting not formula function. */
class NotFormulaFunctionException : public std::runtime_error {
  public:
    /** @brief Construct a \c NotFormulaFunctionException from str. */
    NotFormulaFunctionException(const String& str) : std::runtime_error(str) {}
};

//! \brief Returns \a true if the function \a f is syntactically constant in the indices \a is.
template<class Y>
Bool is_constant_in(const ScalarMultivariateFormulaFunction<Y>& f, const Set<Nat>& is) { return is_constant_in(f._formula, is); }
//! \brief Returns \a true if the function \a f is syntactically affine in the indices \a is.
template<class Y>
Bool is_affine_in(const ScalarMultivariateFormulaFunction<Y>& f, const Set<Nat>& is) { return is_affine_in(f._formula, is); }
//! \brief Returns \a true if the vector function \a f is syntactically affine in the indices \a is.
template<class Y>
Bool is_affine_in(const VectorMultivariateFormulaFunction<Y>& f, const Set<Nat>& is) { return is_affine_in(f._formulae, is); }
//! \brief Returns \a true if the vector function \a f is syntactically additive (possibly with multipliers) in the indices \a is.
template<class Y>
Bool is_additive_in(const VectorMultivariateFormulaFunction<Y>& f, const Set<Nat>& is) { return is_additive_in(f._formulae, is); }

/** @brief Test affine dependence on the selected arguments. */
template<class Y>
Bool is_affine_in(const VectorMultivariateFunction<Y>& f, const Set<Nat>& is) {
    auto ff = dynamic_pointer_extract<const EffectiveVectorMultivariateFormulaFunction>(f.managed_pointer());
    if (ff == nullptr)
        ARIADNE_THROW(NotFormulaFunctionException, "is_affine_in(f,is)", "Affinity checking currently available only for formula functions.");
    return is_affine_in(ff->formulae(), is);
}
/** @brief Test additive dependence on the selected arguments. */
template<class Y>
Bool is_additive_in(const VectorMultivariateFunction<Y>& f, const Set<Nat>& is) {
    auto ff = dynamic_pointer_extract<const EffectiveVectorMultivariateFormulaFunction>(f.managed_pointer());
    if (ff == nullptr)
        ARIADNE_THROW(NotFormulaFunctionException, "is_additive_in(f,is)", "Additivity checking currently available only for formula functions.");
    return is_additive_in(ff->formulae(), is);
}

/** @brief Return the component independent of the noise variables. */
inline EffectiveVectorMultivariateFunction noise_independent_component(const EffectiveVectorMultivariateFunction& function, SizeType num_inputs) {
    auto ff = dynamic_pointer_extract<const EffectiveVectorMultivariateFormulaFunction>(function.managed_pointer());
    if (ff == nullptr)
        ARIADNE_THROW(NotFormulaFunctionException, "noise_independent_component(f,num_inputs)", "Noise independent component extraction currently available only for formula functions.");

    CoordinateFormulaPairs substitutions;
    for (auto i: range(ff->result_size(), ff->result_size() + num_inputs)) {
        substitutions.append({static_cast<Nat>(i), EffectiveFormula::zero()});
    }

    return EffectiveVectorMultivariateFormulaFunction(function.argument_size(), simplify(substitute(ff->formulae(), substitutions)));
}

/** @brief Return derivatives with respect to the input variables. */
inline Vector<EffectiveVectorMultivariateFunction> input_derivatives(const EffectiveVectorMultivariateFunction& function, SizeType num_inputs) {
    Vector<EffectiveVectorMultivariateFunction> result(num_inputs);

    auto ff = dynamic_pointer_extract<const EffectiveVectorMultivariateFormulaFunction>(function.managed_pointer());
    if (ff == nullptr)
        ARIADNE_THROW(NotFormulaFunctionException, "input_derivatives(f,num_inputs)", "Input derivatives extraction currently available only for formula functions.");

    SizeType n = function.result_size();

    for (auto j: range(num_inputs)) {
        Vector<EffectiveFormula> derivative_formulae(n);
        for (auto i: range(n)) {
            derivative_formulae[i] = simplify(derivative(ff->formulae()[i], n + j));
        }
        result[j] = EffectiveVectorMultivariateFormulaFunction(function.argument_size(), derivative_formulae);
    }

    return result;
}

//------------------------ Arithmetic scalar functions  -----------------------------------//

//! A constant function f(x)=c
template<class Y, class... ARGS>
class ConstantFunction {
    using P = InformationTag<Y>;
    using SIG = Real(ARGS...);
    using D = typename SignatureTraits<SIG>::DomainType;

    D _domain;
    Y _value;

  public:
    typedef D DomainType;
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;

    // ConstantFunction(SizeType as, const Y& c) : _argument_size(as), _value(c) { }
    /** @brief Construct a \c ConstantFunction from dom, c. */
    ConstantFunction(DomainType dom, const Y& c) : _domain(dom), _value(c) {}
    /** @brief Convert to y. */
    operator Y() const { return _value; }

    /** @brief Return the domain. */
    DomainType domain() const { return _domain; }
    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const { return _domain.dimension(); }
    /** @brief Return the number of result components. */
    SizeOne result_size() const { return SizeOne(); }
    /** @brief Return the stored value. */
    const Y& value() const { return _value; }

    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline X operator()(const ElementType<D, X>& x) const {
        return _make_constant(this->_value, x);
    }
    /** @brief Return the derivative with respect to the selected argument. */
    friend ConstantFunction<Y, ARGS...> derivative(const ConstantFunction<Y, ARGS...>& f, ElementIndexType<D>) {
        return ConstantFunction<Y, ARGS...>(f._domain, Y(0));
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ConstantFunction<Y, ARGS...>& f) {
        return os << f._value;
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Representation<ConstantFunction<Y, ARGS...>>& f) {
        return os << "ConstantFunction(" << f.reference().domain() << "," << f.reference().value() << ")";
    }

  private:
    /** @brief Construct constant from the supplied arguments. */
    template<class X>
    inline static X _make_constant(Y y, const X& x) { return make_constant(y, x); }
    /** @brief Construct constant from the supplied arguments. */
    template<class X>
    inline static X _make_constant(Y y, const Vector<X>& vx) { return make_constant(y, vx.zero_element()); }
};

template<class P>
using ConstantUnivariateFunction = ConstantFunction<P, RealScalar>;
template<class P>
using ConstantMultivariateFunction = ConstantFunction<P, RealVector>;

//! A coordinate function \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=x_i\f$.
template<class P, class... ARGS>
class CoordinateFunction {
    typedef Number<P> Y;
    using SIG = Real(ARGS...);
    using D = typename SignatureTraits<SIG>::DomainType;

  public:
    typedef D DomainType;
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;
    typedef typename SignatureTraits<SIG>::ArgumentIndexType ArgumentIndexType;

  private:
    DomainType _domain;
    ArgumentIndexType _index;

  public:
    /** @brief Construct a \c CoordinateFunction from dom, i. */
    CoordinateFunction(DomainType dom, ArgumentIndexType i) : _domain(dom), _index(i) {}
    /** @brief Return the selected index. */
    ArgumentIndexType index() const { return _index; }

    /** @brief Return the domain. */
    DomainType domain() const { return _domain; }
    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const { return _domain.dimension(); }
    /** @brief Return the number of result components. */
    SizeOne result_size() const { return SizeOne(); }

    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline X operator()(const Vector<X>& x) const { return x[this->_index]; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline X operator()(const Scalar<X>& x) const { return x; }

    /** @brief Return the derivative with respect to the selected argument. */
    friend ConstantFunction<Y, ARGS...> derivative(const CoordinateFunction<P, ARGS...>& f, ArgumentIndexType j) {
        if (j == f._index) {
            return ConstantFunction<Y, ARGS...>(f.domain(), Y(1));
        } else {
            return ConstantFunction<Y, ARGS...>(f.domain(), Y(0));
        }
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const CoordinateFunction<P, ARGS...>& f) {
        if constexpr (Same<Tuple<ARGS...>, Tuple<RealScalar>>) {
            return os << "x";
        } else {
            return os << "x[" << f._index << "]";
        }
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Representation<CoordinateFunction<P, ARGS...>>& f) {
        return os << "CoordinateFunction(" << f.reference().argument_size() << "," << f.reference().index() << ")";
    }
};

//! \brief The identity function \f$ x\mapsto x\f$ in \f$\R^n\f$.
template<class P, class... ARGS>
class UnaryFunction {
    typedef Number<P> Y;
    using SIG = RealVector(ARGS...);
    using D = typename SignatureTraits<SIG>::DomainType;

  public:
    typedef D DomainType;
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;
    typedef typename SignatureTraits<SIG>::ArgumentIndexType ArgumentIndexType;

    /** @brief Construct a \c UnaryFunction from op, arg. */
    UnaryFunction(const UnaryElementaryOperator& op, const ScalarFunction<P, ARGS...>& arg)
        : _op(op), _arg(arg) {}
    /** @brief Return the domain. */
    DomainType domain() const { return this->_arg.domain(); }
    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const { return this->_arg.argument_size(); }
    /** @brief Return the number of result components. */
    SizeOne result_size() const { return SizeOne(); }

    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline X operator()(const ElementType<D, X>& x) const {
        return this->_op(this->_arg(x));
    }
    /** @brief Return the derivative with respect to the selected argument. */
    friend ScalarFunction<P, ARGS...> derivative(const UnaryFunction<P, ARGS...>& f, ArgumentIndexType j) {
        return f._op.accept([&](auto op) -> ScalarFunction<P, ARGS...> {
            if constexpr(Same<decltype(op),Abs>) { ARIADNE_FAIL_MSG("Cannot differentiate a non-smooth absolute-value function."); }
            else { return op.derivative(f._arg,f._arg.derivative(j)); } });
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const UnaryFunction<P, ARGS...>& f) {
        return os << f._op << '(' << f._arg << ')';
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Representation<UnaryFunction<P, ARGS...>>& f) {
        return os << "UnaryFunction(" << f.reference().argument_size() << "," << f.reference() << ")";
    }

  private:
    UnaryElementaryOperator _op;
    ScalarFunction<P, ARGS...> _arg;
};

/** @brief Return the square. */
template<class P, class... ARGS>
ScalarFunction<P, ARGS...> sqr(const ScalarFunction<P, ARGS...>& f) {
    return ScalarFunction<P, ARGS...>(new UnaryFunction<P, ARGS...>(OperatorCode::SQR, f));
}
/** @brief Return the square root. */
template<class P, class... ARGS>
ScalarFunction<P, ARGS...> sqrt(const ScalarFunction<P, ARGS...>& f) {
    return ScalarFunction<P, ARGS...>(new UnaryFunction<P, ARGS...>(OperatorCode::SQRT, f));
}
/** @brief Return the sine. */
template<class P, class... ARGS>
ScalarFunction<P, ARGS...> sin(const ScalarFunction<P, ARGS...>& f) {
    return ScalarFunction<P, ARGS...>(new UnaryFunction<P, ARGS...>(OperatorCode::SIN, f));
}
/** @brief Return the cosine. */
template<class P, class... ARGS>
ScalarFunction<P, ARGS...> cos(const ScalarFunction<P, ARGS...>& f) {
    return ScalarFunction<P, ARGS...>(new UnaryFunction<P, ARGS...>(OperatorCode::COS, f));
}

/** @brief Concrete function applying a binary operation to two operand functions. */
template<class P, class... ARGS>
class BinaryFunction {
    typedef Number<P> Y;
    using SIG = RealVector(ARGS...);
    using D = typename SignatureTraits<SIG>::DomainType;

  public:
    typedef D DomainType;
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;

    /** @brief Construct a \c BinaryFunction from op, arg1, arg2. */
    BinaryFunction(BinaryArithmeticOperator op, const ScalarFunction<P, ARGS...>& arg1, const ScalarFunction<P, ARGS...>& arg2)
        : BinaryFunction(BinaryElementaryOperator(op.code()), arg1, arg2) {}
    /** @brief Construct a \c BinaryFunction from op, arg1, arg2. */
    BinaryFunction(BinaryElementaryOperator op, const ScalarFunction<P, ARGS...>& arg1, const ScalarFunction<P, ARGS...>& arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) { ARIADNE_ASSERT_MSG(arg1.argument_size() == arg2.argument_size(), "op='" << op << "', arg1=" << arg1 << ", arg2=" << arg2); }
    /** @brief Return the domain. */
    DomainType domain() const {
        return intersection(this->_arg1.domain(), this->_arg2.domain());
    }
    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const {
        return this->_arg1.argument_size();
    }
    /** @brief Return the number of result components. */
    SizeOne result_size() const {
        return SizeOne();
    }

    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline X operator()(const ElementType<D, X>& x) const {
        return this->_op(this->_arg1(x), this->_arg2(x));
    }

    /** @brief Return the derivative with respect to the selected argument. */
    friend ScalarFunction<P, ARGS...> derivative(const BinaryFunction<P, ARGS...>& f, ElementIndexType<D> j) {
        return f._op.accept([&](auto op) -> ScalarFunction<P, ARGS...> {
            if constexpr (Same<decltype(op), Max> || Same<decltype(op), Min>) {
                ARIADNE_FAIL_MSG("Cannot differentiate a non-smooth maximum or minimum function.");
            } else {
                return op.derivative(f._arg1, f._arg1.derivative(j), f._arg2, f._arg2.derivative(j));
            }
        });
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const BinaryFunction<P, ARGS...>& f) {
        if (f._op.code() == OperatorCode::ADD || f._op.code() == OperatorCode::SUB) {
            return os << '(' << f._arg1 << symbol(f._op.code()) << f._arg2 << ')';
        } else {
            return os << f._arg1 << symbol(f._op.code()) << f._arg2;
        }
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Representation<BinaryFunction<P, ARGS...>>& f) {
        return os << "BinaryFunction(" << f.reference().argument_size() << "," << f.reference() << ")";
    }

  private:
    BinaryElementaryOperator _op;
    ScalarFunction<P, ARGS...> _arg1;
    ScalarFunction<P, ARGS...> _arg2;
};

// \brief The power function \f$(x,n)\mapsto x^n\f$.
/** @brief Concrete function applying an integer-graded operation to an operand function. */
template<class P, class... ARGS>
class GradedFunction {
    typedef Number<P> Y;
    using SIG = Real(ARGS...);
    using D = typename SignatureTraits<SIG>::DomainType;

  public:
    typedef D DomainType;
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;
    /** @brief Construct a \c GradedFunction from op, arg1, arg2. */
    GradedFunction(GradedElementaryOperator op, const ScalarFunction<P, ARGS...>& arg1, const Int& arg2)
        : _op(op), _arg1(arg1), _arg2(arg2) {}
    /** @brief Return the domain. */
    DomainType domain() const {
        return this->_arg1.domain();
    }
    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const {
        return this->_arg1.argument_size();
    }
    /** @brief Return the number of result components. */
    SizeOne result_size() const {
        return SizeOne();
    }

    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline X operator()(const ElementType<D, X>& x) const {
        return this->_op(this->_arg1.evaluate(x), this->_arg2);
    }
    /** @brief Return the derivative with respect to the selected argument. */
    friend Function<P, SIG> derivative(const GradedFunction<P, ARGS...>& f, ElementIndexType<D> j) {
        return f._op.accept([&](auto op) { return op.derivative(f._arg1, f._arg1.derivative(j), f._arg2); });
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const GradedFunction<P, ARGS...>& f) {
        return os << f._op.code() << "(" << f._arg1 << "," << f._arg2 << ")";
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Representation<GradedFunction<P, ARGS...>>& f) {
        return os << "GradedFunction(" << f.reference().argument_size() << "," << f.reference() << ")";
    }

  private:
    GradedElementaryOperator _op;
    ScalarFunction<P, ARGS...> _arg1;
    Int _arg2;
};

template<class P>
using UnaryMultivariateFunction = UnaryFunction<P, RealVector>;
template<class P>
using BinaryMultivariateFunction = BinaryFunction<P, RealVector>;
template<class P>
using GradedMultivariateFunction = GradedFunction<P, RealVector>;

//------------------------ Vector of Scalar functions  -----------------------------------//

/** @brief Scalar-function base with a fixed argument dimension. */
template<class P, class... ARGS>
class NonResizableScalarFunction : public CheckedAssignable<ScalarFunction<P, ARGS...>> {
  public:
    /** @brief Assign the supplied value and return this object. */
    NonResizableScalarFunction<P, ARGS...>& operator=(const ScalarFunction<P, ARGS...>& f) {
        ARIADNE_PRECONDITION_MSG(this->domain() == f.domain(), "this->domain()=" << this->domain() << ", f.domain()=" << f.domain() << "\n\n*this=" << *this << "\nf=" << f << "\n");
        this->ScalarFunction<P, ARGS...>::operator=(f);
        return *this;
    }
};

/** @brief Vector-valued function represented by scalar component functions. */
template<class P, class... ARGS>
class VectorOfScalarFunction {
    using SIG = RealVector(ARGS...);
    using D = typename SignatureTraits<SIG>::DomainType;
    friend class FunctionConstructors<P>;

  public:
    typedef D DomainType;
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;
    typedef typename SignatureTraits<SIG>::ArgumentIndexType ArgumentIndexType;

    /** @brief Construct a \c VectorOfScalarFunction from rs, as. */
    VectorOfScalarFunction(SizeType rs, SizeType as)
        : VectorOfScalarFunction(rs, ScalarFunction<P, ARGS...>(as)) {}
    /** @brief Construct a \c VectorOfScalarFunction from rs, dom. */
    VectorOfScalarFunction(SizeType rs, DomainType dom)
        : VectorOfScalarFunction(rs, ScalarFunction<P, ARGS...>(dom)) {}
    /** @brief Construct a \c VectorOfScalarFunction from rs, f. */
    VectorOfScalarFunction(SizeType rs, const ScalarFunction<P, ARGS...>& f)
        : _dom(f.domain()), _vec(rs, f) {}
    /** @brief Construct a \c VectorOfScalarFunction from vsf. */
    VectorOfScalarFunction(const Vector<ScalarFunction<P, ARGS...>>& vsf)
        : _dom(vsf.zero_element().domain()), _vec(vsf) {}

    /** @brief Replace the selected value. */
    Void set(SizeType i, ScalarFunction<P, ARGS...> f) {
        if (this->argument_size() == 0u) {
            this->_dom = f.domain();
        }
        ARIADNE_ASSERT(f.argument_size() == this->argument_size());
        this->_vec[i] = f;
    }
    /** @brief Return the requested element or stored value. */
    ScalarFunction<P, ARGS...> get(SizeType i) const {
        return this->_vec[i];
    }

    /** @brief Return the number of result components. */
    SizeType result_size() const {
        return _vec.size();
    }
    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const {
        return _dom.dimension();
    }
    /** @brief Return the domain. */
    DomainType domain() const {
        return _dom;
    }

    /** @brief Return the selected element or component. */
    const ScalarFunction<P, ARGS...> operator[](SizeType i) const {
        return this->_vec[i];
    }
    /** @brief Return the selected element or component. */
    NonResizableScalarFunction<P, ARGS...>& operator[](SizeType i) {
        return static_cast<NonResizableScalarFunction<P, ARGS...>&>(this->_vec[i]);
    }

    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline Vector<X> operator()(const ElementType<D, X>& x) const {
        Vector<X> r(this->_vec.size(), zero_element(x));
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = _vec[i].evaluate(x);
        }
        return r;
    }

    /** @brief Return the derivative with respect to the selected argument. */
    friend VectorOfScalarFunction<P, ARGS...> derivative(const VectorOfScalarFunction<P, ARGS...>& f, ArgumentIndexType k) {
        auto chrs = get_characteristics(derivative(f._vec.zero_element(), k));
        auto r = VectorOfScalarFunction<P, ARGS...>(Vector<ScalarFunction<P, ARGS...>>(f._vec.size(), [&](SizeType i) { return derivative(f._vec[i], k); }, chrs));
        return r;
    }

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const VectorOfScalarFunction<P, ARGS...>& f) {
        os << "[";
        for (SizeType i = 0; i != f._vec.size(); ++i) {
            if (i != 0) {
                os << ",";
            }
            os << f[i];
        }
        return os << "]";
    }

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Representation<VectorOfScalarFunction<P, ARGS...>>& fr) {
        VectorOfScalarFunction<P, ARGS...> f = fr.reference();
        // os << "VoSF[R" << this->argument_size() << "->R" << this->result_size() << "]";
        os << "[";
        for (SizeType i = 0; i != f.result_size(); ++i) {
            if (i != 0) {
                os << ",";
            }
            os << representation(f[i]);
        }
        return os << "]";
    }

  private:
    DomainType _dom;
    Vector<ScalarFunction<P, ARGS...>> _vec;
};

template<class P>
using VectorOfScalarUnivariateFunction = VectorOfScalarFunction<P, RealScalar>;
template<class P>
using VectorOfScalarMultivariateFunction = VectorOfScalarFunction<P, RealVector>;

/** @brief Scalar view of one component of a vector-valued function. */
template<class P, class... ARGS>
class FunctionElement {
    using SIG = RealScalar(ARGS...);
    using D = typename SignatureTraits<SIG>::DomainType;

  public:
    typedef D DomainType;
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;
    typedef typename SignatureTraits<SIG>::ArgumentIndexType ArgumentIndexType;

    /** @brief Construct a \c FunctionElement from vf, i. */
    FunctionElement(const VectorFunction<P, ARGS...>& vf, SizeType i)
        : _vf(vf), _i(i) { ARIADNE_PRECONDITION(i < vf.result_size()); }

    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const { return _vf.argument_size(); }
    /** @brief Return the number of result components. */
    SizeOne result_size() const { return SizeOne(); }
    /** @brief Return the domain. */
    DomainType domain() const { return _vf.domain(); }

    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline X operator()(const ElementType<D, X>& x) const {
        return this->_vf.evaluate(x)[_i];
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FunctionElement<P, ARGS...>& f) {
        return os << f._vf << "[" << f._i << "]";
    }

    /** @brief Return the derivative with respect to the selected argument. */
    friend FunctionElement<P, ARGS...> derivative(const FunctionElement<P, ARGS...>& f, ArgumentIndexType j) {
        return FunctionElement<P, ARGS...>(f._vf.derivative(j), f._i);
    }

  private:
    VectorFunction<P, ARGS...> _vf;
    SizeType _i;
};

//------------------------ Results of functional operations  -----------------------------------//

/** @brief Function lifted to a larger argument space. */
template<class P, class D1, class D2, class D3, class C>
class EmbeddedFunction {
    typedef CartesianProductType<D1, D2, D3> D;
    using ARG1 = ElementKind<D1>;
    using ARG2 = ElementKind<D2>;
    using ARG3 = ElementKind<D3>;
    using ARG = ElementKind<D>;
    using RES = ElementKind<C>;
    using SIG = RES(ARG);

  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;
    typedef typename SignatureTraits<SIG>::ResultSizeType ResultSizeType;

    /** @brief Construct a \c EmbeddedFunction from as1, f2, as3. */
    EmbeddedFunction(ElementSizeType<D1> as1, const Function<P, RES(ARG2)>& f2, ElementSizeType<D3> as3)
        : EmbeddedFunction(D1(as1), f2, D3(as3)) {}
    /** @brief Construct a \c EmbeddedFunction from dom1, f2, dom3. */
    EmbeddedFunction(D1 dom1, const Function<P, RES(ARG2)>& f2, D3 dom3)
        : _dom1(dom1), _f2(f2), _dom3(dom3) {}
    /** @brief Return the domain. */
    DomainType domain() const { return product(_dom1, _f2.domain(), _dom3); }
    /** @brief Return the codomain. */
    const CodomainType codomain() const { return _f2.codomain(); }
    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const { return _dom1.dimension() + _f2.argument_size() + _dom3.dimension(); }
    /** @brief Return the number of result components. */
    ResultSizeType result_size() const { return _f2.result_size(); }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline ElementType<C, X> operator()(const Vector<X>& x) const {
        Vector<X> px = project(x, Range(_dom1.dimension(), _dom1.dimension() + _f2.argument_size()));
        return _f2.evaluate(px);
    }
    /** @brief Return the selected element or component. */
    template<class I>
    inline decltype(auto) operator[](I i) const { return EmbeddedFunction<P, D1, D2, D3, RealDomain>(_dom1, _f2[i], _dom3); }
    /** @brief Return the derivative with respect to the selected argument. */
    friend EmbeddedFunction<P, D1, D2, D3, C> derivative(const EmbeddedFunction<P, D1, D2, D3, C>& f, ArgumentSizeType k) {
        ARIADNE_PRECONDITION(k < f.argument_size());
        const SizeType n1 = f._dom1.dimension();
        const SizeType n2 = f._f2.argument_size();
        if (k < n1 || k >= n1 + n2) {
            Function<P, RES(ARG2)> z(f._f2.result_size(), f._f2.domain());
            return EmbeddedFunction<P, D1, D2, D3, C>(f._dom1, z, f._dom3);
        }
        return EmbeddedFunction<P, D1, D2, D3, C>(f._dom1, f._f2.derivative(k - n1), f._dom3);
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const EmbeddedFunction<P, D1, D2, D3, C>& f) {
        return os << "EmbeddedFunction( dom1=" << f._dom1 << ", f2=" << f._f2 << ", dom3=" << f._dom3 << " )";
    }

  private:
    D1 _dom1;
    Function<P, RES(ARG2)> _f2;
    D3 _dom3;
};

/** @brief Function represented as composition of outer and inner functions. */
template<class P, class R, class T, class... AS>
class ComposedFunction {
    using SIG = R(AS...);

    using D = typename DomainTraits<AS...>::EntireDomainType;
    using E = typename DomainTraits<T>::EntireDomainType;
    using C = typename DomainTraits<R>::EntireDomainType;

  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;
    typedef typename SignatureTraits<SIG>::ResultSizeType ResultSizeType;
    typedef typename SignatureTraits<SIG>::ArgumentIndexType ArgumentIndexType;
    typedef ElementIndexType<C> ResultIndexType;

    /** @brief Construct a \c ComposedFunction from f, g. */
    ComposedFunction(const Function<P, R(T)>& f, const Function<P, T(AS...)>& g)
        : _f(f), _g(g) { ARIADNE_PRECONDITION(f.argument_size() == g.result_size()); }
    /** @brief Return the domain. */
    DomainType domain() const { return _g.domain(); }
    /** @brief Return the codomain. */
    const CodomainType codomain() const { return _f.codomain(); }
    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const { return _g.argument_size(); }
    /** @brief Return the number of result components. */
    ResultSizeType result_size() const { return _f.result_size(); }
    /** @brief Return the selected element or component. */
    ScalarFunction<P, AS...> operator[](ResultIndexType i) const { return compose(_f[i], _g); }

    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline ElementType<C, X> operator()(const ElementType<D, X>& x) const {
        return _f.evaluate(_g.evaluate(x));
    }
    /** @brief Return the derivative with respect to the selected argument. */
    Function<P, R(AS...)> derivative(ArgumentIndexType j) const {
        if constexpr (Same<R, RealScalar>) {
            if constexpr (Same<T, RealScalar>) {
                return compose(this->_f.derivative(IndexZero()), this->_g) * this->_g.derivative(j);
            } else {
                ScalarFunction<P, AS...> r = ScalarFunction<P, AS...>::zero(this->_g.domain());
                for (SizeType k = 0; k != this->_g.result_size(); ++k) {
                    r = r + compose(this->_f.derivative(k), this->_g) * this->_g[k].derivative(j);
                }
                return r;
            }
        } else {
            Function<P, R(AS...)> r(this->result_size(), this->_g.domain());
            for (SizeType i = 0; i != this->result_size(); ++i) {
                ScalarFunction<P, AS...> ri;
                if constexpr (Same<T, RealScalar>) {
                    ri = compose(this->_f[i].derivative(IndexZero()), this->_g) * this->_g.derivative(j);
                } else {
                    ri = ScalarFunction<P, AS...>::zero(this->_g.domain());
                    for (SizeType k = 0; k != this->_g.result_size(); ++k) {
                        ri = ri + compose(this->_f[i].derivative(k), this->_g) * this->_g[k].derivative(j);
                    }
                }
                r.set(i, ri);
            }
            return r;
        }
    }
    /** @brief Return the derivative with respect to the selected argument. */
    friend Function<P, R(AS...)> derivative(const ComposedFunction<P, R, T, AS...>& f, ArgumentIndexType j) {
        return f.derivative(j);
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ComposedFunction<P, R, T, AS...>& f) {
        return os << "ComposedFunction( f=" << f._f << ", g=" << f._g << " )";
    }

  private:

  private:
    Function<P, R(T)> _f;
    Function<P, T(AS...)> _g;
};

/** @brief Vector function concatenating the results of its operand functions. */
template<class P, class D, class C1, class C2>
class JoinedFunction {
    typedef CartesianProductType<C1, C2> C;

    using ARG = ElementKind<D>;
    using RES1 = ElementKind<C2>;
    using RES2 = ElementKind<C1>;
    using RES = ElementKind<C>;
    using SIG = RES(ARG);
    static_assert(Same<typename SignatureTraits<SIG>::CodomainType, CartesianProductType<C1, C2>>);

  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;
    typedef typename SignatureTraits<SIG>::ResultSizeType ResultSizeType;
    typedef typename SignatureTraits<SIG>::ArgumentIndexType ArgumentIndexType;

    /** @brief Construct a \c JoinedFunction from f1, f2. */
    JoinedFunction(Function<P, RES1(ARG)> f1, Function<P, RES2(ARG)> f2)
        : _f1(f1), _f2(f2) { ARIADNE_PRECONDITION(f1.argument_size() == f2.argument_size()); }
    /** @brief Return the selected element or component. */
    ScalarFunction<P, ARG> operator[](SizeType i) const {
        return (i < _f1.result_size()) ? _f1[i] : _f2[i - _f1.result_size()];
    }

    /** @brief Return the domain. */
    DomainType domain() const { return intersection(_f1.domain(), _f2.domain()); }
    /** @brief Return the codomain. */
    const CodomainType codomain() const { return product(_f1.codomain(), _f2.codomain()); }
    /** @brief Return the number of result components. */
    SizeType result_size() const { return _f1.result_size() + _f2.result_size(); }
    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const { return _f1.argument_size(); }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline ElementType<C, X> operator()(const ElementType<D, X>& x) const {
        return join(_f1.evaluate(x), _f2.evaluate(x));
    }
    /** @brief Return the derivative with respect to the selected argument. */
    friend JoinedFunction<P, D, C1, C2> derivative(const JoinedFunction<P, D, C1, C2>& f, ArgumentIndexType j) {
        return JoinedFunction<P, D, C1, C2>(f._f1.derivative(j), f._f2.derivative(j));
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const JoinedFunction<P, D, C1, C2>& f) {
        return os << "JoinedFunction( f1=" << f._f1 << ", f2=" << f._f2 << " )";
    }

  private:
    Function<P, RES1(ARG)> _f1;
    Function<P, RES2(ARG)> _f2;
};

/** @brief Function combining operands acting on separate argument groups. */
template<class P, class D1, class D2, class C1, class C2>
class CombinedFunction {
    typedef CartesianProductType<D1, D2> D;
    typedef CartesianProductType<C1, C2> C;
    using RES = ElementKind<C>;
    using ARG = ElementKind<C>;
    using SIG = RES(ARG);
    template<class PP, class DD, class CC>
    using FunctionType = Function<PP, ElementKind<CC>(ElementKind<DD>)>;
    template<class PP, class DD, class CC>
    using FunctionInterfaceType = typename FunctionType<PP, DD, CC>::Interface;

  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;
    typedef typename SignatureTraits<SIG>::ResultSizeType ResultSizeType;
    typedef typename SignatureTraits<SIG>::ArgumentIndexType ArgumentIndexType;

    /** @brief Construct a \c CombinedFunction from f1, f2. */
    CombinedFunction(FunctionType<P, D1, C1> f1, FunctionType<P, D2, C2> f2)
        : _f1(f1), _f2(f2) {}
    /** @brief Return the domain. */
    DomainType domain() const { return product(_f1.domain(), _f2.domain()); }
    /** @brief Return the codomain. */
    CodomainType codomain() const { return product(_f1.codomain(), _f2.codomain()); }
    /** @brief Return the number of result components. */
    ResultSizeType result_size() const { return _f1.result_size() + _f2.result_size(); }
    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const { return _f1.argument_size() + _f2.argument_size(); }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline ElementType<C, X> operator()(const ElementType<D, X>& x) const {
        return join(_f1.evaluate(project(x, range(0, _f1.argument_size()))),
                    _f2.evaluate(project(x, range(_f1.argument_size(), this->argument_size()))));
    }
    /** @brief Return the derivative with respect to the selected argument. */
    CombinedFunction<P, D1, D2, C1, C2> derivative(const CombinedFunction<P, D1, D2, C1, C2>&, ArgumentIndexType) {
        ARIADNE_NOT_IMPLEMENTED;
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const CombinedFunction<P, D1, D2, C1, C2>& f) {
        return os << "CombinedFunction( f1=" << f._f1 << ", f2=" << f._f2 << " )";
    }

  private:
    FunctionType<P, D1, C1> _f1;
    FunctionType<P, D2, C2> _f2;
};

/** @brief Function selecting result coordinates of another function. */
template<class P, class D>
class ProjectedFunction {
    using ARG = ElementKind<D>;
    using RES = RealVector;
    using SIG = RES(ARG);

  public:
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;
    typedef typename SignatureTraits<SIG>::ArgumentIndexType ArgumentIndexType;

    /** @brief Construct a \c ProjectedFunction from f, prj. */
    ProjectedFunction(VectorFunction<P, D> f, Projection prj)
        : _f(f), _prj(prj) { ARIADNE_PRECONDITION(f.result_size() == prj.argument_size()); }
    /** @brief Return the number of result components. */
    SizeType result_size() const { return _prj.result_size(); }
    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const { return _f.argument_size(); }

    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline Vector<X> operator()(const ElementType<D, X>& x) const {
        return _prj(_f(x));
    }
    /** @brief Return the derivative with respect to the selected argument. */
    friend ProjectedFunction<P, D> derivative(const ProjectedFunction<P, D>& f, ArgumentIndexType j) {
        return ProjectedFunction<P, D>(f._f.derivative(j), f._prj);
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ProjectedFunction<P, D>& f) {
        return os << "ProjectedFunction( f=" << f._f << ", prj=" << f._prj << " )";
    }

  private:
    VectorFunction<P, D> _f;
    Projection _prj;
};

// A Lie deriviative \f$\nabla g\cdot f\f$.
/** @brief Function representing differentiation along a vector field. */
template<class P>
class LieDerivativeFunction {
  public:
    //! \brief Construct the identity function in dimension \a n.
    LieDerivativeFunction(const ScalarMultivariateFunction<P>& g, const VectorMultivariateFunction<P>& f)
        : _g(g), _dg(g.argument_size(), [&g](SizeType j) { return g.derivative(j); }), _f(f) {
        ARIADNE_ASSERT(g.argument_size() == f.argument_size());
        ARIADNE_ASSERT(f.result_size() == f.argument_size());
    }

    /** @brief Return the number of result components. */
    SizeOne result_size() const { return _g.result_size(); }
    /** @brief Return the number of scalar arguments. */
    SizeType argument_size() const { return _g.argument_size(); }
    /** @brief Return the domain. */
    EuclideanDomain domain() const { return _g.domain(); }

    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    inline X operator()(const Vector<X>& x) const {
        // const Vector<R> fx=_f.evaluate(x); r=0; for(SizeType i=0; i!=_dg.size(); ++i) { r+=fx[i]+_dg[i].evaluate(x); } }
        Vector<X> fx = _f.evaluate(x);
        X r = x.zero_element();
        for (SizeType i = 0; i != _dg.size(); ++i) {
            r += fx[i] * _dg[i].evaluate(x);
        }
        return r;
    }
    /** @brief Return the derivative with respect to the selected argument. */
    friend LieDerivativeFunction<P> derivative(const LieDerivativeFunction<P>&, SizeType) {
        ARIADNE_NOT_IMPLEMENTED;
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const LieDerivativeFunction<P>& f) {
        return os << "LieDerivativeFunction( g=" << f._g << ", f=" << f._f << " )";
    }

  private:
    ScalarMultivariateFunction<P> _g;
    Vector<ScalarMultivariateFunction<P>> _dg;
    VectorMultivariateFunction<P> _f;
};

} // namespace Ariadne

#endif
