/***************************************************************************
 *            formula.cpp
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

#include "function/formula.hpp"

#include "function/formula.tpl.hpp"
#include "numeric/operators.tpl.hpp"

namespace Ariadne {

template<class Y>
inline Operator Formula<Y>::op() const {
    return this->_root->accept([](auto fn) { return static_cast<Operator>(fn._op); });
}
template<class Y>
inline OperatorCode Formula<Y>::code() const {
    return this->op().code();
}
template<class Y>
inline OperatorKind Formula<Y>::kind() const {
    return this->op().kind();
}
template<class Y>
inline const Y& Formula<Y>::val() const {
    return std::get<ConstantFormulaNode<Y>>(*this->_root)._val;
}
template<class Y>
inline const Index& Formula<Y>::ind() const {
    return std::get<IndexFormulaNode<Y>>(*this->_root)._ind;
}
template<class Y>
inline const Formula<Y>& Formula<Y>::arg() const {
    if (std::holds_alternative<GradedFormulaNode<Y>>(*this->_root)) {
        return std::get<GradedFormulaNode<Y>>(*this->_root)._arg;
    }
    if (std::holds_alternative<ScalarFormulaNode<Y>>(*this->_root)) {
        return std::get<ScalarFormulaNode<Y>>(*this->_root)._arg;
    }
    return std::get<UnaryFormulaNode<Y>>(*this->_root)._arg;
}
template<class Y>
inline const Int& Formula<Y>::num() const {
    return std::get<GradedFormulaNode<Y>>(*this->_root)._num;
}
template<class Y>
inline const Y& Formula<Y>::cnst() const {
    return std::get<ScalarFormulaNode<Y>>(*this->_root)._cnst;
}
template<class Y>
inline const Formula<Y>& Formula<Y>::arg1() const {
    return std::get<BinaryFormulaNode<Y>>(*this->_root)._arg1;
}
template<class Y>
inline const Formula<Y>& Formula<Y>::arg2() const {
    return std::get<BinaryFormulaNode<Y>>(*this->_root)._arg2;
}

template<class Y>
inline Formula<Y>::Formula() : Formula(Y()) {}
template<class Y>
inline Formula<Y>::Formula(const Y& c) : _root(new FormulaNode<Y>(ConstantFormulaNode<Y>(Cnst(), c))) {}
template<class Y>
inline Formula<Y>::Formula(const Index& i) : _root(new FormulaNode<Y>(IndexFormulaNode<Y>(Var(), i))) {}
template<class Y>
inline Formula<Y>& Formula<Y>::operator=(const Y& c) { return *this = Formula<Y>::constant(c); }

template<class Y>
inline Formula<Y> Formula<Y>::create_zero() const { return Formula<Y>::constant(static_cast<Y>(0)); }
template<class Y>
inline Formula<Y> Formula<Y>::create_constant(const Y& c) const { return Formula<Y>::constant(c); }

template<class Y>
inline Formula<Y> Formula<Y>::zero() {
    return Formula<Y>(new FormulaNode<Y>(ConstantFormulaNode<Y>(Cnst(), static_cast<Y>(0))), PointerTag());
}
template<class Y>
inline Formula<Y> Formula<Y>::constant(const Y& c) {
    return Formula<Y>(new FormulaNode<Y>(ConstantFormulaNode<Y>(Cnst(), c)), PointerTag());
}
template<class Y>
inline Formula<Y> Formula<Y>::coordinate(SizeType j) {
    return Formula<Y>(new FormulaNode<Y>(IndexFormulaNode<Y>(Var(), Index(static_cast<Nat>(j)))), PointerTag());
}
template<class Y>
inline Vector<Formula<Y>> Formula<Y>::coordinates(SizeType n) {
    return Vector<Formula<Y>>(n, [](SizeType i) { return Formula<Y>::coordinate(i); });
}
template<class Y>
inline Vector<Formula<Y>> Formula<Y>::identity(SizeType n) {
    return Formula<Y>::coordinates(n);
}
template<class Y>
inline Formula<Y> Formula<Y>::unary(const UnaryElementaryOperator& op, const Formula<Y>& a) {
    return Formula<Y>(new FormulaNode<Y>(UnaryFormulaNode<Y>(op, a)), PointerTag());
}
template<class Y>
inline Formula<Y> Formula<Y>::binary(const BinaryElementaryOperator& op, const Formula<Y>& a1, const Formula<Y>& a2) {
    return Formula<Y>(new FormulaNode<Y>(BinaryFormulaNode<Y>(op, a1, a2)), PointerTag());
}
template<class Y>
inline Formula<Y> Formula<Y>::graded(const GradedElementaryOperator& op, const Formula<Y>& a1, Int n2) {
    return Formula<Y>(new FormulaNode<Y>(GradedFormulaNode<Y>(op, a1, n2)), PointerTag());
}
template<class Y>
inline Formula<Y> Formula<Y>::scalar(const BinaryElementaryOperator& op, const Y& c1, const Formula<Y>& a2) {
    return Formula<Y>(new FormulaNode<Y>(ScalarFormulaNode<Y>(op, c1, a2)), PointerTag());
}
// template<class Y> inline Formula<Y> Formula<Y>::scalar(const BinaryElementaryOperator& op, Formula<Y> const& a1, Y const& c2) {
//     return Formula<Y>(new FormulaNode<Y>(ScalarFormulaNode<Y>(op,a1,c2)),PointerTag()); }

template<class Y, Bool UNIVARIATE = false>
class FormulaWriter {
    Void write_formula(OutputStream& os, const Formula<Y>& f) const {
        f.node_ref().accept([this, &os](const auto& node) { this->_write(os, node); });
    }

    template<class E>
    Void write_parenthesised(OutputStream& os, const E& e) const {
        os << '(';
        write_formula(os, e);
        os << ')';
    }

    template<class E>
    Void write_unary(OutputStream& os, Pos, const E& e) const {
        os << '+';
        write_formula(os, e);
    }
    template<class E>
    Void write_unary(OutputStream& os, Neg, const E& e) const {
        os << '-';
        switch (e.op().code()) {
            case Cnst::code():
            case Var::code():
                write_formula(os, e);
                break;
            default:
                write_parenthesised(os, e);
        }
    }
    template<class OP, class E>
    Void write_unary(OutputStream& os, OP op, const E& e) const {
        os << op << '(';
        write_formula(os, e);
        os << ')';
    }

    template<class E>
    Void write_binary(OutputStream& os, Add, const E& a, const E& b) const {
        write_formula(os, a);
        os << '+';
        write_formula(os, b);
    }
    template<class E>
    Void write_binary(OutputStream& os, Sub, const E& a, const E& b) const {
        write_formula(os, a);
        os << '-';
        if (b.op().code() == Add::code() || b.op().code() == Sub::code())
            write_parenthesised(os, b);
        else
            write_formula(os, b);
    }
    template<class E>
    Void write_binary(OutputStream& os, Mul, const E& a, const E& b) const {
        if (a.op().code() == Add::code() || a.op().code() == Sub::code() || a.op().code() == Div::code())
            write_parenthesised(os, a);
        else
            write_formula(os, a);
        os << '*';
        if (b.op().code() == Add::code() || b.op().code() == Sub::code())
            write_parenthesised(os, b);
        else
            write_formula(os, b);
    }
    template<class E>
    Void write_binary(OutputStream& os, Div, const E& a, const E& b) const {
        if (a.op().code() == Add::code() || a.op().code() == Sub::code())
            write_parenthesised(os, a);
        else
            write_formula(os, a);
        os << '/';
        if (b.op().code() == Add::code() || b.op().code() == Sub::code() || b.op().code() == Mul::code() || b.op().code() == Div::code())
            write_parenthesised(os, b);
        else
            write_formula(os, b);
    }
    template<class E>
    Void write_binary(OutputStream& os, Max, const E& a, const E& b) const {
        os << "max(";
        write_formula(os, a);
        os << ',';
        write_formula(os, b);
        os << ')';
    }
    template<class E>
    Void write_binary(OutputStream& os, Min, const E& a, const E& b) const {
        os << "min(";
        write_formula(os, a);
        os << ',';
        write_formula(os, b);
        os << ')';
    }
    template<class OP, class E>
    Void write_binary(OutputStream& os, OP op, const E& a, const E& b) const {
        os << op << '(';
        write_formula(os, a);
        os << ',';
        write_formula(os, b);
        os << ')';
    }

  public:
    Void _write(OutputStream& os, const ConstantFormulaNode<Y>& node) const { os << node._val; }
    Void _write(OutputStream& os, const IndexFormulaNode<Y>& node) const {
        if constexpr (UNIVARIATE)
            os << 'x';
        else
            os << "x[" << node._ind << ']';
    }
    Void _write(OutputStream& os, const UnaryFormulaNode<Y>& node) const {
        node._op.accept([&](auto op) { write_unary(os, op, node._arg); });
    }
    Void _write(OutputStream& os, const BinaryFormulaNode<Y>& node) const {
        node._op.accept([&](auto op) { write_binary(os, op, node._arg1, node._arg2); });
    }
    Void _write(OutputStream& os, const GradedFormulaNode<Y>& node) const {
        os << node._op << '(';
        write_formula(os, node._arg);
        os << ',' << node._num << ')';
    }
    Void _write(OutputStream& os, const ScalarFormulaNode<Y>& node) const {
        os << '(';
        Formula<Y> cnst = Formula<Y>::constant(node._cnst);
        node._op.accept([&](auto op) { write_binary(os, op, cnst, node._arg); });
        os << ')';
    }

    Void write(OutputStream& os, const Formula<Y>& f) const { write_formula(os, f); }
};

template<class Y>
Formula<Y> Formula<Y>::_derivative(SizeType j) const {
    return this->_root->accept([j](const auto& fn) { return FormulaNodeDetail::derivative_node(fn, j); });
}

template<class Y>
OutputStream& Formula<Y>::_write_multivariate(OutputStream& os, const Formula<Y>& f) {
    FormulaWriter<Y> w;
    w.write(os, f);
    return os;
}

template<class Y>
OutputStream& Formula<Y>::_write_multivariate(OutputStream& os, const Vector<Formula<Y>>& f) {
    os << "[";
    for (SizeType i = 0; i != f.size(); ++i) {
        if (i != 0) {
            os << ",";
        }
        _write_multivariate(os, f[i]);
    }
    os << "]";
    return os;
}

template<class Y>
OutputStream& Formula<Y>::_write_univariate(OutputStream& os, const Formula<Y>& f) {
    FormulaWriter<Y, true> w;
    w.write(os, f);
    return os;
}

template<class Y>
OutputStream& Formula<Y>::_write_univariate(OutputStream& os, const Vector<Formula<Y>>& f) {
    os << "[";
    for (SizeType i = 0; i != f.size(); ++i) {
        if (i != 0) {
            os << ",";
        }
        _write_univariate(os, f[i]);
    }
    os << "]";
    return os;
}

template class Formula<ApproximateNumber>;
template class Formula<ValidatedNumber>;
template class Formula<EffectiveNumber>;
template class Formula<ExactNumber>;

template class Formula<Real>;

template Formula<EffectiveNumber> simplify(const Formula<EffectiveNumber>&);
template Vector<Formula<EffectiveNumber>> simplify(const Vector<Formula<EffectiveNumber>>&);
template Formula<EffectiveNumber> substitute(const Formula<EffectiveNumber>& a, const Nat& i, const Formula<EffectiveNumber>& is);

template Bool identical(const Formula<EffectiveNumber>&, const Formula<EffectiveNumber>&);

template Bool identical(const Vector<Formula<EffectiveNumber>>&, const Vector<Formula<EffectiveNumber>>&);
template Bool is_affine_in(const Vector<Formula<EffectiveNumber>>&, const Set<Nat>&);
template Bool is_additive_in(const Vector<Formula<EffectiveNumber>>&, const Set<Nat>&);

} // namespace Ariadne
