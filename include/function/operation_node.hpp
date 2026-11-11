/***************************************************************************
 *            function/operation_node.hpp
 *
 *  Copyright  2013-26  Luca Geretti
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

/** @file function/operation_node.hpp
 *  @brief Expression nodes representing function operations.
 */

#ifndef ARIADNE_FUNCTION_OPERATION_NODE_HPP
#define ARIADNE_FUNCTION_OPERATION_NODE_HPP

#include "numeric/operators.hpp"

namespace Ariadne {

/** @brief Expression node storing an operation and its operands. */
template<class OP, class... AS>
struct OperationNode;

/** @brief Expression node storing an operation and its operands. */
template<class C>
struct OperationNode<Cnst, C> {
    Cnst _op;
    C _val;
    /** @brief Construct a \c OperationNode from c. */
    explicit OperationNode(C c) : _op(), _val(c) {}
    /** @brief Construct a \c OperationNode from o, c. */
    OperationNode(Cnst o, C c) : _op(o), _val(c) {}
    /** @brief Return the stored value. */
    C val() const { return _val; }
    /** @brief Convert to t. */
    template<class T>
    explicit operator T() const { return static_cast<T>(_val); }
    /** @brief Apply the stored operation to the supplied arguments. */
    template<class... AS>
    auto operator()(AS...) const -> C { return _val; }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const OperationNode& node) { return os << node._val; }
};

/** @brief Expression node storing an operation and its operands. */
template<class I>
struct OperationNode<Var, I> {
    Var _op;
    I _ind;
    /** @brief Construct a \c OperationNode from i. */
    explicit OperationNode(I i) : _op(), _ind(i) {}
    /** @brief Construct a \c OperationNode from o, i. */
    OperationNode(Var o, I i) : _op(o), _ind(i) {}
    /** @brief Return the coordinate index. */
    I ind() const { return _ind; }
    /** @brief Apply the stored operation to the supplied arguments. */
    template<class AS>
    auto operator()(AS vals) const -> decltype(vals[_ind]) { return vals[_ind]; }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const OperationNode& node) { return os << node._ind; }
};

/** @brief Expression node storing an operation and its operands. */
template<class O, class A>
struct OperationNode<O, A> {
    O _op;
    A _arg;
    /** @brief Return the stored operation. */
    O op() const { return _op; }
    /** @brief Return the operand formula. */
    A arg() const { return _arg; }
    /** @brief Construct a \c OperationNode from o, a. */
    OperationNode(O o, A a) : _op(o), _arg(a) {}
    /** @brief Apply the stored operation to the supplied arguments. */
    template<class... AS>
    auto operator()(AS... vals) const -> decltype(_op(_arg(vals...))) {
        return _op(_arg(vals...));
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const OperationNode& node) {
        return os << node._op.code() << "(" << node._arg << ")";
    }
};

/** @brief Expression node storing an operation and its operands. */
template<class O, class A1, class A2>
struct OperationNode<O, A1, A2> {
    O _op;
    A1 _arg1;
    A2 _arg2;
    /** @brief Return the stored operation. */
    O op() const { return _op; }
    /** @brief Return the first operand. */
    A1 arg1() const { return _arg1; }
    /** @brief Return the second operand. */
    A2 arg2() const { return _arg2; }
    /** @brief Construct a \c OperationNode from o, a1, a2. */
    OperationNode(O o, A1 a1, A2 a2) : _op(o), _arg1(a1), _arg2(a2) {}
    /** @brief Apply the stored operation to the supplied arguments. */
    template<class... AS>
    auto operator()(AS... vals) const -> decltype(_op(_arg1(vals...), _arg2(vals...))) {
        return _op(_arg1(vals...), _arg2(vals...));
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const OperationNode& node) {
        return os << node._op.code() << "(" << node._arg1 << "," << node._arg2 << ")";
    }
};

template<class O>
concept AGraded = Same<O, Pow> or Same<O, GradedElementaryOperator>;

/** @brief Expression node storing an operation and its operands. */
template<class O, class A, class N>
    requires AGraded<O>
struct OperationNode<O, A, N> {
    O _op;
    A _arg;
    N _num;
    /** @brief Return the stored operation. */
    O op() const { return _op; }
    /** @brief Return the operand formula. */
    A arg() const { return _arg; }
    /** @brief Return the stored numeric operand. */
    N num() const { return _num; }
    /** @brief Return the first operand. */
    A arg1() const { return _arg; }
    /** @brief Return the second operand. */
    N arg2() const { return _num; }
    /** @brief Construct a \c OperationNode from o, a, n. */
    OperationNode(O o, A a, N n) : _op(o), _arg(a), _num(n) {}
    /** @brief Apply the stored operation to the supplied arguments. */
    template<class... AS>
    auto operator()(AS... vals) const -> decltype(_op(_arg(vals...), _num)) {
        return _op(_arg(vals...), _num);
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const OperationNode& node) {
        return os << node._op.code() << "(" << node._arg << "," << node._num << ")";
    }
};

/** @brief Expression node storing an operation and its operands. */
template<class O, class A, template<class> class E>
struct OperationNode<O, A, E<A>> {
    O _op;
    A _cnst;
    E<A> _arg;
    /** @brief Construct a \c OperationNode from o, c, a. */
    OperationNode(O o, A c, E<A> a) : _op(o), _cnst(c), _arg(a) {}
    /** @brief Return the stored operation. */
    O op() const { return _op; }
    /** @brief Return the stored constant. */
    A cnst() const { return _cnst; }
    /** @brief Return the operand formula. */
    E<A> arg() const { return _arg; }
    /** @brief Convert to t. */
    template<class T>
    operator T() const { return _op(static_cast<T>(_cnst), static_cast<T>(_arg)); }
    /** @brief Apply the stored operation to the supplied arguments. */
    template<class... AS>
    auto operator()(AS... vals) const -> decltype(_op(_cnst, _arg(vals...))) {
        return _op(_cnst, _arg(vals...));
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const OperationNode& node) {
        return os << node._op.code() << "(" << node._cnst << "," << node._arg << ")";
    }
};

/** @brief Expression node storing an operation and its operands. */
template<class O, class A1, class A2, class A3>
struct OperationNode<O, A1, A2, A3> {
    O _op;
    A1 _arg1;
    A2 _arg2;
    A3 _arg3;
    /** @brief Construct a \c OperationNode from o, a1, a2, a3. */
    OperationNode(O o, A1 a1, A2 a2, A3 a3) : _op(o), _arg1(a1), _arg2(a2), _arg3(a3) {}
    /** @brief Convert to t. */
    template<class T>
    explicit operator T() const {
        return _op(static_cast<T>(_arg1), static_cast<T>(_arg2), static_cast<T>(_arg3));
    }
    /** @brief Apply the stored operation to the supplied arguments. */
    template<class... AS>
    auto operator()(AS... vals) const -> decltype(_op(_arg1(vals...), _arg2(vals...), _arg3(vals...))) {
        return _op(_arg1(vals...), _arg2(vals...), _arg3(vals...));
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const OperationNode& node) {
        return os << node._op.code() << "(" << node._arg1 << "," << node._arg2 << "," << node._arg3 << ")";
    }
};

} // namespace Ariadne

#endif
