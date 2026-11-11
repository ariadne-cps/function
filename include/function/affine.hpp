/***************************************************************************
 *            function/affine.hpp
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

/*! @file function/affine.hpp
 *  \brief Affine scalar and vector functions
 */

#ifndef ARIADNE_AFFINE_HPP
#define ARIADNE_AFFINE_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "algebra/covector.hpp"
#include "algebra/matrix.hpp"
#include "function/declarations.hpp"
#include "utility/macros.hpp"
#include "utility/pointer.hpp"

namespace Ariadne {

template<class PR>
class FloatFactory;

template<class X>
class Affine;

/** @brief Provide algebra operations for the supported operand types. */
template<class A, class X = typename A::NumericType>
struct ProvideAlgebraOperations;

/** @brief Provide algebra operations for the supported operand types. */
template<class X>
struct ProvideAlgebraOperations<Affine<X>, X> {
    //! \relates Affine
    //! \brief Test equality of two affine expressions.
    friend inline auto operator==(const Affine<X>& f1, const Affine<X>& f2) {
        return f1._c == f2._c && f1._g == f2._g;
    }
    //! \relates Affine
    //! \brief Negation of an affine expression.
    friend inline Affine<X> operator-(const Affine<X>& f) {
        return Affine<X>(Covector<X>(-f._g), -f._c);
    }
    //! \relates Affine
    //! \brief Addition of two affine expressions.
    friend inline Affine<X> operator+(const Affine<X>& f1, const Affine<X>& f2) {
        return Affine<X>(Covector<X>(f1._g + f2._g), f1._c + f2._c);
    }
    //! \relates Affine
    //! \brief Subtraction of two affine expressions.
    friend inline Affine<X> operator-(const Affine<X>& f1, const Affine<X>& f2) {
        return Affine<X>(Covector<X>(f1._g - f2._g), f1._c - f2._c);
    }
    //! \relates Affine
    //! \brief Addition of a constant to an affine expression.
    friend inline Affine<X> operator+(const Affine<X>& f1, const X& c2) {
        return Affine<X>(Covector<X>(f1._g), f1._c + c2);
    }
    //! \relates Affine
    //! \brief Addition of a constant to an affine expression.
    friend inline Affine<X> operator+(const X& c1, const Affine<X>& f2) {
        return Affine<X>(Covector<X>(f2._g), c1 + f2._c);
    }
    //! \relates Affine
    //! \brief Subtraction of a constant to an affine expression.
    friend inline Affine<X> operator-(const Affine<X>& f1, const X& c2) {
        return Affine<X>(Covector<X>(f1._g), f1._c - c2);
    }
    //! \relates Affine
    //! \brief Subtraction of an affine expression from a constant.
    friend inline Affine<X> operator-(const X& c1, const Affine<X>& f2) {
        return Affine<X>(Covector<X>(-f2._g), c1 - f2._c);
    }
    //! \relates Affine
    //! \brief Scalar multiplication of an affine expression.
    friend inline Affine<X> operator*(const X& c, const Affine<X>& f) {
        return Affine<X>(Covector<X>(c * f._g), c * f._c);
    }
    //! \relates Affine
    //! \brief Scalar multiplication of an affine expression.
    friend inline Affine<X> operator*(const Affine<X>& f, const X& c) { return c * f; }
    //! \relates Affine
    //! \brief Scalar division of an affine expression.
    friend inline Affine<X> operator/(const Affine<X>& f, const X& c) { return (1 / c) * f; }
    //! \relates Affine
    //! \brief The derivative of an affine expression gives a constant.
    friend inline X derivative(const Affine<X>& f, Nat k) { return f.derivative(k); }

    /** @brief Return the sum or unary positive value. */
    template<GenericNumber Y>
    friend decltype(auto) operator+(const Affine<X>& x, const Y& y) { return x + factory(x).create(y); }
    /** @brief Return the difference or additive inverse. */
    template<GenericNumber Y>
    friend decltype(auto) operator-(const Affine<X>& x, const Y& y) { return x - factory(x).create(y); }
    /** @brief Return the product. */
    template<GenericNumber Y>
    friend decltype(auto) operator*(const Affine<X>& x, const Y& y) { return x * factory(x).create(y); }
    /** @brief Return the quotient. */
    template<GenericNumber Y>
    friend decltype(auto) operator/(const Affine<X>& x, const Y& y) { return x / factory(x).create(y); }
    /** @brief Return the sum or unary positive value. */
    template<GenericNumber Y>
    friend decltype(auto) operator+(const Y& y, const Affine<X>& x) { return factory(x).create(y) + x; }
    /** @brief Return the difference or additive inverse. */
    template<GenericNumber Y>
    friend decltype(auto) operator-(const Y& y, const Affine<X>& x) { return factory(x).create(y) - x; }
    /** @brief Return the product. */
    template<GenericNumber Y>
    friend decltype(auto) operator*(const Y& y, const Affine<X>& x) { return factory(x).create(y) * x; }
};

/** @brief Return a factory for compatible objects. */
template<class X>
FloatFactory<PrecisionType<X>> factory(const Affine<X>& a) {
    return FloatFactory<PrecisionType<X>>(a.value().precision());
}

/*
template<ConcreteNumber X, GenericNumber Y>
decltype(auto) operator+(Affine<X> const& x, Y const& y) { return x+factory(x).create(y); }
template<ConcreteNumber X, GenericNumber Y>
decltype(auto) operator-(Affine<X> const& x, Y const& y) { return x-factory(x).create(y); }
template<ConcreteNumber X, GenericNumber Y>
decltype(auto) operator*(Affine<X> const& x, Y const& y) { return x*factory(x).create(y); }
template<ConcreteNumber X, GenericNumber Y>
decltype(auto) operator/(Affine<X> const& x, Y const& y) { return x/factory(x).create(y); }
template<ConcreteNumber X, GenericNumber Y>
decltype(auto) operator+(Y const& y, Affine<X> const& x) { return factory(x).create(y)+y; }
template<ConcreteNumber X, GenericNumber Y>
decltype(auto) operator*(Y const& y, Affine<X> const& x) { return factory(x).create(y)*y; }
*/

//! \ingroup FunctionModule
//! \brief An affine function \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b\f$.
template<class X>
class Affine
    : public ProvideAlgebraOperations<Affine<X>, X> {
  public:
    typedef X NumericType;
    typedef Ariadne::CharacteristicsType<X> CoefficientCharacteristicsType;
    typedef Tuple<SizeType, CoefficientCharacteristicsType> CharacteristicsType;

  public:
    /** @brief Construct a \c Affine from n, prs. */
    template<class... PRS>
        requires Constructible<X, PRS...>
    explicit Affine(SizeType n, PRS... prs) : _c(0, prs...), _g(n, prs...) {}
    /** @brief Construct a \c Affine from g, c. */
    explicit Affine(const Covector<X>& g, const X& c) : _c(c), _g(g) {}
    /** @brief Construct a \c Affine from c, g. */
    explicit Affine(X c, InitializerList<X> g) : _c(c), _g(g) {}
    /** @brief Copy-construct the object from the supplied instance. */
    template<class XX>
    explicit Affine(const Affine<XX>& aff)
        : _c(aff.b()), _g(aff.a()) {}

    /** @brief Assign the supplied value and return this object. */
    Affine<X>& operator=(const X& c) {
        this->_c = c;
        for (SizeType i = 0; i != this->_g.size(); ++i) {
            this->_g[i] = 0;
        }
        return *this;
    }
    /** @brief Construct a constant with the supplied value. */
    static Affine<X> constant(SizeType n, X c) {
        return Affine<X>(Covector<X>(n, nul(c)), c);
    }
    /** @brief Construct or return the selected coordinate. */
    template<class... PRS>
        requires Constructible<X, Nat, PRS...>
    static Affine<X> coordinate(SizeType n, SizeType j, PRS... prs) {
        return Affine<X>(Covector<X>::unit(n, j, prs...), X(0u, prs...));
    }
    /** @brief Construct or return the coordinate family. */
    template<class... PRS>
        requires Constructible<X, Nat, PRS...>
    static Vector<Affine<X>> coordinates(SizeType n, PRS... prs) {
        return Vector<Affine<X>>(n, [&](SizeType i) { return Affine<X>::coordinate(n, i, prs...); }, Affine<X>(n, prs...));
    }
    /** @brief Construct the selected independent variable. */
    template<class... PRS>
        requires Constructible<X, Nat, PRS...>
    static Affine<X> variable(SizeType n, SizeType j, PRS... prs) { return coordinate(n, j, prs...); }
    /** @brief Construct the independent variables. */
    template<class... PRS>
        requires Constructible<X, Nat, PRS...>
    static Vector<Affine<X>> variables(SizeType n, PRS... prs) { return coordinates(n, prs...); }

    /** @brief Return the selected element or component. */
    const X& operator[](SizeType i) const { return this->_g[i]; }
    /** @brief Return the selected element or component. */
    X& operator[](Nat i) { return this->_g[i]; }

    /** @brief Return the linear coefficients of the affine object. */
    const Covector<X>& a() const { return this->_g; }
    /** @brief Return the constant part of the affine object. */
    const X& b() const { return this->_c; }

    /** @brief Return the vector of first partial derivatives. */
    const Covector<X>& gradient() const { return this->_g; }
    /** @brief Return the vector of first partial derivatives. */
    const X& gradient(SizeType i) const { return this->_g[i]; }
    /** @brief Return the stored value. */
    const X& value() const { return this->_c; }

    /** @brief Return the number of scalar arguments. */
    SizeType argument_size() const { return this->_g.size(); }
    /** @brief Return the common characteristics of the coefficients. */
    CoefficientCharacteristicsType coefficient_characteristics() const {
        return get_characteristics(this->_c);
    }
    /** @brief Return the characteristics required to construct compatible values. */
    CharacteristicsType characteristics() const {
        return std::make_tuple(this->argument_size(), this->coefficient_characteristics());
    }

    /** @brief Evaluate at the supplied argument. */
    template<class Y>
    Y evaluate(const Vector<Y>& x) const;

    /** @brief Return the derivative with respect to the selected argument. */
    const X& derivative(SizeType j) const { return this->_g[j]; }

  private:

  public:
    X _c;
    Covector<X> _g;
};

template<class X>
template<class Y>
Y Affine<X>::evaluate(const Vector<Y>& x) const {
    Y r = x.zero_element();
    r += this->_c;
    for (SizeType j = 0; j != this->_g.size(); ++j) {
        r += this->_g[j] * x[j];
    }
    return r;
}

/** @brief Write the value to the output stream. */
template<class X>
OutputStream& operator<<(OutputStream& os, const Affine<X>& f) {
    os << f.b();
    for (SizeType j = 0; j != f.argument_size(); ++j) {
        os << "+" << "(" << f.a()[j] << ")*x" << j;
    }
    return os;
}

} // namespace Ariadne

#endif /* ARIADNE_AFFINE_HPP */
