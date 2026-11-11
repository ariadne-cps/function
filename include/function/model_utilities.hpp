/***************************************************************************
 *            function/model_utilities.hpp
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

/** @file function/model_utilities.hpp
 *  @brief Utility operations for function models.
 */

#include "numeric/rounding.hpp"

namespace Ariadne {

namespace {

/** @brief Approximate centre accompanied by validated lower and upper bounds. */
template<class F>
struct ValidatedApproximation {
    Bounds<F> _v;
    Approximation<F> _a;
    /** @brief Construct a \c ValidatedApproximation from x. */
    ValidatedApproximation(const Bounds<F>& x) : _v(x), _a(x) {}
    /** @brief Convert to const bounds< f > &. */
    operator const Bounds<F>&() const { return _v; }
    /** @brief Return the lower endpoint. */
    LowerBound<F> lower() const { return _v.lower(); }
    /** @brief Return the middle value. */
    Approximation<F> middle() const { return _a; }
    /** @brief Return the upper endpoint. */
    UpperBound<F> upper() const { return _v.upper(); }
    /** @brief Return the underlying lower-endpoint representation. */
    const F& lower_raw() const { return _v.lower_raw(); }
    /** @brief Return the underlying middle-value representation. */
    const F& middle_raw() const { return _a.raw(); }
    /** @brief Return the underlying upper-endpoint representation. */
    const F& upper_raw() const { return _v.upper_raw(); }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ValidatedApproximation<F>& x) {
        return os << "{" << x._v.lower() << ":" << x._a << ":" << x._v.upper() << "}";
    }
};
/** @brief Construct validated approximation from the supplied arguments. */
template<class F>
const Approximation<F>& make_validated_approximation(const Approximation<F>& x) { return x; }
/** @brief Construct validated approximation from the supplied arguments. */
template<class F>
ValidatedApproximation<F> make_validated_approximation(const Bounds<F>& x) { return ValidatedApproximation<F>(x); }

/** @brief Convert to the rounded representation. */
template<ARawFloat F>
Rounded<F> cast_rounded(const F& x) { return Rounded<F>(x); }

/** @brief Accumulate a value while accounting for rounding error. */
template<class F>
inline Void acc_err(const F& ml, const F& u, F& e) {
    e = add(rounded, e, hlf(add(rounded, ml, u)));
}

/** @brief Return the sum of the operands and account for rounding error. */
template<class F>
F add_err(const F& x1, const F& x2, Error<F>& e) {
    F mx1 = -x1;
    RoundingSession<F> rounding;
    rounding.set(near);
    F r(add(rounded, x1, x2));
    rounding.set(up);
    F u = add(rounded, x1, x2);
    F ml = sub(rounded, mx1, x2);
    acc_err(ml, u, e.raw());
    return r;
}

/** @brief Return the sum of the operands and account for rounding error. */
template<class F>
F add_err(const F& x, const ValidatedApproximation<F>& c, Error<F>& e) {
    const F& xv = x;
    const F& cl = c.lower_raw();
    const F& cm = c.middle_raw();
    const F& cu = c.upper_raw();
    F& re = e.raw();
    RoundingSession<F> rounding;
    rounding.set(near);
    F rv = add(rounded, xv, cm);
    rounding.set(up);
    F u = add(rounded, xv, cu);
    F ml = sub(rounded, (-xv), cl);
    acc_err(ml, u, re);
    return F(rv);
}

/** @brief Return the sum of the operands and account for rounding error. */
template<class F>
F add_err(const F& x, const Bounds<F>& c, Error<F>& e) {
    return add_err(x, ValidatedApproximation<F>(c), e);
}

/** @brief Return the sum of the operands and account for rounding error. */
template<class F>
Bounds<F> add_err(const Bounds<F>& x1, const Bounds<F>& x2, Error<F>&) {
    return add(x1, x2);
}

/** @brief Return the sum of the operands and account for rounding error. */
template<class F>
UpperInterval<F> add_err(const UpperInterval<F>& x1, const UpperInterval<F>& x2, Error<F>&) {
    return add(x1, x2);
}

/** @brief Return the sum of the operands and account for rounding error. */
template<class F>
UpperInterval<F> add_err(const UpperInterval<F>& x1, Nat n2, Error<F>&) {
    return add(x1, n2);
}

/** @brief Return the sum of the operands and account for rounding error. */
template<class F>
Approximation<F> add_err(const Approximation<F>& x1, const Approximation<F>& x2, UnknownError<F>&) {
    return add(x1, x2);
}

/** @brief Return the difference of the operands and account for rounding error. */
template<class F>
F sub_err(const F& x1, const F& x2, Error<F>& e) {
    F mx1 = -x1;
    RoundingSession<F> rounding;
    rounding.set(near);
    F r(sub(rounded, x1.raw(), x2.raw()));
    rounding.set(up);
    F u = sub(rounded, x1.raw(), x2.raw());
    F ml = add(rounded, mx1.raw(), x2.raw());
    acc_err(ml, u, e.raw());
    return r;
}

/** @brief Return the difference of the operands and account for rounding error. */
template<class F>
Bounds<F> sub_err(const Bounds<F>& x1, const Bounds<F>& x2, Error<F>&) {
    return sub(x1, x2);
}

/** @brief Return the difference of the operands and account for rounding error. */
template<class F>
UpperInterval<F> sub_err(const UpperInterval<F>& x1, const UpperInterval<F>& x2, Error<F>&) {
    return sub(x1, x2);
}

/** @brief Return the difference of the operands and account for rounding error. */
template<class F>
UpperInterval<F> sub_err(const UpperInterval<F>& x1, Nat n2, Error<F>&) {
    return sub(x1, n2);
}

/** @brief Return the difference of the operands and account for rounding error. */
template<class F>
Approximation<F> sub_err(const Approximation<F>& x1, const Approximation<F>& x2, UnknownError<F>&) {
    return sub(x1, x2);
}

/** @brief Return the product of the operands and account for rounding error. */
template<class F>
F mul_err(const F& x1, const F& x2, Error<F>& e) {
    F mx1 = -x1;
    RoundingSession<F> rounding;
    rounding.set(near);
    F r(mul(rounded, x1, x2));
    rounding.set(up);
    F u = mul(rounded, x1, x2);
    F ml = mul(rounded, mx1, x2);
    acc_err(ml, u, e.raw());
    return r;
}

/** @brief Return the product of the operands and account for rounding error. */
template<class F>
F mul_err(const F& x, const ValidatedApproximation<F>& c, Error<F>& e) {
    const F& xv = x;
    const F& cu = c.upper_raw();
    const F& cm = c.middle_raw();
    const F& cl = c.lower_raw();
    F& re = e.raw();
    RoundingSession<F> rounding;
    rounding.set(near);
    F rv = mul(rounded, xv, cm);
    rounding.set(up);
    if (xv >= 0) {
        F mcl = -cl;
        F u = mul(rounded, xv, cu);
        F ml = mul(rounded, xv, mcl);
        acc_err(ml, u, re);
    } else {
        F mcu = -cu;
        F u = mul(rounded, xv, cl);
        F ml = mul(rounded, xv, mcu);
        acc_err(ml, u, re);
    }
    return F(rv);
}

/** @brief Return the product of the operands and account for rounding error. */
template<class F>
F mul_err(const F& x, const Bounds<F>& c, Error<F>& e) {
    return mul_err(x, ValidatedApproximation<F>(c), e);
}

/** @brief Return the product of the operands and account for rounding error. */
template<class F>
F mul_err(const F& x1, Nat n2, Error<F>& e) {
    return mul_err(x1, F(n2, x1.precision()), e);
}

/** @brief Return the product of the operands and account for rounding error. */
template<class F>
Bounds<F> mul_err(const Bounds<F>& x1, const Bounds<F>& x2, Error<F>&) {
    return mul(x1, x2);
}

/** @brief Return the product of the operands and account for rounding error. */
template<class F>
Bounds<F> mul_err(const Bounds<F>& x1, const ValidatedApproximation<F>& x2, Error<F>&) {
    return mul(x1, static_cast<Bounds<F>>(x2));
}

/** @brief Return the product of the operands and account for rounding error. */
template<class F>
Bounds<F> mul_err(const Bounds<F>& x1, const Nat& n2, Error<F>&) {
    return mul(x1, n2);
}

/** @brief Return the product of the operands and account for rounding error. */
template<class F>
UpperInterval<F> mul_err(const UpperInterval<F>& x1, const UpperInterval<F>& x2, Error<F>&) {
    return mul(x1, x2);
}

/** @brief Return the product of the operands and account for rounding error. */
template<class F>
UpperInterval<F> mul_err(const UpperInterval<F>& x1, Nat n2, Error<F>&) {
    return mul(x1, n2);
}

/** @brief Return the product of the operands and account for rounding error. */
template<class F>
Approximation<F> mul_err(const Approximation<F>& x1, const Approximation<F>& x2, UnknownError<F>&) {
    return mul(x1, x2);
}

/** @brief Return the product of the operands and account for rounding error. */
template<class F>
Approximation<F> mul_err(const Approximation<F>& x1, Nat n2, UnknownError<F>& e) {
    return mul_err(x1, Approximation<F>(n2, x1.precision()), e);
}

/** @brief Return the quotient of the operands and account for rounding error. */
template<class F>
F div_err(const F& x1, const F& x2, Error<F>& e) {
    F mx1 = -x1;
    RoundingSession<F> rounding;
    rounding.set(near);
    F r(div(rounded, x1.raw(), x2.raw()));
    rounding.set(up);
    F u = div(rounded, x1.raw(), x2.raw());
    F ml = div(rounded, mx1.raw(), x2.raw());
    acc_err(ml, u, e.raw());
    return r;
}

/** @brief Return the quotient of the operands and account for rounding error. */
template<class F>
F div_err(const F& x1, Nat n2, Error<F>& e) {
    return div_err(x1, F(n2, x1.precision()), e);
}

/** @brief Return the quotient of the operands and account for rounding error. */
template<class F>
Bounds<F> div_err(const Bounds<F>& x1, const Bounds<F>& x2, Error<F>&) {
    return div(x1, x2);
}

/** @brief Return the quotient of the operands and account for rounding error. */
template<class F>
Bounds<F> div_err(const Bounds<F>& x1, Nat n2, Error<F>&) {
    return div(x1, n2);
}

/** @brief Return the quotient of the operands and account for rounding error. */
template<class F>
UpperInterval<F> div_err(const UpperInterval<F>& x1, const UpperInterval<F>& x2, Error<F>&) {
    return div(x1, x2);
}

/** @brief Return the quotient of the operands and account for rounding error. */
template<class F>
UpperInterval<F> div_err(const UpperInterval<F>& x1, Nat n2, Error<F>&) {
    return div(x1, n2);
}

/** @brief Return the quotient of the operands and account for rounding error. */
template<class F>
Approximation<F> div_err(const Approximation<F>& x1, const Approximation<F>& x2, UnknownError<F>&) {
    return div(x1, x2);
}

/** @brief Return the quotient of the operands and account for rounding error. */
template<class F>
Approximation<F> div_err(const Approximation<F>& x1, Nat n2, UnknownError<F>& e) {
    return div_err(x1, Approximation<F>(n2, x1.precision()), e);
}

/** @brief Compute the product of the first two operands plus the third and account for rounding error. */
template<class F>
F fma_err(const F& x, const F& y, F z, Error<F>& e) {
    const F& xv = x.raw();
    const F& yv = y.raw();
    const F& zv = z.raw();
    F& re = e.raw();
    RoundingSession<F> rounding;
    rounding.set(near);
    F rv = fma(rounded, xv, yv, zv);
    rounding.set(up);
    F myv = -yv;
    F mzv = -zv;
    F u = fma(rounded, xv, yv, zv);
    F ml = fma(rounded, xv, myv, mzv);
    acc_err(ml, u, re);
    return F(rv);
}

/** @brief Compute the product of the first two operands plus the third and account for rounding error. */
template<class F>
F fma_err(const ValidatedApproximation<F>& c, const F& x, F y, Error<F>& e) {
    const F& xv = x.raw();
    const F& cu = c.upper_raw();
    const F& cm = c.middle_raw();
    const F& cl = c.lower_raw();
    const F& yv = y.raw();
    F& re = e.raw();
    RoundingSession<F> rounding;
    rounding.set(near);
    F rv = fma(rounded, xv, cm, yv);
    rounding.set(up);
    F u, ml;
    if (xv >= 0) {
        F mcl = -cl;
        F myv = -yv;
        u = fma(xv, cu, yv);
        ml = fma(xv, mcl, myv);
    } else {
        F mcu = -cu;
        F myv = -yv;
        u = fma(xv, cl, yv);
        ml = fma(xv, mcu, myv);
    }
    acc_err(ml, u, e);
    return F(rv);
}

/** @brief Compute the product of the first two operands plus the third and account for rounding error. */
template<class F>
F fma_err(const Bounds<F>& c, const F& x, F y, Error<F>& e) {
    return fma_err(ValidatedApproximation<F>(c), x, y, e);
}

/** @brief Compute the product of the first two operands plus the third and account for rounding error. */
template<class F>
Bounds<F> fma_err(const Bounds<F>& x, const Bounds<F>& y, Bounds<F> z, Error<F>&) {
    return fma(x, y, z);
}

/** @brief Compute the product of the first two operands plus the third and account for rounding error. */
template<class F>
UpperInterval<F> fma_err(const UpperInterval<F>& x, const UpperInterval<F>& y, UpperInterval<F> z, Error<F>&) {
    return fma(x, y, z);
}

/** @brief Compute the product of the first two operands plus the third and account for rounding error. */
template<class F>
Bounds<F> fma_err(const ValidatedApproximation<F>& x, const Bounds<F>& y, Bounds<F> z, Error<F>&) {
    return fma(static_cast<Bounds<F>>(x), y, z);
}

/** @brief Compute the product of the first two operands plus the third and account for rounding error. */
template<class F>
Approximation<F> fma_err(const Approximation<F>& x, const Approximation<F>& y, Approximation<F> z, UnknownError<F>&) {
    return fma(x, y, z);
}

// Returns an approximation to a1*b2+a2*b2, adding error to e
/** @brief Evaluate the linear combination while accounting for rounding error. */
template<class F>
F lin_err(const F& a1, const F& b1, const F& a2, const F& b2, Error<F>& e) {
    RoundingSession<F> rounding;
    rounding.set(near);
    Rounded<F> rr = cast_rounded(a1) * cast_rounded(b2) + cast_rounded(b1) * cast_rounded(a2);
    F r = rr.raw();
    rounding.set(up);
    F mb1 = -b1;
    F mb2 = -b2;
    Rounded<F> u = cast_rounded(a1) * cast_rounded(b2) + cast_rounded(b1) * cast_rounded(a2);
    Rounded<F> ml = cast_rounded(a1) * cast_rounded(mb2) + cast_rounded(mb1) * cast_rounded(a2);
    Rounded<F> re(e.raw());
    re += hlf(ml + u);
    e = Error<F>(re.raw());
    return r;
}

} // namespace

} // namespace Ariadne

