/***************************************************************************
 *            function/scaling.hpp
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

/*! @file function/scaling.hpp
 *  \brief Scaling functions.
 */

#ifndef ARIADNE_SCALING_HPP
#define ARIADNE_SCALING_HPP

#include "function/domain.hpp"
#include "numeric/dyadic.hpp"

namespace Ariadne {

/** @brief Exception reporting unscaling. */
struct UnscalingException : std::runtime_error {
    IntervalDomainType domain;
    /** @brief Construct a \c UnscalingException from msg, dom. */
    inline UnscalingException(String msg, IntervalDomainType dom) : std::runtime_error(msg + to_str(dom)), domain(dom) {}
};

/** @brief Return an approximate interval midpoint. */
template<class F>
inline Approximation<F> med_apprx(const Interval<F>& ivl) {
    return Approximation<F>(hlf(add(approx, ivl.lower_bound().raw(), ivl.upper_bound().raw())));
}

/** @brief Return an approximate interval radius. */
template<class F>
inline Approximation<F> rad_apprx(const Interval<F>& ivl) {
    return Approximation<F>(hlf(sub(approx, ivl.upper_bound().raw(), ivl.lower_bound().raw())));
}

/** @brief Return a validated interval midpoint. */
template<class F>
inline Bounds<F> med_val(const Interval<F>& ivl) {
    return hlf(ivl.lower_bound() + ivl.upper_bound());
}

/** @brief Return a validated interval radius. */
template<class F>
inline Bounds<F> rad_val(const Interval<F>& ivl) {
    return hlf(ivl.upper_bound() - ivl.lower_bound());
}

/** @brief Return the midpoint value. */
inline Dyadic med(const IntervalDomainType& ivl) {
    return hlf(add(Dyadic(ivl.lower_bound()), Dyadic(ivl.upper_bound())));
}

/** @brief Return the radius. */
inline Dyadic rad(const IntervalDomainType& ivl) {
    return hlf(sub(Dyadic(ivl.upper_bound()), Dyadic(ivl.lower_bound())));
}

/** @brief Convert the interval to the requested endpoint representation. */
inline Interval<FloatDP> convert_interval(const Interval<FloatDP>& ivl, DoublePrecision) {
    return ivl;
}

/** @brief Convert the interval to the requested endpoint representation. */
inline Interval<FloatMP> convert_interval(const Interval<FloatDP>& ivl, MultiplePrecision pr) {
    return Interval<FloatMP>(FloatMP(ivl.lower_bound().raw(), pr), FloatMP(ivl.upper_bound().raw(), pr));
}

/** @brief Map from the normalized domain to the physical domain. */
template<class T>
inline T scale(T x, const IntervalDomainType& cd) {
    return std::move(x) * rad(cd) + med(cd);
}

/** @brief Map from the physical domain to the normalized domain. */
template<class T>
inline T unscale(T x, const IntervalDomainType& d) {
    if (d.lower_bound() == d.upper_bound()) {
        return std::move(x) * 0;
        // throw UnscalingException{"Cannot unscale "+to_str(x)+" over empty interval ",d};
    } else {
        return (std::move(x) - med(d)) / rad(d);
    }
}

/** @brief Map from the normalized domain to the physical domain. */
template<class X>
inline Vector<X> scale(const Vector<X>& v, const BoxDomainType& cd) {
    ARIADNE_PRECONDITION(v.size() == cd.dimension());
    return Vector<X>(v.size(), [&](SizeType i) { return scale(v[i], cd[i]); }, v.element_characteristics());
}

/** @brief Conversion probe used to constrain scaling overloads. */
template<class F, class T>
struct MyConvertible {
    /** @brief Convert the supplied value to the requested representation. */
    template<class TT>
    static Void convert(TT const&);
    /** @brief Evaluate the compile-time or runtime predicate for the supplied arguments. */
    template<class FF, class TT, class = decltype(convert<TT>(std::declval<FF>()))>
    static std::true_type test(int);
    /** @brief Evaluate the compile-time or runtime predicate for the supplied arguments. */
    template<class FF, class TT>
    static std::false_type test(...);
    static const bool value = decltype(test<F, T>(1))::value;
};

/** @brief Map from the physical domain to the normalized domain. */
template<class X>
inline Vector<X> unscale(const Vector<X>& v, const BoxDomainType& d) {
    ARIADNE_PRECONDITION(v.size() == d.dimension());
    return Vector<X>(v.size(), [&](SizeType i) { return unscale(v[i], d[i]); }, v.element_characteristics());
}

/** @brief Affine map from a normalized interval to a physical interval. */
class Scaling {
    IntervalDomainType _codom;

  public:
    /** @brief Construct a \c Scaling from codom. */
    Scaling(IntervalDomainType codom) : _codom(codom) {}
    /** @brief Return the domain. */
    UnitInterval domain() const { return UnitInterval(); }
    /** @brief Return the codomain. */
    IntervalDomainType codomain() const { return _codom; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    X operator()(X x) const { return scale(std::move(x), _codom); }
};

/** @brief Coordinatewise affine map from a normalized box to a physical box. */
class VectorScaling {
    Box<IntervalDomainType> _codom;

  public:
    /** @brief Construct a \c VectorScaling from codom. */
    VectorScaling(Box<IntervalDomainType> codom) : _codom(codom) {}
    /** @brief Return the number of elements. */
    SizeType size() const { return _codom.dimension(); }
    /** @brief Return the selected element or component. */
    Scaling operator[](SizeType i) const { return Scaling(_codom[i]); }
    /** @brief Return the codomain. */
    const Box<IntervalDomainType>& codomain() const { return _codom; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    Vector<X> operator()(const Vector<X>& v) const { return scale(v, _codom); }
};

/** @brief Affine map from a physical interval to a normalized interval. */
class Unscaling {
    IntervalDomainType _dom;

  public:
    /** @brief Construct a \c Unscaling from dom. */
    Unscaling(IntervalDomainType dom) : _dom(dom) {}
    /** @brief Return the domain. */
    IntervalDomainType domain() const { return _dom; }
    /** @brief Return the codomain. */
    UnitInterval codomain() const { return UnitInterval(); }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    X operator()(X x) const { return unscale(std::move(x), _dom); }
};

/** @brief Coordinatewise affine map from a physical box to a normalized box. */
class VectorUnscaling {
    Box<IntervalDomainType> _dom;

  public:
    /** @brief Construct a \c VectorUnscaling from dom. */
    VectorUnscaling(Box<IntervalDomainType> dom) : _dom(dom) {}
    /** @brief Return the number of elements. */
    SizeType size() const { return _dom.dimension(); }
    /** @brief Return the selected element or component. */
    Unscaling operator[](SizeType i) const { return Unscaling(_dom[i]); }
    /** @brief Return the domain. */
    const Box<IntervalDomainType>& domain() const { return _dom; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    Vector<X> operator()(const Vector<X>& v) const { return unscale(v, _dom); }
};

} // namespace Ariadne

#endif // ARIADNE_SCALING_HPP
