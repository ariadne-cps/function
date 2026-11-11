/***************************************************************************
 *            function/measurable_function.hpp
 *
 *  Copyright  2020  Pieter Collins
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

/*! @file function/measurable_function.hpp
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_MEASURABLE_FUNCTION_HPP
#define ARIADNE_MEASURABLE_FUNCTION_HPP

#include "domain.hpp"
#include "function_model.hpp"
#include "geometry/measurable_set.hpp"
#include "geometry/set.hpp"
#include "geometry/set_wrapper.hpp"
#include "numeric/number.hpp"

namespace Ariadne {

/** @brief Polymorphic function exposing measurable-set preimages. */
template<class P, class SIG>
class MeasurableFunction;
template<class SIG>
using EffectiveMeasurableFunction = MeasurableFunction<EffectiveTag, SIG>;
template<class SIG>
using ValidatedMeasurableFunction = MeasurableFunction<ValidatedTag, SIG>;
using EffectiveScalarUnivariateMeasurableFunction = MeasurableFunction<EffectiveTag, Real(Real)>;
using ValidatedScalarUnivariateMeasurableFunction = MeasurableFunction<ValidatedTag, Real(Real)>;

template<class P, class SIG>
using ContinuousFunction = Function<P, SIG>;
template<class SIG>
using EffectiveContinuousFunction = ContinuousFunction<EffectiveTag, SIG>;
template<class SIG>
using ValidatedContinuousFunction = ContinuousFunction<ValidatedTag, SIG>;
using EffectiveScalarUnivariateContinuousFunction = ContinuousFunction<EffectiveTag, Real(Real)>;
using ValidatedScalarUnivariateContinuousFunction = ContinuousFunction<ValidatedTag, Real(Real)>;

/** @brief Polymorphic interface for measurable function. */
template<class P, class SIG>
class MeasurableFunctionInterface;

/** @brief Polymorphic interface for measurable function. */
template<class P, class RES, class ARG>
class MeasurableFunctionInterface<P, RES(ARG)> {
    using SIG = RES(ARG);

  public:
    /** @brief Destroy the object and release its owned resources. */
    virtual ~MeasurableFunctionInterface() = default;
    /** @brief Return the preimage of the supplied set. */
    virtual LowerMeasurableSet<P, ARG> _preimage(OpenSet<P, RES>) const = 0;
};

template<class P, class SIG>
class MeasurableFunction;

/** @brief Compile-time predicate testing whether a type is measurable function. */
template<class F, class P, class SIG>
struct IsMeasurableFunction;

/** @brief Compile-time predicate testing whether a type is measurable function. */
template<class F, class P, class R, class A>
struct IsMeasurableFunction<F, P, R(A)> {
    using LMS = LowerMeasurableSet<P, A>;
    using OPS = OpenSet<P, R>;
    /** @brief Evaluate the compile-time or runtime predicate for the supplied arguments. */
    template<class FF, class = decltype(declval<LMS>() = preimage(declval<OPS>(), declval<FF>()))>
    static std::true_type test(int);
    /** @brief Evaluate the compile-time or runtime predicate for the supplied arguments. */
    template<class FF>
    static std::false_type test(...);
    static const bool value = decltype(test<F>(1))::value;
};

template<class F, class P, class SIG>
concept AMeasurableFunction = IsMeasurableFunction<F, P, SIG>::value;

/** @brief Polymorphic function exposing measurable-set preimages. */
template<class P, class RES, class ARG>
class MeasurableFunction<P, RES(ARG)>
    : public Handle<const MeasurableFunctionInterface<P, RES(ARG)>> {
    using SIG = RES(ARG);
    using typename Handle<const MeasurableFunctionInterface<P, SIG>>::Interface;

  public:
    using Handle<const Interface>::Handle;
    /** @brief Return the number of references to the shared representation. */
    SizeType ref_count() const { return this->_ptr.ref_count(); }
    /** @brief Construct a \c MeasurableFunction with the specified characteristics. */
    template<AMeasurableFunction<P, SIG> F>
    MeasurableFunction(const F&);
    /** @brief Return the preimage of the supplied set. */
    LowerMeasurableSet<P, ARG> preimage(const OpenSet<P, RES>& u) const { return this->_ptr->_preimage(u); }
};

/** @brief Adapter exposing measurable function. */
template<class F, class P, class SIG>
class MeasurableFunctionWrapper;

/** @brief Adapter exposing measurable function. */
template<class F, class P, class RES, class ARG>
class MeasurableFunctionWrapper<F, P, RES(ARG)>
    : public MeasurableFunctionInterface<P, RES(ARG)>,
      public F {
    using SIG = RES(ARG);

  public:
    /** @brief Return the preimage of the supplied set. */
    virtual LowerMeasurableSet<P, ARG> _preimage(const OpenSet<P, RES>& ops) const final override {
        return preimage(ops, static_cast<const F&>(*this));
    }
};

/** @brief Sequence of finite function approximations indexed by accuracy. */
template<class SIG>
class FastFanSequence;

/** @brief Sequence of finite function approximations indexed by accuracy. */
template<class RES, class ARG>
class FastFanSequence<RES(ARG)>
    : public Sequence<ValidatedContinuousFunction<RES(ARG)>> {
    using P = EffectiveTag;
    using SIG = RES(ARG);

  public:
    using Sequence<ValidatedContinuousFunction<SIG>>::Sequence;
    /** @brief Return the preimage of the supplied set. */
    LowerMeasurableSet<P, ARG> preimage(const OpenSet<P, RES>& ops) const;
    /** @brief Return the preimage of the supplied set. */
    friend LowerMeasurableSet<P, ARG> preimage(const OpenSet<P, RES>& ops, const FastFanSequence<RES(ARG)>& ffs) {
        return ffs.preimage(ops);
    }
};

static_assert(AMeasurableFunction<FastFanSequence<Real(Real)>, EffectiveTag, Real(Real)>);

/** @brief Finite approximation of a measurable function with a measure-error bound. */
template<class SIG, class PR, class PRE = PR>
class FanModel;

/** @brief Finite approximation of a measurable function with a measure-error bound. */
template<class PR, class PRE>
class FanModel<Real(Real), PR, PRE> {
    using P = ValidatedTag;
    using RES = Real;
    using ARG = Real;
    using SIG = RES(ARG);

  public:
    typedef IntervalDomainType DomainType;
    typedef FloatError<PRE> ErrorType;

    DomainType _domain;
    ValidatedContinuousFunction<SIG> _function;
    FloatError<PRE> _measure_error;

  public:
    /** @brief Construct a \c FanModel from dom, fn, err. */
    FanModel(DomainType dom, ValidatedContinuousFunction<SIG> fn, FloatError<PRE> err)
        : _domain(dom), _function(fn), _measure_error(err) {}
    /** @brief Return the domain. */
    const DomainType& domain() const { return this->_domain; }
    /** @brief Return the constrained function. */
    const ValidatedContinuousFunction<SIG>& function() const { return this->_function; }
    /** @brief Return the error bound on the set measure. */
    const ErrorType& measure_error() const { return this->_measure_error; }
    /** @brief Return the preimage of the supplied set. */
    ValidatedLowerMeasurableSet<ARG> preimage(const ValidatedOpenSet<RES>&) const;
    /** @brief Return the preimage of the supplied set. */
    friend LowerMeasurableSet<P, ARG> preimage(const OpenSet<P, RES>& ops, const FanModel<Real(Real), PR, PRE>& fm) {
        return fm.preimage(ops);
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FanModel<SIG, PR, PRE>& fm) {
        return os << "FanModel(domain=" << fm._domain << ",function=" << fm._function << ",measure_error=" << fm._measure_error << ")";
    }
};

static_assert(AMeasurableFunction<FanModel<Real(Real), DoublePrecision, DoublePrecision>, ValidatedTag, Real(Real)>);

// TODO: Move to Interval
/** @brief Return the image of the supplied domain under the function. */
template<class PR>
Interval<UpperBound<PR>> image(const Interval<UpperBound<PR>>& ivl, const ValidatedContinuousFunction<Real(Real)>& f) {
    return f(cast_singleton(ivl));
}

ValidatedOpenSet<Real> preimage(const ValidatedContinuousFunction<Real(Real)>& f, const ValidatedOpenSet<Real>& ops, IntervalDomainType dom, Accuracy acc);

template<class PR, class PRE>
auto FanModel<Real(Real), PR, PRE>::preimage(const ValidatedOpenSet<Real>& rng) const -> ValidatedLowerMeasurableSet<Real> {
    // FIXME: Make accuracy dynamic
    Accuracy acc(exp2(-4));

    const ValidatedOpenSet<Real> prim = Ariadne::preimage(this->_function, rng, this->_domain, acc);
    // Both paths of Ariadne::preimage construct this exact wrapper type.
    auto ivls_ptr = dynamic_pointer_cast<const OpenSetWrapper<UnionOfIntervals<FloatDP>, ValidatedTag, Real>>(prim.managed_pointer());
    const UnionOfIntervals<FloatDP>& ivls = *ivls_ptr;
    return LowerMeasurableSetModel(ivls, this->_measure_error);
}

/** @brief Return interval components of the preimage. */
template<class PR>
UnionOfIntervals<Float<PR>> preimage_intervals(const ValidatedContinuousFunction<Real(Real)>& f, ValidatedRegularSet<Real> rgs, Interval<Float<PR>> ivl, Accuracy acc) {
    Interval<Float<PR>> rng = cast_exact(image(ivl, f));
    if (definitely(rgs.covers(rng))) {
        return {cast_exact_interval(ivl)};
    } else if (definitely(rgs.separated(rng))) {
        return {{}, ivl.precision()};
    } else if (definitely(ivl.width() < acc.error())) {
        return {{}, ivl.precision()};
    } else {
        Pair<Interval<Float<PR>>, Interval<Float<PR>>> subivls = split(ivl);
        return join(preimage_intervals(f, rgs, subivls.first, acc), preimage_intervals(f, rgs, subivls.second, acc));
    }
}

/** @brief Return interval components of the preimage. */
template<class PR>
UnionOfIntervals<Float<PR>> preimage_intervals(const ValidatedContinuousFunction<Real(Real)>& f, ValidatedOpenSet<Real> ops, Interval<Float<PR>> ivl, Accuracy acc) {
    if (definitely(ops.covers(cast_exact(image(ivl, f))))) {
        return {ivl};
    } else if (definitely(ivl.width() < acc.error())) {
        return {{}, ivl.precision()};
    } else {
        Pair<Interval<Float<PR>>, Interval<Float<PR>>> subivls = split(ivl);
        return join(preimage_intervals(f, ops, subivls.first, acc), preimage_intervals(f, ops, subivls.second, acc));
    }
}

/** @brief Return the preimage of the supplied set. */
inline ValidatedOpenSet<Real> preimage(const ValidatedContinuousFunction<Real(Real)>& f, const ValidatedOpenSet<Real>& ops, IntervalDomainType dom, Accuracy acc) {
    DoublePrecision pr;
    Interval<FloatDP> ivl = cast_exact(Interval<UpperBound<FloatDP>>(dom, pr));
    auto rgsp = dynamic_pointer_cast<ValidatedRegularSet<Real>::Interface>(ops.managed_pointer());
    if (rgsp) {
        auto rgs = ValidatedRegularSet<Real>(rgsp);
        return ValidatedOpenSet<Real>(wrap_open(preimage_intervals(f, rgs, ivl, acc)));
    } else {
        return ValidatedOpenSet<Real>(wrap_open(preimage_intervals(f, ops, ivl, acc)));
    }
}

} // namespace Ariadne

#endif
