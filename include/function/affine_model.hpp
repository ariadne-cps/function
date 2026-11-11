/***************************************************************************
 *            function/affine_model.hpp
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

/*! @file function/affine_model.hpp
 *  \brief Affine models defined on the unit box
 */

#ifndef ARIADNE_AFFINE_MODEL_HPP
#define ARIADNE_AFFINE_MODEL_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "algebra/covector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/operations.hpp"
#include "algebra/vector.hpp"
#include "function/declarations.hpp"
#include "function/function.hpp"
#include "numeric/numeric.hpp"
#include "utility/macros.hpp"
#include "utility/pointer.hpp"

namespace Ariadne {

template<class X>
class Affine;

/** @brief Affine approximation together with its uniform error bound. */
template<class P, class F>
class AffineModel;

template<class P, class F>
class TaylorModel;

//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x) \approx \sum_{i=0}^{n-1} a_i x_i + b\f$.
template<class F>
class AffineModel<ApproximateTag, F>
    : DispatchAlgebraOperations<AffineModel<ApproximateTag, F>, FloatApproximation<typename F::PrecisionType>> {
    typedef ApproximateTag P;
    typedef typename F::PrecisionType PR;
    friend struct AlgebraOperations<AffineModel<ApproximateTag, F>, FloatApproximation<PR>>;

  public:
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef Tuple<SizeType, PrecisionType> CharacteristicsType;
    typedef Approximation<F> NumericType;
    typedef Approximation<F> CoefficientType;
    typedef Interval<FloatApproximation<PR>> RangeType;

    /** @brief Construct a \c AffineModel from n, prec. */
    explicit AffineModel(SizeType n, PrecisionType prec) : _c(0, prec), _g(n, CoefficientType(0, prec)) {}
    /** @brief Construct a \c AffineModel from n, c. */
    explicit AffineModel(SizeType n, const CoefficientType& c) : _c(c), _g(n, nul(c)) {}
    /** @brief Construct a \c AffineModel from c, g. */
    explicit AffineModel(const CoefficientType& c, const Covector<CoefficientType>& g) : _c(c), _g(g) {}
    /** @brief Construct a \c AffineModel from c, g. */
    explicit AffineModel(CoefficientType c, InitializerList<CoefficientType> g) : _c(c), _g(g) {}

    /** @brief Assign the supplied value and return this object. */
    AffineModel<ApproximateTag, F>& operator=(const CoefficientType& c) {
        this->_c = c;
        for (SizeType i = 0; i != this->_g.size(); ++i) {
            this->_g[i] = 0.0;
        }
        return *this;
    }

    /** @brief Construct a zero object with compatible characteristics. */
    AffineModel<ApproximateTag, F> create_zero() const {
        return AffineModel<ApproximateTag, F>(this->argument_size(), this->precision());
    }
    /** @brief Construct a constant with the supplied value. */
    static AffineModel<ApproximateTag, F> constant(SizeType n, const CoefficientType& c) {
        AffineModel<ApproximateTag, F> res(n, c.precision());
        res._c = c;
        return res;
    }
    /** @brief Construct or return the selected coordinate. */
    static AffineModel<ApproximateTag, F> coordinate(SizeType n, SizeType j, PrecisionType pr) {
        AffineModel<ApproximateTag, F> res(n, pr);
        res._g[j] = 1;
        return res;
    }
    /** @brief Return the affine scaling for the selected coordinate. */
    static AffineModel<ApproximateTag, F> scaling(SizeType n, SizeType j, const IntervalDomainType& codom, PrecisionType precision);
    /** @brief Return the affine scalings of the domain coordinates. */
    static Vector<AffineModel<ApproximateTag, F>> scalings(const BoxDomainType& codom, PrecisionType precision);

    /** @brief Return the number of scalar arguments. */
    SizeType argument_size() const { return this->_g.size(); }
    /** @brief Return the working precision. */
    PrecisionType precision() const { return this->_c.precision(); }
    /** @brief Return the characteristics required to construct compatible values. */
    CharacteristicsType characteristics() const { return std::make_pair(this->argument_size(), this->_c.precision()); }
    /** @brief Return the linear coefficients of the affine object. */
    const Covector<CoefficientType>& a() const { return this->_g; }
    /** @brief Return the constant part of the affine object. */
    const CoefficientType& b() const { return this->_c; }

    /** @brief Return the stored value. */
    const CoefficientType& value() const { return this->_c; }
    /** @brief Return the vector of first partial derivatives. */
    const Covector<CoefficientType>& gradient() const { return this->_g; }
    /** @brief Return the vector of first partial derivatives. */
    const CoefficientType& gradient(SizeType i) const { return this->_g[i]; }
    /** @brief Return the selected element or component. */
    CoefficientType& operator[](SizeType i) { return this->_g[i]; }
    /** @brief Return the selected element or component. */
    const CoefficientType& operator[](SizeType i) const { return this->_g[i]; }

    /** @brief Return the range enclosure. */
    RangeType range() const;

    /** @brief Set the value. */
    Void set_value(const CoefficientType& c) { _c = c; }
    /** @brief Set the gradient. */
    Void set_gradient(SizeType j, const CoefficientType& g) { _g[j] = g; }
    /** @brief Set the gradient. */
    Void set_gradient(SizeType j, const ApproximateNumber& g) { _g[j] = g; }
    /** @brief Evaluate at the supplied argument. */
    template<class X>
    X evaluate(const Vector<X>& v) const;
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const AffineModel<ApproximateTag, F>& f) { return f._write(os); }

  private:
    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream&) const;

  private:

  public:
    /** @brief Compose the outer object with the inner object. */
    static AffineModel<P, F> _compose(const ScalarMultivariateFunction<P>& f, const Vector<AffineModel<P, F>>& g);
    /** @brief Compose the outer object with the inner object. */
    static Vector<AffineModel<P, F>> _compose(const VectorMultivariateFunction<P>& f, const Vector<AffineModel<P, F>>& g);

  private:
    CoefficientType _c;
    Covector<CoefficientType> _g;
};

//! An affine expression \f$f:[-1,+1]^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b \pm e\f$.
template<class F>
class AffineModel<ValidatedTag, F>
    : DispatchAlgebraOperations<AffineModel<ValidatedTag, F>, FloatBounds<typename F::PrecisionType>> {
    typedef ValidatedTag P;
    typedef typename F::PrecisionType PR;
    friend struct AlgebraOperations<AffineModel<ValidatedTag, F>, FloatBounds<PR>>;

  public:
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef Tuple<SizeType, PrecisionType> CharacteristicsType;
    typedef Float<PR> CoefficientType;
    typedef FloatError<PR> ErrorType;
    typedef FloatBounds<PR> NumericType;
    typedef AffineModel<ValidatedTag, F> AffineModelType;
    typedef Interval<FloatUpperBound<PR>> RangeType;

    /** @brief Construct a \c AffineModel from n, prec. */
    explicit AffineModel(SizeType n, PrecisionType prec) : _c(0, prec), _g(n, CoefficientType(0, prec)), _e(0u, prec) {}
    /** @brief Construct a \c AffineModel from n, c. */
    explicit AffineModel(SizeType n, const CoefficientType& c) : _c(c), _g(n, nul(c)), _e(nul(c)) {}
    /** @brief Construct a \c AffineModel from c, g, e. */
    explicit AffineModel(const CoefficientType& c, const Covector<CoefficientType>& g, const ErrorType& e) : _c(c), _g(g), _e(e) {}
    /** @brief Construct a \c AffineModel from c, g. */
    explicit AffineModel(CoefficientType c, InitializerList<CoefficientType> g) : _c(c), _g(g), _e(0u, c.precision()) {}

    /** @brief Construct a \c AffineModel from affine. */
    explicit AffineModel(const Affine<FloatBounds<PR>>& affine);

    /** @brief Construct a \c AffineModel with the specified characteristics. */
    AffineModel(const TaylorModel<ValidatedTag, F>&);

    /** @brief Assign the supplied value and return this object. */
    AffineModelType& operator=(const CoefficientType& c) {
        this->_c = c;
        for (SizeType i = 0; i != this->_g.size(); ++i) {
            this->_g[i] = 0;
        }
        this->_e = 0u;
        return *this;
    }

    /** @brief Return the characteristics required to construct compatible values. */
    CharacteristicsType characteristics() const { return std::make_pair(this->argument_size(), this->_c.precision()); }
    /** @brief Construct a zero object with compatible characteristics. */
    AffineModel<ValidatedTag, F> create_zero() const {
        return AffineModel<ValidatedTag, F>(this->argument_size(), this->precision());
    }
    /** @brief Construct a constant with the supplied value. */
    static AffineModel<ValidatedTag, F> constant(SizeType n, const CoefficientType& c) {
        AffineModel<ValidatedTag, F> res(n, c.precision());
        res._c = c;
        return res;
    }
    /** @brief Construct or return the selected coordinate. */
    static AffineModel<ValidatedTag, F> coordinate(SizeType n, SizeType j, PrecisionType prec) {
        AffineModel<ValidatedTag, F> res(n, prec);
        res._g[j] = 1;
        return res;
    }

    /** @brief Return the affine scaling for the selected coordinate. */
    static AffineModel<ValidatedTag, F> scaling(SizeType n, SizeType j, const IntervalDomainType& codom, PrecisionType precision);
    /** @brief Return the affine scalings of the domain coordinates. */
    static Vector<AffineModel<ValidatedTag, F>> scalings(const BoxDomainType& codom, PrecisionType precision);

    /** @brief Return the number of scalar arguments. */
    SizeType argument_size() const { return this->_g.size(); }
    /** @brief Return the working precision. */
    PrecisionType precision() const { return this->_c.precision(); }
    /** @brief Return the linear coefficients of the affine object. */
    const Covector<CoefficientType>& a() const { return this->_g; }
    /** @brief Return the constant part of the affine object. */
    const CoefficientType& b() const { return this->_c; }
    /** @brief Return the error bound of the affine model. */
    const ErrorType& e() const { return this->_e; }

    /** @brief Return the stored value. */
    const CoefficientType& value() const { return this->_c; }
    /** @brief Return the vector of first partial derivatives. */
    const Covector<CoefficientType>& gradient() const { return this->_g; }
    /** @brief Return the vector of first partial derivatives. */
    const CoefficientType& gradient(SizeType i) const { return this->_g[i]; }
    /** @brief Return the error bound. */
    const ErrorType& error() const { return this->_e; }
    /** @brief Return the stored value. */
    CoefficientType& value() { return this->_c; }
    /** @brief Return the vector of first partial derivatives. */
    CoefficientType& gradient(SizeType i) { return this->_g[i]; }
    /** @brief Return the error bound. */
    ErrorType& error() { return this->_e; }
    /** @brief Return the selected element or component. */
    CoefficientType& operator[](SizeType i) { return this->_g[i]; }
    /** @brief Return the selected element or component. */
    const CoefficientType& operator[](SizeType i) const { return this->_g[i]; }

    /** @brief Return the range enclosure. */
    RangeType range() const;

    /** @brief Set the value. */
    Void set_value(const CoefficientType& c) { _c = c; }
    /** @brief Set the gradient. */
    Void set_gradient(SizeType j, const CoefficientType& g) { _g[j] = g; }
    /** @brief Set the gradient. */
    Void set_gradient(SizeType j, const Dyadic& g) { _g[j] = g; }
    /** @brief Set the error. */
    Void set_error(const ErrorType& e) { _e = e; }
    /** @brief Set the error. */
    Void set_error(Nat m) { _e = m; }

    /** @brief Evaluate at the supplied argument. */
    template<class X>
    X evaluate(const Vector<X>& v) const;
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const AffineModel<ValidatedTag, F>& f) { return f._write(os); }

  private:
    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream&) const;

  private:

  public:
    /** @brief Compose the outer object with the inner object. */
    static AffineModel<P, F> _compose(const ScalarMultivariateFunction<P>& f, const Vector<AffineModel<P, F>>& g);
    /** @brief Compose the outer object with the inner object. */
    static Vector<AffineModel<P, F>> _compose(const VectorMultivariateFunction<P>& f, const Vector<AffineModel<P, F>>& g);

  private:
    CoefficientType _c;
    Covector<CoefficientType> _g;
    ErrorType _e;
};

template<class F>
template<class X>
X AffineModel<ValidatedTag, F>::evaluate(const Vector<X>& v) const {
    X r = v.zero_element() + static_cast<X>(this->_c);
    for (SizeType i = 0; i != this->_g.size(); ++i) {
        r += X(this->_g[i]) * v[i];
    }
    r += FloatBounds<PR>(-_e, +_e);
    return r;
}

template<class F>
template<class X>
X AffineModel<ApproximateTag, F>::evaluate(const Vector<X>& v) const {
    X r = v.zero_element() + static_cast<X>(this->_c);
    for (SizeType i = 0; i != this->_g.size(); ++i) {
        r += X(this->_g[i]) * v[i];
    }
    return r;
}

/** @brief Construct an affine model from the supplied function or coefficients. */
template<class P, class PR>
inline AffineModel<P, RawFloatType<PR>> affine_model(const Affine<Number<P>>& affine, PR precision) {
    return AffineModel<P, RawFloatType<PR>>(affine, precision);
}
/** @brief Construct an affine model from the supplied function or coefficients. */
template<class F>
inline AffineModel<ApproximateTag, F> affine_model(const Affine<Approximation<F>>& affine) {
    return AffineModel<ApproximateTag, F>(affine);
}
/** @brief Construct an affine model from the supplied function or coefficients. */
template<class F>
inline AffineModel<ValidatedTag, F> affine_model(const Affine<Bounds<F>>& affine) {
    return AffineModel<ValidatedTag, F>(affine);
}

template<class P, class F>
inline AffineModel<P, F> affine_model(const TaylorModel<P, F>& taylor_model);
template<class P, class F>
inline Vector<AffineModel<P, F>> affine_models(const Vector<TaylorModel<P, F>>& taylor_models);

// DEPRECATED
template<class P, class PR>
AffineModel<P, RawFloatType<PR>> affine_model(const BoxDomainType& domain, const ScalarMultivariateFunction<P>& function, PR precision);
template<class P, class PR>
inline Vector<AffineModel<P, RawFloatType<PR>>> affine_models(const BoxDomainType& domain, const VectorMultivariateFunction<P>& function, PR precision);

/** @brief Construct an affine model from the supplied function or coefficients. */
template<class PR>
inline AffineModel<ValidatedTag, RawFloatType<PR>> affine_model(const BoxDomainType& domain, const ScalarMultivariateFunction<EffectiveTag>& function, PR precision) {
    return affine_model(domain, ScalarMultivariateFunction<ValidatedTag>(function), precision);
}
/** @brief Construct affine models for the component functions. */
template<class PR>
inline Vector<AffineModel<ValidatedTag, RawFloatType<PR>>> affine_models(const BoxDomainType& domain, const VectorMultivariateFunction<EffectiveTag>& function, PR precision) {
    return affine_models(domain, VectorMultivariateFunction<ValidatedTag>(function), precision);
}

/** @brief Construct an affine model from the supplied function or coefficients. */
template<class P, class F>
inline AffineModel<P, F> affine_model(const TaylorModel<P, F>& taylor_model) {
    return AffineModel<P, F>(taylor_model);
}

/** @brief Construct affine models for the component functions. */
template<class P, class F>
Vector<AffineModel<P, F>> affine_models(const Vector<TaylorModel<P, F>>& models) {
    SizeType rs = models.size();
    SizeType as = models.zero_element().argument_size();
    auto zero = models.zero_element().value();

    Vector<AffineModel<ValidatedTag, F>> result(rs, AffineModel<ValidatedTag, F>(as, zero));
    for (SizeType i = 0; i != rs; ++i) {
        result[i] = affine_model(models[i]);
    }
    return result;
}

/** @brief Construct affine models for the component functions. */
template<class P, class PR>
Vector<AffineModel<P, RawFloatType<PR>>> affine_models(const BoxDomainType& domain, const VectorMultivariateFunction<P>& function, PR precision) {
    typedef RawFloat<PR> F;
    SizeType rs = function.result_size();
    SizeType as = function.argument_size();

    Vector<AffineModel<P, F>> result(rs, AffineModel<P, F>(as, precision));
    for (SizeType i = 0; i != rs; ++i) {
        result[i] = affine_model(domain, function[i], precision);
    }
    return result;
}

/** @brief Return the component range enclosures. */
template<class P, class F>
Vector<typename AffineModel<P, F>::RangeType> ranges(const Vector<AffineModel<P, F>>& f) {
    return Vector<typename AffineModel<P, F>::RangeType>(f.size(), [&](SizeType i) { return f[i].range(); });
    // Vector<typename AffineModel<P,F>::RangeType> r(f.size(),f.zero_element().range()); for(SizeType i=0; i!=f.size(); ++i) { r[i]=f[i].range(); } return r;
}

/** @brief Compose the outer object with the inner object. */
template<class P, class F>
Vector<AffineModel<P, F>> compose(const VectorMultivariateFunction<P>& f, const Vector<AffineModel<P, F>>& g) {
    return AffineModel<P, F>::_compose(f, g);
}

/** @brief Compose the outer object with the inner object. */
template<class P, class F>
AffineModel<P, F> compose(const ScalarMultivariateFunction<P>& f, const Vector<AffineModel<P, F>>& g) {
    return AffineModel<P, F>::_compose(f, g);
}

} // namespace Ariadne

#endif /* ARIADNE_AFFINE_MODEL_HPP */
