/***************************************************************************
 *            function/function_interface.hpp
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

/*! @file function/function_interface.hpp
 *  \brief Interface for functions for which derivatives can be computed.
 */

#ifndef ARIADNE_FUNCTION_INTERFACE_HPP
#define ARIADNE_FUNCTION_INTERFACE_HPP

#include <iosfwd>

#include "foundation/representation.hpp"
#include "function/declarations.hpp"
#include "function/formula.decl.hpp"
#include "function/function.decl.hpp"
#include "function/function_traits.hpp"

namespace Ariadne {

static const Int SMOOTH = 255;

template<class P, class SIG>
OutputStream& operator<<(OutputStream& os, const Representation<Function<P, SIG>>& f);

template<class P, class SIG>
class FunctionInterface;

template<class P, class... ARGS>
using ScalarFunctionInterface = FunctionInterface<P, RealScalar(ARGS...)>;
template<class P, class... ARGS>
using VectorFunctionInterface = FunctionInterface<P, RealVector(ARGS...)>;
template<class P>
using ScalarUnivariateFunctionInterface = ScalarFunctionInterface<P, RealScalar>;
template<class P>
using VectorUnivariateFunctionInterface = VectorFunctionInterface<P, RealScalar>;
template<class P>
using ScalarMultivariateFunctionInterface = ScalarFunctionInterface<P, RealVector>;
template<class P>
using VectorMultivariateFunctionInterface = VectorFunctionInterface<P, RealVector>;

/** @brief Polymorphic interface for vector of function. */
template<class P, class... ARGS>
class VectorOfFunctionInterface {
};

/** @brief Polymorphic interface for vector of function. */
template<class... ARGS>
class VectorOfFunctionInterface<ApproximateTag, ARGS...> {
  public:
    /** @brief Destroy the object and release its owned resources. */
    virtual ~VectorOfFunctionInterface() = default;
    /** @brief Return the requested element or stored value. */
    virtual ScalarFunctionInterface<ApproximateTag, ARGS...>* _get(SizeType i) const = 0;
};

/** @brief Polymorphic interface for vector of function. */
template<class... ARGS>
class VectorOfFunctionInterface<ValidatedTag, ARGS...>
    : public virtual VectorOfFunctionInterface<ApproximateTag, ARGS...> {
  public:
    /** @brief Return the requested element or stored value. */
    virtual ScalarFunctionInterface<ValidatedTag, ARGS...>* _get(SizeType i) const override = 0;
};

/** @brief Polymorphic interface for vector of function. */
template<class... ARGS>
class VectorOfFunctionInterface<EffectiveTag, ARGS...>
    : public virtual VectorOfFunctionInterface<ValidatedTag, ARGS...> {
  public:
    /** @brief Return the requested element or stored value. */
    virtual ScalarFunctionInterface<EffectiveTag, ARGS...>* _get(SizeType i) const override = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for functions \f$\R^n\rightarrow\R^m\f$ which can be evaluated.
template<class P, class SIG>
class FunctionInterface;

template<class T>
struct Representation;
/** @brief Write the value to the output stream. */
template<class P, class S>
OutputStream& operator<<(OutputStream& os, const Representation<Function<P, S>>& f);

/** @brief Polymorphic interface for function. */
template<class SIG>
class FunctionInterface<Void, SIG> {
    template<class X>
    using Argument = typename SignatureTraits<SIG>::template Argument<X>;
    template<class X>
    using Result = typename SignatureTraits<SIG>::template Result<X>;

  public:
    typedef typename SignatureTraits<SIG>::DomainType DomainType;
    typedef typename SignatureTraits<SIG>::CodomainType CodomainType;
    typedef typename SignatureTraits<SIG>::ArgumentSizeType ArgumentSizeType;
    typedef typename SignatureTraits<SIG>::ResultSizeType ResultSizeType;
    typedef typename SignatureTraits<SIG>::ArgumentIndexType ArgumentIndexType;
    typedef typename SignatureTraits<SIG>::ResultIndexType ResultIndexType;

    /** @brief Destroy the object and release its owned resources. */
    virtual ~FunctionInterface() = default;
    /** @brief Return the number of scalar arguments. */
    virtual ArgumentSizeType argument_size() const = 0;
    /** @brief Return the number of result components. */
    virtual ResultSizeType result_size() const = 0;
    //    virtual DomainType const domain() const = 0;
    //    virtual CodomainType const codomain() const = 0;

    /** @brief Write the representation of the object. */
    virtual OutputStream& _repr(OutputStream& os) const = 0;
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const = 0;

  public:
    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    virtual FunctionInterface<Void, SIG>* _clone() const = 0;
    /** @brief Return a copy of the object. */
    inline FunctionInterface<Void, SIG>* _copy() const { return this->_clone(); }

  public:
    /** @brief Write the value to the output stream. */
    friend inline OutputStream& operator<<(OutputStream& os, const FunctionInterface<Void, SIG>& f) { return f._write(os); }
    /** @brief Write the value to the output stream. */
    template<class P, class S>
    friend OutputStream& operator<<(OutputStream& os, const Representation<Function<P, S>>& f);
};

//! \ingroup FunctionModule
//! \brief Interface for functions \f$\R^n\rightarrow\R^m\f$ which can only be evaluated approximately.
//! \sa \ref FunctionInterface "FunctionInterface<P,SIG>"
template<class SIG>
class FunctionInterface<ApproximateTag, SIG>
    : public virtual FunctionInterface<Void, SIG> {
    using P = ApproximateTag;

  public:
    template<class X>
    using Argument = typename SignatureTraits<SIG>::template Argument<X>;
    template<class X>
    using Result = typename SignatureTraits<SIG>::template Result<X>;
    using typename FunctionInterface<Void, SIG>::ArgumentIndexType;

  public:
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<ApproximateNumber> _call(const Argument<ApproximateNumber>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<FloatDPApproximation> _call(const Argument<FloatDPApproximation>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<FloatMPApproximation> _call(const Argument<FloatMPApproximation>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<Differential<FloatDPApproximation>> _call(const Argument<Differential<FloatDPApproximation>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<Differential<FloatMPApproximation>> _call(const Argument<Differential<FloatMPApproximation>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<TaylorModel<ApproximateTag, FloatDP>> _call(const Argument<TaylorModel<ApproximateTag, FloatDP>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<TaylorModel<ApproximateTag, FloatMP>> _call(const Argument<TaylorModel<ApproximateTag, FloatMP>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<Formula<ApproximateNumber>> _call(const Argument<Formula<ApproximateNumber>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<ElementaryAlgebra<ApproximateNumber>> _call(const Argument<ElementaryAlgebra<ApproximateNumber>>& x) const = 0;

    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    virtual FunctionInterface<P, SIG>* _clone() const = 0;
    /** @brief Return a copy of the object. */
    inline FunctionInterface<P, SIG>* _copy() const { return this->_clone(); }
    /** @brief Return the derivative with respect to the selected argument. */
    virtual FunctionInterface<P, SIG>* _derivative(ArgumentIndexType i) const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for functions \f$\R^n\rightarrow\R\f$ which can be evaluated with guaranteed bounds on the error.
//! \sa \ref FunctionInterface "FunctionInterface<P,SIG>"
template<class SIG>
class FunctionInterface<ValidatedTag, SIG>
    : public virtual FunctionInterface<ApproximateTag, SIG> {
    using P = ValidatedTag;
    using AP = ApproximateTag;

  public:
    template<class X>
    using Argument = typename SignatureTraits<SIG>::template Argument<X>;
    template<class X>
    using Result = typename SignatureTraits<SIG>::template Result<X>;
    using typename FunctionInterface<Void, SIG>::ArgumentIndexType;

  public:
    using FunctionInterface<AP, SIG>::_call;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<ValidatedNumber> _call(const Argument<ValidatedNumber>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<FloatDPBounds> _call(const Argument<FloatDPBounds>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<FloatMPBounds> _call(const Argument<FloatMPBounds>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<Differential<FloatDPBounds>> _call(const Argument<Differential<FloatDPBounds>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<Differential<FloatMPBounds>> _call(const Argument<Differential<FloatMPBounds>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<TaylorModel<ValidatedTag, FloatDP>> _call(const Argument<TaylorModel<ValidatedTag, FloatDP>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<TaylorModel<ValidatedTag, FloatMP>> _call(const Argument<TaylorModel<ValidatedTag, FloatMP>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<TaylorModel<ValidatedTag, FloatDPBounds>> _call(const Argument<TaylorModel<ValidatedTag, FloatDPBounds>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<TaylorModel<ValidatedTag, FloatMPBounds>> _call(const Argument<TaylorModel<ValidatedTag, FloatMPBounds>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<TaylorModel<ValidatedTag, FloatDPUpperInterval>> _call(const Argument<TaylorModel<ValidatedTag, FloatDPUpperInterval>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<TaylorModel<ValidatedTag, FloatMPUpperInterval>> _call(const Argument<TaylorModel<ValidatedTag, FloatMPUpperInterval>>& x) const = 0;

    /** @brief Evaluate with the supplied arguments. */
    virtual Result<Formula<ValidatedNumber>> _call(const Argument<Formula<ValidatedNumber>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<ElementaryAlgebra<ValidatedNumber>> _call(const Argument<ElementaryAlgebra<ValidatedNumber>>& x) const = 0;

    /** @brief Evaluate with the supplied arguments. */
    virtual Result<ScalarMultivariateFunction<ValidatedTag>> _call(const Argument<ScalarMultivariateFunction<ValidatedTag>>& x) const = 0;

    /** @brief Evaluate with the supplied arguments. */
    inline Result<FloatDPBounds> _call(const Argument<FloatDP>& x) const {
        return this->_call(Argument<FloatDPBounds>(x));
    }
    /** @brief Evaluate with the supplied arguments. */
    inline Result<FloatMPBounds> _call(const Argument<FloatMP>& x) const {
        return this->_call(Argument<FloatMPBounds>(x));
    }

    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    virtual FunctionInterface<P, SIG>* _clone() const = 0;
    /** @brief Return a copy of the object. */
    inline FunctionInterface<P, SIG>* _copy() const { return this->_clone(); }
    /** @brief Return the derivative with respect to the selected argument. */
    virtual FunctionInterface<P, SIG>* _derivative(ArgumentIndexType i) const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\R^n\rightarrow\R\f$ which can be evaluated exactly.
//! \sa \ref FunctionInterface "FunctionInterface<P,SIG>"
template<class SIG>
class FunctionInterface<EffectiveTag, SIG>
    : public virtual FunctionInterface<ValidatedTag, SIG> {
    using WP = ValidatedTag;
    using P = EffectiveTag;
    using D = typename SignatureTraits<SIG>::DomainType;
    using C = typename SignatureTraits<SIG>::CodomainType;

  public:
    template<class X>
    using Argument = typename SignatureTraits<SIG>::template Argument<X>;
    template<class X>
    using Result = typename SignatureTraits<SIG>::template Result<X>;
    using typename FunctionInterface<Void, SIG>::ArgumentIndexType;

  public:
    using FunctionInterface<WP, SIG>::_call;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<EffectiveNumber> _call(const Argument<EffectiveNumber>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<Real> _call(const Argument<Real>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<ElementaryAlgebra<Real>> _call(const Argument<ElementaryAlgebra<Real>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<Formula<Real>> _call(const Argument<Formula<Real>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<ElementaryAlgebra<EffectiveNumber>> _call(const Argument<ElementaryAlgebra<EffectiveNumber>>& x) const = 0;
    /** @brief Evaluate with the supplied arguments. */
    virtual Result<Formula<EffectiveNumber>> _call(const Argument<Formula<EffectiveNumber>>& x) const = 0;

    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    virtual FunctionInterface<P, SIG>* _clone() const = 0;
    /** @brief Return a copy of the object. */
    inline FunctionInterface<P, SIG>* _copy() const { return this->_clone(); }
    /** @brief Return the derivative with respect to the selected argument. */
    virtual FunctionInterface<P, SIG>* _derivative(ArgumentIndexType i) const = 0;
};

/** @brief Polymorphic interface for vector. */
template<class I>
class VectorInterface {
    /** @brief Destroy the object and release its owned resources. */
    ~VectorInterface() = default;
    /** @brief Return the requested element or stored value. */
    virtual I* _get(SizeType i) const = 0;
};

/** @brief Polymorphic interface for function factory. */
template<class P>
class FunctionFactoryInterface;

/** @brief Polymorphic interface for function factory. */
template<>
class FunctionFactoryInterface<ValidatedTag> {
    using P = ValidatedTag;
    typedef BoxDomainType DomainType;

  public:
    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    virtual FunctionFactoryInterface<P>* clone() const = 0;
    /** @brief Return a copy of the object. */
    inline FunctionFactoryInterface<P>* _copy() const { return this->clone(); }
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const = 0;
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    inline ScalarMultivariateFunction<P> create(const BoxDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const;
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    inline VectorMultivariateFunction<P> create(const BoxDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const;
    /** @brief Construct a zero object with compatible characteristics. */
    inline ScalarMultivariateFunction<P> create_zero(const BoxDomainType& domain) const;
    /** @brief Construct the identity function. */
    inline VectorMultivariateFunction<P> create_identity(const BoxDomainType& domain) const;

  private:
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual ScalarMultivariateFunctionInterface<P>* _create(const BoxDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const = 0;
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual VectorMultivariateFunctionInterface<P>* _create(const BoxDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const = 0;

  public:
    /** @brief Write the value to the output stream. */
    friend inline OutputStream& operator<<(OutputStream& os, const FunctionFactoryInterface<P>& factory) {
        return factory._write(os);
    }

  public:
    /** @brief Destroy the object and release its owned resources. */
    virtual ~FunctionFactoryInterface() = default;
};

} // namespace Ariadne

#endif /* ARIADNE_FUNCTION_INTERFACE_HPP */
