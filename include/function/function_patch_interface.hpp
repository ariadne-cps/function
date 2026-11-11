/***************************************************************************
 *            function/function_patch_interface.hpp
 *
 *  Copyright  2011-20  Pieter Collins
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

/*! @file function/function_patch_interface.hpp
 *  \brief Interface for functions on bounded sets.
 */

#ifndef ARIADNE_FUNCTION_PATCH_INTERFACE_HPP
#define ARIADNE_FUNCTION_PATCH_INTERFACE_HPP

#include "../function/domain.hpp"
#include "../function/function.decl.hpp"
#include "../function/function_interface.hpp"
#include "algebra/algebra_interface.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/numeric.hpp"
#include "numeric/operators.hpp"

namespace Ariadne {

typedef Interval<FloatDPUpperBound> IntervalRangeType;

template<class P, class SIG>
class FunctionPatch;

template<class P>
class FunctionPatchFactoryInterface;
/** @brief Polymorphic interface for function patch creator. */
template<class P, class... ARGS>
class FunctionPatchCreatorInterface;

template<class P, class SIG>
class FunctionPatchInterface;
template<class P, class... ARGS>
using ScalarFunctionPatchInterface = FunctionPatchInterface<P, RealScalar(ARGS...)>;
template<class P, class... ARGS>
using VectorFunctionPatchInterface = FunctionPatchInterface<P, RealVector(ARGS...)>;
template<class P>
using ScalarUnivariateFunctionPatchInterface = FunctionPatchInterface<P, RealScalar(RealScalar)>;
template<class P>
using VectorUnivariateFunctionPatchInterface = FunctionPatchInterface<P, RealVector(RealScalar)>;
template<class P>
using ScalarMultivariateFunctionPatchInterface = FunctionPatchInterface<P, RealScalar(RealVector)>;
template<class P>
using VectorMultivariateFunctionPatchInterface = FunctionPatchInterface<P, RealVector(RealVector)>;

using ValidatedFunctionPatchFactoryInterface = FunctionPatchFactoryInterface<ValidatedTag>;

/** @brief Polymorphic interface for function patch algebra. */
template<class P, class SIG>
class FunctionPatchAlgebraInterface;

/** @brief Polymorphic interface for function patch algebra. */
template<class P, class... ARGS>
class FunctionPatchAlgebraInterface<P, RealScalar(ARGS...)>
    : public virtual ElementaryAlgebraInterface<Number<P>> {};

/** @brief Polymorphic interface for function patch algebra. */
template<class P, class... ARGS>
class FunctionPatchAlgebraInterface<P, RealVector(ARGS...)> {
    using RES = RealVector;
    using SIG = RES(ARGS...);

  public:
    typedef SizeType ResultIndexType;
    /** @brief Return the selected component as a function patch. */
    virtual ScalarFunctionPatchInterface<P, ARGS...>* _patch_get(ResultIndexType) const = 0;
    /** @brief Replace the selected value. */
    virtual Void _set(ResultIndexType, const ScalarFunctionPatchInterface<P, ARGS...>&) = 0;
    /** @brief Add the supplied value or component. */
    virtual Void _adjoin(const ScalarFunctionPatchInterface<P, ARGS...>& f2) = 0;
    /** @brief Join the operands in their natural order. */
    virtual VectorFunctionPatchInterface<P, ARGS...>* _join(const VectorFunctionPatchInterface<P, ARGS...>& f2) const = 0;
    /** @brief Combine the operands over their respective arguments. */
    virtual VectorFunctionPatchInterface<P, ARGS...>* _combine(const VectorFunctionPatchInterface<P, ARGS...>& f2) const = 0;
};

/** @brief Polymorphic interface for function patch. */
template<class P, class SIG>
class FunctionPatchInterface
    : public virtual FunctionInterface<P, SIG>,
      public virtual FunctionPatchAlgebraInterface<P, SIG> {
    using RES = typename SignatureTraits<SIG>::ResultKind;
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    static_assert(Same<RES, RealScalar> or Same<RES, RealVector>);
    static_assert(Same<ARG, RealScalar> or Same<ARG, RealVector>);

  public:
    typedef typename SignatureTraits<SIG>::BoundedCodomainType CodomainType;
    typedef typename SignatureTraits<SIG>::BoundedDomainType DomainType;
    typedef typename SignatureTraits<SIG>::BoundedRangeType RangeType;
    typedef PositiveValidatedUpperNumber NormType;
    typedef PositiveValidatedUpperNumber ErrorType;
    typedef ExactNumber CoefficientType;
    typedef ElementSizeType<DomainType> ArgumentSizeType;
    typedef ElementSizeType<CodomainType> ResultSizeType;
    typedef ElementIndexType<DomainType> ArgumentIndexType;

    template<class X>
    using Argument = typename SignatureTraits<SIG>::template Argument<X>;
    template<class X>
    using Result = typename SignatureTraits<SIG>::template Result<X>;

  public:
    /** @brief Return the component error bounds. */
    virtual const Result<ErrorType> _errors() const = 0;
    //    virtual CoefficientType const _gradient_value(ArgumentIndexType) const = 0;
    /** @brief Return the error bound. */
    virtual const ErrorType _error() const = 0;
    /** @brief Return the norm. */
    virtual const NormType _norm() const = 0;
    /** @brief Return the range enclosure. */
    virtual const RangeType _range() const = 0;

    /** @brief Discard the error bound of the model. */
    virtual Void _clobber() = 0;
    /** @brief Test whether the first enclosure is contained in the second. */
    virtual Bool _refines(const FunctionPatchInterface<P, SIG>& fp) const = 0;

    /** @brief Return a factory for compatible objects. */
    virtual FunctionPatchFactoryInterface<P>* _factory() const = 0;

    /** @brief Return the domain. */
    virtual const DomainType domain() const = 0;
    /** @brief Return the codomain. */
    virtual const CodomainType codomain() const = 0;

    /** @brief Evaluate without checking that the argument belongs to the domain. */
    virtual Result<Number<P>> _unchecked_evaluate(const Argument<Number<P>>& x) const = 0;
    /** @brief Evaluate without checking that the argument belongs to the domain. */
    virtual Result<CanonicalNumericType<P, DP>> _unchecked_evaluate(const Argument<CanonicalNumericType<P, DP>>& x) const = 0;
    /** @brief Evaluate without checking that the argument belongs to the domain. */
    virtual Result<CanonicalNumericType<P, MultiplePrecision>> _unchecked_evaluate(const Argument<CanonicalNumericType<P, MultiplePrecision>>& x) const = 0;
    /** @brief Fix the selected argument to the supplied value. */
    virtual FunctionPatchInterface<P, SIG>* _partial_evaluate(SizeType j, const Number<P>& c) const = 0;

    /** @brief Compose the outer object with the inner object. */
    virtual ScalarFunctionPatchInterface<P, ARG>* _compose(const ScalarFunction<P, RES>& f) const = 0;
    /** @brief Compose the outer object with the inner object. */
    virtual VectorFunctionPatchInterface<P, ARG>* _compose(const VectorFunction<P, RES>& f) const = 0;
    /** @brief Compose without checking domain compatibility. */
    virtual ScalarFunctionPatchInterface<P, ARG>* _unchecked_compose(const ScalarFunction<P, RES>& f) const = 0;
    /** @brief Compose without checking domain compatibility. */
    virtual VectorFunctionPatchInterface<P, ARG>* _unchecked_compose(const VectorFunction<P, RES>& f) const = 0;

    /** @brief Return a dynamically allocated copy preserving the concrete type using the concrete implementation. */
    virtual FunctionPatchInterface<P, SIG>* _concrete_clone() const = 0;
    /** @brief Create a function patch with compatible characteristics. */
    virtual FunctionPatchInterface<P, SIG>* _patch_create() const = 0;
    /** @brief Embed in the larger argument space. */
    virtual FunctionPatchInterface<P, SIG>* _embed(const DomainType& d1, const DomainType& d2) const = 0;
    /** @brief Return the object restricted to the supplied domain. */
    virtual FunctionPatchInterface<P, SIG>* _restriction(const DomainType& d) const = 0;

    /** @brief Return the derivative with respect to the selected argument using the concrete implementation. */
    virtual FunctionPatchInterface<P, SIG>* _concrete_derivative(ArgumentIndexType j) const = 0;
    /** @brief Return an antiderivative with respect to the selected argument. */
    virtual FunctionPatchInterface<P, SIG>* _antiderivative(ArgumentIndexType j) const = 0;
    /** @brief Return an antiderivative with respect to the selected argument. */
    virtual FunctionPatchInterface<P, SIG>* _antiderivative(ArgumentIndexType j, Number<P> c) const = 0;

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FunctionPatchInterface<P, SIG>& f) {
        return os << static_cast<const FunctionInterface<P, SIG>&>(f);
    }
};

/** @brief Polymorphic interface for function patch factory. */
template<class P>
class FunctionPatchFactoryInterface {
    typedef RealScalar SARG;
    typedef RealVector VARG;

  public:
    typedef IntervalDomainType ScalarDomainType;
    typedef BoxDomainType VectorDomainType;
    /** @brief Destroy the object and release its owned resources. */
    virtual ~FunctionPatchFactoryInterface() = default;
    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    virtual FunctionPatchFactoryInterface<P>* clone() const = 0;
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const = 0;
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FunctionPatchFactoryInterface<P>& factory) {
        factory._write(os);
        return os;
    }

  private:

  public:
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual ScalarFunctionPatchInterface<P, SARG>* _create(const ScalarDomainType& domain, const ScalarUnivariateFunctionInterface<P>& function) const = 0;
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual VectorFunctionPatchInterface<P, SARG>* _create(const ScalarDomainType& domain, const VectorUnivariateFunctionInterface<P>& function) const = 0;
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual ScalarFunctionPatchInterface<P, VARG>* _create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const = 0;
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual VectorFunctionPatchInterface<P, VARG>* _create(const VectorDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const = 0;

    /** @brief Construct a zero object with compatible characteristics. */
    virtual ScalarFunctionPatchInterface<P, SARG>* _create_zero(const ScalarDomainType& domain) const = 0;
    /** @brief Construct a constant object with the supplied value. */
    virtual ScalarFunctionPatchInterface<P, SARG>* _create_constant(const ScalarDomainType& domain, const Number<P>& value) const = 0;
    /** @brief Construct zero components with compatible characteristics. */
    virtual VectorFunctionPatchInterface<P, SARG>* _create_zeros(SizeType rsize, const ScalarDomainType& domain) const = 0;
    /** @brief Construct constant components with the supplied values. */
    virtual VectorFunctionPatchInterface<P, SARG>* _create_constants(const ScalarDomainType& domain, const Vector<Number<P>>& values) const = 0;
    /** @brief Construct the identity function. */
    virtual ScalarFunctionPatchInterface<P, SARG>* _create_identity(const ScalarDomainType& domain) const = 0;

    /** @brief Construct a zero object with compatible characteristics. */
    virtual ScalarFunctionPatchInterface<P, VARG>* _create_zero(const VectorDomainType& domain) const = 0;
    /** @brief Construct a constant object with the supplied value. */
    virtual ScalarFunctionPatchInterface<P, VARG>* _create_constant(const VectorDomainType& domain, const Number<P>& value) const = 0;
    /** @brief Construct the function returning the selected coordinate. */
    virtual ScalarFunctionPatchInterface<P, VARG>* _create_coordinate(const VectorDomainType& domain, SizeType index) const = 0;
    /** @brief Construct zero components with compatible characteristics. */
    virtual VectorFunctionPatchInterface<P, VARG>* _create_zeros(SizeType rsize, const VectorDomainType& domain) const = 0;
    /** @brief Construct constant components with the supplied values. */
    virtual VectorFunctionPatchInterface<P, VARG>* _create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const = 0;
    /** @brief Construct the function selecting the supplied coordinates. */
    virtual VectorFunctionPatchInterface<P, VARG>* _create_projection(const VectorDomainType& domain, Range indices) const = 0;
    /** @brief Construct the identity function. */
    virtual VectorFunctionPatchInterface<P, VARG>* _create_identity(const VectorDomainType& domain) const = 0;

  public:
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    ScalarFunctionPatch<P, VARG> create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const {
        return ScalarFunctionPatch<P, VARG>(this->_create(domain, function));
    }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    VectorFunctionPatch<P, VARG> create(const VectorDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const {
        return VectorFunctionPatch<P, VARG>(this->_create(domain, function));
    }

    /** @brief Construct a zero object with compatible characteristics. */
    ScalarFunctionPatch<P, VARG> create_zero(const VectorDomainType& domain) const {
        return ScalarFunctionPatch<P, VARG>(_create_zero(domain));
    }
    /** @brief Construct a constant object with the supplied value. */
    ScalarFunctionPatch<P, VARG> create_constant(const VectorDomainType& domain, const Number<P>& value) const {
        return ScalarFunctionPatch<P, VARG>(_create_constant(domain, value));
    }
    /** @brief Construct the function returning the selected coordinate. */
    ScalarFunctionPatch<P, VARG> create_coordinate(const VectorDomainType& domain, SizeType index) const {
        return ScalarFunctionPatch<P, VARG>(_create_coordinate(domain, index));
    }
    /** @brief Construct zero components with compatible characteristics. */
    VectorFunctionPatch<P, VARG> create_zeros(SizeType rsize, const VectorDomainType& domain) const {
        return VectorFunctionPatch<P, VARG>(_create_zeros(rsize, domain));
    }
    /** @brief Construct constant components with the supplied values. */
    VectorFunctionPatch<P, VARG> create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const {
        return VectorFunctionPatch<P, VARG>(_create_constants(domain, values));
    }
    /** @brief Construct the function selecting the supplied coordinates. */
    VectorFunctionPatch<P, VARG> create_projection(const VectorDomainType& domain, Range indices) const {
        return VectorFunctionPatch<P, VARG>(_create_projection(domain, indices));
    }
    /** @brief Construct the identity function. */
    VectorFunctionPatch<P, VARG> create_identity(const VectorDomainType& domain) const {
        return VectorFunctionPatch<P, VARG>(_create_identity(domain));
    }

    /** @brief Construct the identity function. */
    ScalarFunctionPatch<P, VARG> create_identity(const ScalarDomainType& domain) const;
};

} // namespace Ariadne

#endif
