/***************************************************************************
 *            function/function_model_interface.hpp
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

/*! @file function/function_model_interface.hpp
 *  \brief Interface for functions on bounded sets.
 */

#ifndef ARIADNE_FUNCTION_MODEL_INTERFACE_HPP
#define ARIADNE_FUNCTION_MODEL_INTERFACE_HPP

#include "algebra/algebra_interface.hpp"
#include "function/function.decl.hpp"
#include "function/function_interface.hpp"
#include "function/function_patch_interface.hpp"
#include "numeric/operators.hpp"
#include "utility/range.hpp"
#include "utility/writable.hpp"

namespace Ariadne {

template<class P, class PR, class PRE>
class FunctionModelFactoryInterface;
/** @brief Polymorphic interface for function model creator. */
template<class P, class ARG, class PR, class PRE>
class FunctionModelCreatorInterface;

template<class P, class SIG, class PR, class PRE>
class FunctionModelInterface;

/** @brief Polymorphic interface for function model algebra. */
template<class P, class SIG, class PR, class PRE>
class FunctionModelAlgebraInterface;

/** @brief Polymorphic interface for function model algebra. */
template<class P, class ARG, class PR, class PRE>
class FunctionModelAlgebraInterface<P, RealScalar(ARG), PR, PRE>
    : public virtual ElementaryAlgebraInterface<CanonicalNumericType<P, PR, PRE>> {
    static_assert(Same<ARG, RealScalar> or Same<ARG, RealVector>);
    typedef RealScalar RES;
    typedef RES SIG(ARG);
    typedef BoundedDomainType<ARG> D;
    typedef BoundedDomainType<RES> C;
    typedef RawFloat<PR> F;

  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef typename FunctionModelTraits<P, PR, PRE>::ValueType ValueType;

  public:
    /** @brief Return the stored value. */
    virtual const ValueType _value() const = 0;
};

template<class P, class ARG, class PR, class PRE>
using ScalarFunctionModelInterface = FunctionModelInterface<P, RealScalar(ARG), PR, PRE>;
template<class P, class ARG, class PR, class PRE>
using VectorFunctionModelInterface = FunctionModelInterface<P, RealVector(ARG), PR, PRE>;

template<class P, class ARG, class PR, class PRE>
using ScalarFunctionModelAlgebraInterface = FunctionModelAlgebraInterface<P, RealScalar(ARG), PR, PRE>;
template<class P, class ARG, class PR, class PRE>
using VectorFunctionModelAlgebraInterface = FunctionModelAlgebraInterface<P, RealVector(ARG), PR, PRE>;

/** @brief Polymorphic interface for function model algebra. */
template<class P, class ARG, class PR, class PRE>
class FunctionModelAlgebraInterface<P, RealVector(ARG), PR, PRE>
    : public virtual WritableInterface {
  public:
    typedef typename FunctionModelTraits<P, PR, PRE>::ValueType ValueType;
    typedef typename FunctionModelTraits<P, PR, PRE>::ErrorType ErrorType;

  public:
    /** @brief Replace the selected value using the concrete implementation. */
    virtual Void _concrete_set(SizeType, const ScalarFunctionModelInterface<P, ARG, PR, PRE>&) = 0;
    /** @brief Return the requested element or stored value using the concrete implementation. */
    virtual ScalarFunctionModelInterface<P, ARG, PR, PRE>* _concrete_get(SizeType) const = 0;
    /** @brief Join the operands in their natural order using the concrete implementation. */
    virtual VectorFunctionModelInterface<P, ARG, PR, PRE>* _concrete_join(const VectorFunctionModelInterface<P, ARG, PR, PRE>& f2) const = 0;
    /** @brief Combine the operands over their respective arguments using the concrete implementation. */
    virtual VectorFunctionModelInterface<P, ARG, PR, PRE>* _concrete_combine(const VectorFunctionModelInterface<P, ARG, PR, PRE>& f2) const = 0;
    /** @brief Add the supplied value or component using the concrete implementation. */
    virtual Void _concrete_adjoin(const ScalarFunctionModelInterface<P, ARG, PR, PRE>& f2) = 0;

    /** @brief Return the stored values using the concrete implementation. */
    virtual const Vector<ValueType> _concrete_values() const = 0;
    /** @brief Return the component error bounds using the concrete implementation. */
    virtual const Vector<ErrorType> _concrete_errors() const = 0;
};

/** @brief Polymorphic interface for function model. */
template<class P, class SIG, class PR, class PRE>
class FunctionModelInterface
    : public virtual FunctionPatchInterface<P, SIG>,
      public virtual FunctionModelAlgebraInterface<P, SIG, PR, PRE> {
    using RES = typename SignatureTraits<SIG>::ResultKind;
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    static_assert(Same<ARG, RealScalar> or Same<ARG, RealVector>, "");
    static_assert(Same<RES, RealScalar> or Same<RES, RealVector>, "");
    using C = BoundedDomainType<RES>;
    using D = BoundedDomainType<ARG>;
    static_assert(Same<D, IntervalDomainType> or Same<D, BoxDomainType>, "");
    static_assert(Same<C, IntervalDomainType> or Same<C, BoxDomainType>, "");

    template<class X>
    using Argument = typename FunctionInterface<P, SIG>::template Argument<X>;
    template<class X>
    using Result = typename FunctionInterface<P, SIG>::template Result<X>;

    using FLT = Float<PR>;

  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef typename SignatureTraits<SIG>::template ConcreteRangeType<FLT> RangeType;
    typedef typename FunctionModelTraits<P, PR, PRE>::NormType NormType;
    typedef typename FunctionModelTraits<P, PR, PRE>::ValueType ValueType;
    typedef typename FunctionModelTraits<P, PR, PRE>::ErrorType ErrorType;
    typedef typename FunctionModelTraits<P, PR, PRE>::NumericType NumericType;
    typedef typename FunctionModelTraits<P, PR, PRE>::GenericNumericType GenericNumericType;
    typedef typename SignatureTraits<SIG>::ArgumentIndexType ArgumentIndexType;

  public:
    /** @brief Return the domain. */
    virtual const DomainType domain() const override = 0;
    /** @brief Return the codomain. */
    virtual const CodomainType codomain() const override = 0;
    /** @brief Return the range enclosure. */
    virtual const RangeType range() const = 0;
    /** @brief Return the norm using the concrete implementation. */
    virtual const NormType _concrete_norm() const = 0;
    /** @brief Return the error bound using the concrete implementation. */
    virtual const ErrorType _concrete_error() const = 0;

    /** @brief Discard the error bound of the model. */
    virtual Void clobber() = 0;

    /** @brief Evaluate without checking that the argument belongs to the domain using the concrete implementation. */
    virtual Result<CanonicalNumericType<P, PR, PRE>> _concrete_unchecked_evaluate(const Argument<CanonicalNumericType<P, PR, PRE>>& x) const = 0;
    /** @brief Fix the selected argument to the supplied value using the concrete implementation. */
    virtual FunctionModelInterface<P, SIG, PR, PRE>* _concrete_partial_evaluate(SizeType j, const CanonicalNumericType<P, PR, PRE>& c) const = 0;
    /** @brief Embed in the larger argument space. */
    virtual FunctionModelInterface<P, SIG, PR, PRE>* _embed(const DomainType& d1, const DomainType& d2) const override = 0;

    /** @brief Return a factory for compatible objects. */
    virtual FunctionModelFactoryInterface<P, PR, PRE>* _factory() const override = 0;
    /** @brief Return a dynamically allocated copy preserving the concrete type using the concrete implementation. */
    virtual FunctionModelInterface<P, SIG, PR, PRE>* _concrete_clone() const override = 0;
    /** @brief Construct a compatible object from the supplied value and characteristics using the concrete implementation. */
    virtual FunctionModelInterface<P, SIG, PR, PRE>* _concrete_create() const = 0;
    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    inline FunctionModelInterface<P, SIG, PR, PRE>* clone() const { return this->_concrete_clone(); }

    /** @brief Return the object restricted to the supplied domain. */
    virtual FunctionModelInterface<P, SIG, PR, PRE>* _restriction(const DomainType& d) const override = 0;

    /** @brief Return the derivative with respect to the selected argument using the concrete implementation. */
    virtual FunctionModelInterface<P, SIG, PR, PRE>* _concrete_derivative(ArgumentIndexType j) const override = 0;
    /** @brief Return an antiderivative with respect to the selected argument. */
    virtual FunctionModelInterface<P, SIG, PR, PRE>* _antiderivative(ArgumentIndexType j) const override = 0;
    /** @brief Return an antiderivative with respect to the selected argument. */
    virtual FunctionModelInterface<P, SIG, PR, PRE>* _antiderivative(ArgumentIndexType j, Number<P>) const override = 0;
    /** @brief Return an antiderivative with respect to the selected argument using the concrete implementation. */
    virtual FunctionModelInterface<P, SIG, PR, PRE>* _concrete_antiderivative(ArgumentIndexType j, CanonicalNumericType<P, PR, PRE> c) const = 0;

    using FunctionPatchInterface<P, SIG>::_compose;
    using FunctionPatchInterface<P, SIG>::_unchecked_compose;
    /** @brief Compose the outer object with the inner object. */
    virtual ScalarFunctionModelInterface<P, ARG, PR, PRE>* _compose(const ScalarFunction<P, RES>& f) const override = 0;
    /** @brief Compose the outer object with the inner object. */
    virtual VectorFunctionModelInterface<P, ARG, PR, PRE>* _compose(const VectorFunction<P, RES>& f) const override = 0;
    /** @brief Compose without checking domain compatibility. */
    virtual ScalarFunctionModelInterface<P, ARG, PR, PRE>* _unchecked_compose(const ScalarFunction<P, RES>& f) const override = 0;
    /** @brief Compose without checking domain compatibility. */
    virtual VectorFunctionModelInterface<P, ARG, PR, PRE>* _unchecked_compose(const VectorFunction<P, RES>& f) const override = 0;

    /** @brief Test whether the first enclosure is contained in the second using the concrete implementation. */
    virtual Boolean _concrete_refines(const FunctionModelInterface<P, SIG, PR, PRE>& f) const = 0;
    /** @brief Test whether the enclosures cannot represent a common value using the concrete implementation. */
    virtual Boolean _concrete_inconsistent(const FunctionModelInterface<P, SIG, PR, PRE>& f) const = 0;
    /** @brief Return the common refinement of the enclosures using the concrete implementation. */
    virtual FunctionModelInterface<P, SIG, PR, PRE>* _concrete_refinement(const FunctionModelInterface<P, SIG, PR, PRE>& f) const = 0;

    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream&) const override = 0;

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FunctionModelInterface<P, SIG, PR, PRE>& f) {
        return f._write(os);
    }
};

/** @brief Polymorphic interface for function model factory. */
template<class P, class PR, class PRE>
class FunctionModelFactoryInterface
    : public FunctionPatchFactoryInterface<P> {
    typedef RealScalar SARG;
    typedef RealVector VARG;

    typedef BoxDomainType VD;
    typedef IntervalDomainType SD;
    friend class FunctionModelFactory<P, PR, PRE>;

  public:
    typedef SD ScalarDomainType;
    typedef VD VectorDomainType;
    /** @brief Destroy the object and release its owned resources. */
    virtual ~FunctionModelFactoryInterface() override = default;
    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    virtual FunctionModelFactoryInterface<P, PR, PRE>* clone() const override = 0;
    /** @brief Return a copy of the object. */
    inline FunctionModelFactoryInterface<P, PR, PRE>* _copy() const { return this->clone(); }
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const override = 0;
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FunctionModelFactoryInterface<P, PR, PRE>& factory) {
        factory._write(os);
        return os;
    }

  private:
    using FunctionPatchFactoryInterface<P>::_create;
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual CanonicalNumericType<P, PR, PRE> _create(const Number<P>& number) const = 0;
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual ScalarFunctionModelInterface<P, VARG, PR, PRE>* _create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const override = 0;
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual VectorFunctionModelInterface<P, VARG, PR, PRE>* _create(const VectorDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const override = 0;

    /** @brief Construct a zero object with compatible characteristics. */
    virtual ScalarFunctionModelInterface<P, VARG, PR, PRE>* _create_zero(const VectorDomainType& domain) const override = 0;
    /** @brief Construct a constant object with the supplied value. */
    virtual ScalarFunctionModelInterface<P, VARG, PR, PRE>* _create_constant(const VectorDomainType& domain, const Number<P>& value) const override = 0;
    /** @brief Construct the function returning the selected coordinate. */
    virtual ScalarFunctionModelInterface<P, VARG, PR, PRE>* _create_coordinate(const VectorDomainType& domain, SizeType index) const override = 0;
    /** @brief Construct zero components with compatible characteristics. */
    virtual VectorFunctionModelInterface<P, VARG, PR, PRE>* _create_zeros(SizeType rsize, const VectorDomainType& domain) const override = 0;
    /** @brief Construct constant components with the supplied values. */
    virtual VectorFunctionModelInterface<P, VARG, PR, PRE>* _create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const override = 0;
    /** @brief Construct the function selecting the supplied coordinates. */
    virtual VectorFunctionModelInterface<P, VARG, PR, PRE>* _create_projection(const VectorDomainType& domain, Range indices) const override = 0;
    /** @brief Construct the identity function. */
    virtual VectorFunctionModelInterface<P, VARG, PR, PRE>* _create_identity(const VectorDomainType& domain) const override = 0;

  public:
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    CanonicalNumericType<P, PR, PRE> create(const Number<P>& number) const {
        return CanonicalNumericType<P, PR, PRE>(this->_create(number));
    }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    ScalarMultivariateFunctionModel<P, PR, PRE> create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const {
        return ScalarMultivariateFunctionModel<P, PR, PRE>(this->_create(domain, function));
    }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    VectorMultivariateFunctionModel<P, PR, PRE> create(const VectorDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const {
        return VectorMultivariateFunctionModel<P, PR, PRE>(this->_create(domain, function));
    }

    /** @brief Create a scalar number with characteristics matching the function. */
    CanonicalNumericType<P, PR, PRE> create_number(const Number<P>& number) const {
        return CanonicalNumericType<P, PR, PRE>(this->_create(number));
    }

    /** @brief Construct a zero object with compatible characteristics. */
    ScalarMultivariateFunctionModel<P, PR, PRE> create_zero(const VectorDomainType& domain) const {
        return ScalarMultivariateFunctionModel<P, PR, PRE>(_create_zero(domain));
    }
    /** @brief Construct a constant object with the supplied value. */
    ScalarMultivariateFunctionModel<P, PR, PRE> create_constant(const VectorDomainType& domain, const Number<P>& value) const {
        return ScalarMultivariateFunctionModel<P, PR, PRE>(_create_constant(domain, value));
    }
    /** @brief Construct a constant object with the supplied value. */
    ScalarMultivariateFunctionModel<P, PR, PRE> create_constant(const VectorDomainType& domain, const CanonicalNumericType<P, PR, PRE>& value) const {
        return ScalarMultivariateFunctionModel<P, PR, PRE>(_create_constant(domain, Number<P>(value)));
    }
    /** @brief Construct the function returning the selected coordinate. */
    ScalarMultivariateFunctionModel<P, PR, PRE> create_coordinate(const VectorDomainType& domain, SizeType index) const {
        return ScalarMultivariateFunctionModel<P, PR, PRE>(_create_coordinate(domain, index));
    }
    /** @brief Construct zero components with compatible characteristics. */
    VectorMultivariateFunctionModel<P, PR, PRE> create_zeros(SizeType rsize, const VectorDomainType& domain) const {
        return VectorMultivariateFunctionModel<P, PR, PRE>(_create_zeros(rsize, domain));
    }
    /** @brief Construct constant components with the supplied values. */
    VectorMultivariateFunctionModel<P, PR, PRE> create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const {
        return VectorMultivariateFunctionModel<P, PR, PRE>(_create_constants(domain, values));
    }
    /** @brief Construct the function selecting the supplied coordinates. */
    VectorMultivariateFunctionModel<P, PR, PRE> create_projection(const VectorDomainType& domain, Range indices) const {
        return VectorMultivariateFunctionModel<P, PR, PRE>(_create_projection(domain, indices));
    }
    /** @brief Construct the identity function. */
    VectorMultivariateFunctionModel<P, PR, PRE> create_identity(const VectorDomainType& domain) const {
        return VectorMultivariateFunctionModel<P, PR, PRE>(_create_identity(domain));
    }

    /** @brief Construct the identity function. */
    ScalarMultivariateFunctionModel<P, PR, PRE> create_identity(const ScalarDomainType& domain) const;
};

} // namespace Ariadne

#endif
