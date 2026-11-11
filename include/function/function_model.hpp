/***************************************************************************
 *            function/function_model.hpp
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

/*! @file function/function_model.hpp
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_MODEL_HPP
#define ARIADNE_FUNCTION_MODEL_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "algebra/matrix.hpp"
#include "algebra/operations.hpp"
#include "algebra/vector.hpp"
#include "function/domain.hpp"
#include "function/function.decl.hpp"
#include "function/function.hpp"
#include "function/function_interface.hpp"
#include "function/function_mixin.hpp"
#include "function/function_model_interface.hpp"
#include "function_patch.hpp"
#include "numeric/numeric.hpp"
#include "numeric/operators.hpp"

namespace Ariadne {

template<class P, class A, class PR, class PRE>
struct AlgebraOperations<ScalarFunctionModel<P, Real(A), PR, PRE>>;

// FIXME: Extend with univariate case
/** @brief Factory constructing compatible function model. */
template<class P, class PR, class PRE>
class FunctionModelFactory
    : public Handle<const FunctionModelFactoryInterface<P, PR, PRE>> {
  private:
    typedef RealScalar SARG;
    typedef RealVector VARG;
    typedef IntervalDomainType SD;
    typedef BoxDomainType VD;

  public:
    typedef FunctionModelFactoryInterface<P, PR, PRE> Interface;

    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;
    typedef SD ScalarDomainType;
    typedef VD VectorDomainType;

    using Handle<const Interface>::Handle;

    /** @brief Construct a compatible object from the supplied value and characteristics. */
    CanonicalNumericType<P, PR, PRE> create(const Number<P>& c) const {
        return CanonicalNumericType<P, PR, PRE>(this->_ptr->_create(c));
    }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    ScalarMultivariateFunctionModel<P, PR, PRE> create(const VectorDomainType& dom, const ScalarMultivariateFunction<P>& f) const {
        return ScalarMultivariateFunctionModel<P, PR, PRE>(this->_ptr->_create(dom, f));
    }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    VectorMultivariateFunctionModel<P, PR, PRE> create(const VectorDomainType& dom, const VectorMultivariateFunction<P>& f) const {
        return VectorMultivariateFunctionModel<P, PR, PRE>(this->_ptr->_create(dom, f));
    }

    /** @brief Create a scalar number with characteristics matching the function. */
    CanonicalNumericType<P, PR, PRE> create_number(const Number<P>& c) const {
        return CanonicalNumericType<P, PR, PRE>(this->_ptr->_create(c));
    }
    /** @brief Construct a zero object with compatible characteristics. */
    ScalarMultivariateFunctionModel<P, PR, PRE> create_zero(const VectorDomainType& dom) const {
        return ScalarMultivariateFunctionModel<P, PR, PRE>(this->_ptr->_create_zero(dom));
    }
    /** @brief Construct a constant object with the supplied value. */
    ScalarMultivariateFunctionModel<P, PR, PRE> create_constant(const VectorDomainType& dom, const Number<P>& c) const {
        return ScalarMultivariateFunctionModel<P, PR, PRE>(this->_ptr->_create_constant(dom, c));
    }
    /** @brief Construct a constant object with the supplied value. */
    ScalarMultivariateFunctionModel<P, PR, PRE> create_constant(const VectorDomainType& dom, const CanonicalNumericType<P, PR, PRE>& c) const {
        return ScalarMultivariateFunctionModel<P, PR, PRE>(this->_ptr->_create_constant(dom, c));
    }
    /** @brief Construct the function returning the selected coordinate. */
    ScalarMultivariateFunctionModel<P, PR, PRE> create_coordinate(const VectorDomainType& dom, SizeType index) const {
        return ScalarMultivariateFunctionModel<P, PR, PRE>(this->_ptr->_create_coordinate(dom, index));
    }
    /** @brief Construct zero components with compatible characteristics. */
    VectorMultivariateFunctionModel<P, PR, PRE> create_zeros(SizeType rsize, const VectorDomainType& dom) const {
        return VectorMultivariateFunctionModel<P, PR, PRE>(this->_ptr->_create_zeros(rsize, dom));
    }
    /** @brief Construct constant components with the supplied values. */
    VectorMultivariateFunctionModel<P, PR, PRE> create_constants(const VectorDomainType& dom, const Vector<Number<P>>& c) const {
        return VectorMultivariateFunctionModel<P, PR, PRE>(this->_ptr->_create_constants(dom, c));
    }
    /** @brief Construct the function selecting the supplied coordinates. */
    VectorMultivariateFunctionModel<P, PR, PRE> create_projection(const VectorDomainType& dom, Range indices) const {
        return VectorMultivariateFunctionModel<P, PR, PRE>(this->_ptr->_create_projection(dom, indices));
    }
    /** @brief Construct the identity function. */
    VectorMultivariateFunctionModel<P, PR, PRE> create_identity(const VectorDomainType& dom) const {
        return VectorMultivariateFunctionModel<P, PR, PRE>(this->_ptr->_create_identity(dom));
    }

    // FIXME: Should return a univariate model
    /** @brief Construct the identity function. */
    ScalarFunctionModel<P, VARG, PR, PRE> create_identity(const ScalarDomainType& dom) const {
        return ScalarFunctionModel<P, VARG, PR, PRE>(this->_ptr->_create_coordinate(VectorDomainType(1u, dom), 0u));
    }
    /** @brief Construct zero components with compatible characteristics. */
    VectorFunctionModel<P, SARG, PR, PRE> create_zeros(SizeType, const ScalarDomainType&) const {
        ARIADNE_NOT_IMPLEMENTED;
    }

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FunctionModelFactory<P, PR, PRE>& factory) { return factory._ptr->_write(os); }
};

template<class P, class PR, class PRE>
auto inline FunctionModelFactoryInterface<P, PR, PRE>::create_identity(const IntervalDomainType& dom) const -> ScalarFunctionModel<P, VARG, PR, PRE> {
    return ScalarFunctionModel<P, VARG, PR, PRE>(this->create_coordinate(BoxDomainType(1u, dom), 0u));
}

/** @brief Construction helper for function model. */
template<class FCTRY, class ARG>
class FunctionModelCreator {
    typedef BoundedDomainType<ARG> D;
    typedef typename FCTRY::Paradigm P;
    typedef typename FCTRY::PrecisionType PR;
    typedef typename FCTRY::ErrorPrecisionType PRE;

  public:
    typedef FCTRY FactoryType;
    typedef D DomainType;
    typedef P Paradigm;

    /** @brief Construct a \c FunctionModelCreator from domain, factory. */
    explicit FunctionModelCreator(DomainType domain, FactoryType factory) : _factory(factory), _domain(domain) {}

    /** @brief Construct a compatible object from the supplied value and characteristics. */
    decltype(auto) create(const Number<P>& c) const
        requires requires(const FactoryType& factory) { factory.create(c); }
    { return this->_factory.create(c); }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    decltype(auto) create(const ScalarFunction<P, ARG>& f)
        requires requires(FactoryType& factory, DomainType& domain) { factory.create(domain, f); }
    { return this->_factory.create(this->_domain, f); }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    decltype(auto) create(const VectorFunction<P, ARG>& f)
        requires requires(FactoryType& factory, DomainType& domain) { factory.create(domain, f); }
    { return this->_factory.create(this->_domain, f); }
    /** @brief Construct a zero object with compatible characteristics. */
    decltype(auto) create_zero()
        requires requires(FactoryType& factory, DomainType& domain) { factory.create_zero(domain); }
    { return this->_factory.create_zero(this->_domain); }
    /** @brief Construct zero components with compatible characteristics. */
    decltype(auto) create_zeros(SizeType n)
        requires requires(FactoryType& factory, DomainType& domain) { factory.create_zeros(n, domain); }
    { return this->_factory.create_zeros(n, this->_domain); }
    /** @brief Construct a constant object with the supplied value. */
    decltype(auto) create_constant(const Number<P>& c) const
        requires requires(const FactoryType& factory, const DomainType& domain) { factory.create_constant(domain, c); }
    { return this->_factory.create_constant(this->_domain, c); }
    /** @brief Construct constant components with the supplied values. */
    decltype(auto) create_constants(const Vector<Number<P>>& c) const
        requires requires(const FactoryType& factory, const DomainType& domain) { factory.create_constants(domain, c); }
    { return this->_factory.create_constants(this->_domain, c); }
    /** @brief Construct the identity function. */
    decltype(auto) create_identity()
        requires requires(FactoryType& factory, DomainType& domain) { factory.create_identity(domain); }
    { return this->_factory.create_identity(this->_domain); }

    /** @brief Construct a compatible object from the supplied value and characteristics. */
    decltype(auto) create(const DomainType& dom, const ScalarFunction<P, ARG>& f)
        requires requires(FactoryType& factory) { factory.create(dom, f); }
    { return this->_factory.create(dom, f); }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    decltype(auto) create(const DomainType& dom, const VectorFunction<P, ARG>& f)
        requires requires(FactoryType& factory) { factory.create(dom, f); }
    { return this->_factory.create(dom, f); }
    /** @brief Construct a zero object with compatible characteristics. */
    decltype(auto) create_zero(const DomainType& dom)
        requires requires(FactoryType& factory) { factory.create_zero(dom); }
    { return this->_factory.create_zero(dom); }
    /** @brief Construct zero components with compatible characteristics. */
    decltype(auto) create_zeros(const DomainType& dom, SizeType n)
        requires requires(FactoryType& factory) { factory.create_zeros(n, dom); }
    { return this->_factory.create_zeros(n, dom); }
    /** @brief Construct a constant object with the supplied value. */
    decltype(auto) create_constant(const DomainType& dom, const Number<P>& c) const
        requires requires(const FactoryType& factory) { factory.create_constant(dom, c); }
    { return this->_factory.create_constant(dom, c); }
    /** @brief Construct constant components with the supplied values. */
    decltype(auto) create_constants(const DomainType& dom, const Vector<Number<P>>& c) const
        requires requires(const FactoryType& factory) { factory.create_constants(dom, c); }
    { return this->_factory.create_constants(dom, c); }
    /** @brief Construct the identity function. */
    decltype(auto) create_identity(const DomainType& dom)
        requires requires(FactoryType& factory) { factory.create_identity(dom); }
    { return this->_factory.create_identity(dom); }

    /** @brief Construct a compatible object from the supplied value and characteristics. */
    const CanonicalNumericType<P, PR, PRE>& create(const CanonicalNumericType<P, PR, PRE>& c) const { return c; }

  protected:
    FactoryType _factory;
    DomainType _domain;
};

// FIXME: Merge with multivariate case
/** @brief Construction helper for function model. */
template<class FCTRY>
class FunctionModelCreator<FCTRY, RealScalar> {
    typedef typename FCTRY::Paradigm P;
    typedef typename FCTRY::PrecisionType PR;
    typedef typename FCTRY::ErrorPrecisionType PRE;
    typedef IntervalDomainType D;
    typedef RealScalar ARG;

  public:
    typedef FCTRY FactoryType;
    typedef D DomainType;
    typedef P Paradigm;

    /** @brief Construct a \c FunctionModelCreator from domain, factory. */
    explicit FunctionModelCreator(DomainType domain, FactoryType factory) : _factory(factory), _domain(domain) {}

    /** @brief Construct a compatible object from the supplied value and characteristics. */
    CanonicalNumericType<P, PR, PRE> create(const Number<P>& c) const { return this->_factory.create(c); }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    ScalarFunctionModel<P, ARG, PR, PRE> create(const ScalarFunction<P, ARG>&) { ARIADNE_NOT_IMPLEMENTED; }
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    VectorFunctionModel<P, ARG, PR, PRE> create(const VectorFunction<P, ARG>&) { ARIADNE_NOT_IMPLEMENTED; }
    /** @brief Construct a zero object with compatible characteristics. */
    ScalarFunctionModel<P, ARG, PR, PRE> create_zero() { ARIADNE_NOT_IMPLEMENTED; }
    /** @brief Construct zero components with compatible characteristics. */
    VectorFunctionModel<P, ARG, PR, PRE> create_zeros(SizeType) { ARIADNE_NOT_IMPLEMENTED; }
    /** @brief Construct the identity function. */
    ScalarFunctionModel<P, ARG, PR, PRE> create_identity() { ARIADNE_NOT_IMPLEMENTED; }

    /** @brief Construct a compatible object from the supplied value and characteristics. */
    const CanonicalNumericType<P, PR, PRE>& create(const CanonicalNumericType<P, PR, PRE>& c) const { return c; }

  protected:
    FactoryType _factory;
    DomainType _domain;
};

//! \ingroup FunctionModule
//! \ingroup FunctionModelSubModule
//! \brief Generic class representing approximations to functions on bounded domains.
//!  \tparam P	The information paradigm tag, which can be either ValidatedTag, indicating that the approximation has a known (uniform) error bound, or ApproximateTag, indicating that no error bound is available. See the \ref InformationSubModule for more details.
//!  \tparam SIG The signature, which has the standard C++ form RES(ARG), so signature Real(RealVector) indicates a function \f$f:\R^n\to\R\f$. See \ref function_signature_section for more details.
//!  \tparam PR The precision used for numerical values within the approximation.
//!  \tparam PRE The precision used for error bounds provided for the approximation.
//! \paragraph function_model_template_parameter_note Note:
//!  It is planned to change the precision paramters \c PR and \c PRE to the actual number types used.
template<class P, class SIG, class PR, class PRE>
class FunctionModel;

//! \ingroup FunctionModelSubModule
//! \brief Generic scalar functions on bounded domains.
template<class P, class ARG, class PR, class PRE>
class FunctionModel<P, RealScalar(ARG), PR, PRE>
    : public Handle<FunctionModelInterface<P, RealScalar(ARG), PR, PRE>>
    //    , public DispatchTranscendentalAlgebraOperations<ScalarFunctionModel<P,D,PR,PRE>, CanonicalNumericType<P,PR,PRE>>
    ,
      public DispatchElementaryAlgebraOperations<ScalarFunctionModel<P, ARG, PR, PRE>, CanonicalNumericType<P, PR, PRE>>,
      public ProvideConcreteGenericElementaryOperations<ScalarFunctionModel<P, ARG, PR, PRE>, ScalarMultivariateFunction<P>>,
      public ProvideConcreteGenericElementaryOperations<ScalarFunctionModel<P, ARG, PR, PRE>, Number<P>> {
    static_assert(Same<ARG, RealScalar> or Same<ARG, RealVector>, "");
    using RES = RealScalar;
    using SIG = RES(ARG);
    using C = typename SignatureTraits<SIG>::BoundedCodomainType;
    using D = typename SignatureTraits<SIG>::BoundedDomainType;
    static_assert(Same<D, IntervalDomainType> or Same<D, BoxDomainType>, "");

  public:
    typedef FunctionModelInterface<P, SIG, PR, PRE> Interface;
    typedef ScalarFunction<P, ARG> GenericType;
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef D DomainType;
    typedef C CodomainType;
    typedef typename Interface::ValueType ValueType;
    typedef typename Interface::ErrorType ErrorType;
    typedef typename Interface::NumericType NumericType;
    typedef typename Interface::GenericNumericType GenericNumericType;
    typedef typename Interface::NormType NormType;
    typedef typename Interface::RangeType RangeType;
    typedef typename Interface::ArgumentIndexType ArgumentIndexType;

    template<class Y>
    using Argument = typename SignatureTraits<SIG>::template Argument<Y>;
    template<class Y>
    using Result = typename SignatureTraits<SIG>::template Result<Y>;

  public:
    /** @brief Construct a \c FunctionModel from p. */
    explicit FunctionModel(Interface* p) : Handle<Interface>(p) {}
    /** @brief Construct a \c FunctionModel from p. */
    FunctionModel(UniquePointer<Interface> p) : Handle<Interface>(std::move(p)) {}

    /** @brief Construct a \c FunctionModel with the specified characteristics. */
    FunctionModel() : FunctionModel(nullptr) {}
    /** @brief Construct a \c FunctionModel from p. */
    FunctionModel(const SharedPointer<const Interface> p) : Handle<Interface>(p->_concrete_clone()) {}
    /** @brief Copy-construct the object from the supplied instance. */
    FunctionModel(const FunctionModel<P, SIG, PR, PRE>& f) : Handle<Interface>(f._ptr->_concrete_clone()) {}
    /** @brief Assign the supplied value and return this object. */
    FunctionModel& operator=(const FunctionModel<P, SIG, PR, PRE>& f) {
        this->Handle<FunctionModelInterface<P, RealScalar(ARG), PR, PRE>>::operator=(f);
        return *this;
    }
    /** @brief Construct a \c FunctionModel from f. */
    FunctionModel(const Interface& f) : Handle<Interface>(f._concrete_clone()) {}
    /** @brief Construct a \c FunctionModel from f. */
    FunctionModel(const Function<P, SIG>& f) : FunctionModel(dynamic_cast<Interface*>(f.raw_pointer()->_clone())) {}
    /** @brief Convert to function< p, sig >. */
    operator Function<P, SIG>() const { return Function<P, SIG>(this->_ptr->_clone()); }
    /** @brief Convert to interface &. */
    operator Interface&() { return this->reference(); }
    /** @brief Convert to const interface &. */
    operator const Interface&() const { return this->reference(); }

    /** @brief Assign the supplied value and return this object. */
    ScalarFunctionModel<P, ARG, PR, PRE>& operator=(const Number<P>& c);
    /** @brief Assign the supplied value and return this object. */
    ScalarFunctionModel<P, ARG, PR, PRE>& operator=(const CanonicalNumericType<P, PR, PRE>& c);
    /** @brief Assign the supplied value and return this object. */
    ScalarFunctionModel<P, ARG, PR, PRE>& operator=(const ScalarFunction<P, ARG>& f);
    /** @brief Assign the supplied value and return this object. */
    ScalarFunctionModel<P, ARG, PR, PRE>& operator=(const ScalarFunctionModelInterface<P, ARG, PR, PRE>& f);
    //    ScalarFunctionModel<P,ARG,PR,PRE>& operator=(const ValidatedScalarMultivariateTaylorFunctionModelDP& f);

    /** @brief Return the working precision. */
    inline const PrecisionType precision() const { return this->value().precision(); }
    /** @brief Return the number of scalar arguments. */
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    X operator()(const Argument<X>& x) const {
        return this->_ptr->_call(x);
    }
    /** @brief Evaluate at the supplied argument. */
    template<class X>
    X evaluate(const Argument<X>& x) const {
        return this->_ptr->_call(x);
    }
    /** @brief Return the domain. */
    inline const DomainType domain() const { return this->_ptr->domain(); }
    /** @brief Return the codomain. */
    inline const CodomainType codomain() const { return this->_ptr->codomain(); }
    /** @brief Return the range enclosure. */
    inline const RangeType range() const { return this->_ptr->range(); }

    /** @brief Return the stored value. */
    inline ValueType value() const { return this->_ptr->_value(); }
    /** @brief Return the error bound. */
    inline ErrorType error() const { return this->_ptr->_concrete_error(); }

    /** @brief Discard the error bound of the model. */
    inline Void clobber() { this->pointer()->clobber(); }

    //    inline ScalarFunctionModel<P,ARG,PR,PRE> apply(UnaryElementaryOperator op) const { return ScalarFunctionModel<P,ARG,PR,PRE>(this->_ptr->_apply(op)); }
    /** @brief Restrict to the supplied subset or domain. */
    inline Void restrict(const DomainType& d) { *this = restriction(*this, d); }

  public:
    /** @brief Return a factory for compatible objects. */
    friend FunctionModelCreator<FunctionModelFactory<P, PR, PRE>, ARG> factory(const ScalarFunctionModel<P, ARG, PR, PRE>& f) {
        FunctionModelFactory<P, PR, PRE> factory(f._ptr->_factory());
        return FunctionModelCreator<FunctionModelFactory<P, PR, PRE>, ARG>(f.domain(), factory);
    }

  public:

  public:
    /** @brief Evaluate at the supplied argument. */
    friend CanonicalNumericType<P, PR, PRE> evaluate(const ScalarFunctionModel<P, ARG, PR, PRE>& f, const Argument<CanonicalNumericType<P, PR, PRE>>& x) {
        return f._ptr->_call(x);
    }
    /** @brief Evaluate at the supplied argument. */
    friend Number<P> evaluate(const ScalarFunctionModel<P, ARG, PR, PRE>& f, const Argument<Number<P>>& x) {
        return f._ptr->_call(Argument<CanonicalNumericType<P, PR, PRE>>(x, f.precision()));
    }
    /** @brief Evaluate without checking that the argument belongs to the domain. */
    friend CanonicalNumericType<P, PR, PRE> unchecked_evaluate(const ScalarFunctionModel<P, ARG, PR, PRE>& f, const Argument<CanonicalNumericType<P, PR, PRE>>& x) {
        return f._ptr->_unchecked_evaluate(x);
    }
    /** @brief Evaluate without checking that the argument belongs to the domain. */
    friend Number<P> unchecked_evaluate(const ScalarFunctionModel<P, ARG, PR, PRE>& f, const Argument<Number<P>>& x) {
        return f._ptr->_unchecked_evaluate(Argument<CanonicalNumericType<P, PR, PRE>>(x, f.precision()));
    }

    /** @brief Fix the selected argument to the supplied value. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> partial_evaluate(const ScalarFunctionModel<P, ARG, PR, PRE>& f, SizeType j, const CanonicalNumericType<P, PR, PRE>& c) {
        return ScalarFunctionModel<P, ARG, PR, PRE>(f._ptr->_concrete_partial_evaluate(j, c));
    }
    /** @brief Fix the selected argument to the supplied value. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> partial_evaluate(const ScalarFunctionModel<P, ARG, PR, PRE>& f, SizeType j, const Number<P>& c) {
        return partial_evaluate(f, j, CanonicalNumericType<P, PR, PRE>(c, f.precision()));
    }

    /** @brief Return the norm. */
    friend NormType norm(const ScalarFunctionModel<P, ARG, PR, PRE>& f) {
        return f._ptr->_concrete_norm();
    }
    /** @brief Return the derivative with respect to the selected argument. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> derivative(const ScalarFunctionModel<P, ARG, PR, PRE>& f, ArgumentIndexType j) {
        return ScalarFunctionModel<P, ARG, PR, PRE>(f._ptr->_concrete_derivative(j));
    }
    /** @brief Return an antiderivative with respect to the selected argument. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> antiderivative(const ScalarFunctionModel<P, ARG, PR, PRE>& f, ArgumentIndexType j) {
        return ScalarFunctionModel<P, ARG, PR, PRE>(f._ptr->_antiderivative(j));
    }
    /** @brief Return an antiderivative with respect to the selected argument. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> antiderivative(const ScalarFunctionModel<P, ARG, PR, PRE>& f, ArgumentIndexType j, CanonicalNumericType<P, PR, PRE> c) {
        return ScalarFunctionModel<P, ARG, PR, PRE>(f._ptr->_antiderivative(j, c));
    }
    /** @brief Return an antiderivative with respect to the selected argument. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> antiderivative(const ScalarFunctionModel<P, ARG, PR, PRE>& f, ArgumentIndexType j, const Number<P>& c) {
        return antiderivative(f, j, CanonicalNumericType<P, PR, PRE>(c, f.value().precision()));
    }

    /** @brief Embed in the larger argument space. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> embed(const DomainType& d1, const ScalarFunctionModel<P, ARG, PR, PRE>& f, const DomainType& d2)
        requires Same<ARG, RealVector>
    {
        return ScalarFunctionModel<P, ARG, PR, PRE>(f._ptr->_embed(d1, d2));
    }
    /** @brief Embed in the larger argument space. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> embed(const DomainType& d, const ScalarFunctionModel<P, ARG, PR, PRE>& f)
        requires Same<ARG, RealVector>
    {
        return embed(d, f, DomainType());
    }
    /** @brief Embed in the larger argument space. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> embed(const ScalarFunctionModel<P, ARG, PR, PRE>& f, const BoxDomainType& d)
        requires Same<ARG, RealVector>
    {
        return embed(DomainType(), f, d);
    }
    /** @brief Embed in the larger argument space. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> embed(const ScalarFunctionModel<P, ARG, PR, PRE>& f, const IntervalDomainType& d)
        requires Same<ARG, RealVector>
    {
        return embed(f, DomainType(1, d));
    }
    /** @brief Restrict to the supplied subset or domain. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> restrict(const ScalarFunctionModel<P, ARG, PR, PRE>& f, const DomainType& d) {
        return ScalarFunctionModel<P, ARG, PR, PRE>(f._ptr->_restriction(d));
    }
    /** @brief Return the object restricted to the supplied domain. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> restriction(const ScalarFunctionModel<P, ARG, PR, PRE>& f, const DomainType& d) {
        return ScalarFunctionModel<P, ARG, PR, PRE>(f._ptr->_restriction(d));
    }

    /** @brief Join the operands in their natural order. */
    friend VectorFunctionModel<P, ARG, PR, PRE> join(const ScalarFunctionModel<P, ARG, PR, PRE>& f1, const ScalarFunctionModel<P, ARG, PR, PRE>& f2) {
        return join(VectorFunctionModel<P, ARG, PR, PRE>(1, f1), f2);
    }
    /** @brief Combine the operands over their respective arguments. */
    friend VectorFunctionModel<P, ARG, PR, PRE> combine(const ScalarFunctionModel<P, ARG, PR, PRE>& f1, const ScalarFunctionModel<P, ARG, PR, PRE>& f2);

  public:
    typedef ValidatedTag VP;
    /** @brief Return the common refinement of the enclosures. */
    friend ScalarFunctionModel<VP, ARG, PR, PRE> refinement(const ScalarFunctionModel<VP, ARG, PR, PRE>& f1, const ScalarFunctionModel<VP, ARG, PR, PRE>& f2) {
        return ScalarFunctionModel<VP, ARG, PR, PRE>(f1._ptr->_concrete_refinement(f2));
    }
    /** @brief Test whether the enclosures cannot represent a common value. */
    friend Boolean inconsistent(const ScalarFunctionModel<VP, ARG, PR, PRE>& f1, const ScalarFunctionModel<VP, ARG, PR, PRE>& f2) {
        return f1._ptr->_concrete_inconsistent(f2);
    }
    /** @brief Test whether the first enclosure is contained in the second. */
    friend Boolean refines(const ScalarFunctionModel<VP, ARG, PR, PRE>& f1, const ScalarFunctionModel<VP, ARG, PR, PRE>& f2) {
        return f1._ptr->_concrete_refines(f2);
    }

  public:
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const ScalarFunctionModel<P, ARG, PR, PRE>& f) {
        return os << static_cast<GenericType>(f);
    }
};

/** @brief Operations supported by algebra objects. */
template<class P, class ARG, class PR, class PRE>
struct AlgebraOperations<ScalarFunctionModel<P, ARG, PR, PRE>> {
    /*
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Nul, ScalarFunctionModel<P,ARG,PR,PRE> f) {
            f._ptr->_imul(CanonicalNumericType<P,PR,PRE>(0)); return f; }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Neg, ScalarFunctionModel<P,ARG,PR,PRE> f) {
            f._ptr->_imul(CanonicalNumericType<P,PR,PRE>(-1)); return f; }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Add, ScalarFunctionModel<P,ARG,PR,PRE> f1, const ScalarFunctionModel<P,ARG,PR,PRE>& f2) {
            f1._ptr->_isma(CanonicalNumericType<P,PR,PRE>(+1),f2); return f1; }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Sub, ScalarFunctionModel<P,ARG,PR,PRE> f1, const ScalarFunctionModel<P,ARG,PR,PRE>& f2) {
            f1._ptr->_isma(CanonicalNumericType<P,PR,PRE>(-1),f2); return f1; }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Mul, const ScalarFunctionModel<P,ARG,PR,PRE>& f1, const ScalarFunctionModel<P,ARG,PR,PRE>& f2) {
            ScalarFunctionModel<P,ARG,PR,PRE> r=factory(f1).create_zero(); r._ptr->_ifma(f1,f2); return r; }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Add, ScalarFunctionModel<P,ARG,PR,PRE> f1, const CanonicalNumericType<P,PR,PRE>& c2) {
            f1._ptr->_iadd(c2); return f1; }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Mul, ScalarFunctionModel<P,ARG,PR,PRE> f1, const CanonicalNumericType<P,PR,PRE>& c2) {
            f1._ptr->_imul(c2); return f1; }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Add, const CanonicalNumericType<P,PR,PRE>& c1, ScalarFunctionModel<P,ARG,PR,PRE> f2) {
            f2._ptr->_iadd(c1); return f2; }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Mul, const CanonicalNumericType<P,PR,PRE>& c1, ScalarFunctionModel<P,ARG,PR,PRE> f2) {
            f2._ptr->_imul(c1); return f2; }

        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Sub, ScalarFunctionModel<P,ARG,PR,PRE> f1, const CanonicalNumericType<P,PR,PRE>& c2) {
            return add(std::move(f1),neg(c2)); }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Div, ScalarFunctionModel<P,ARG,PR,PRE> f1, const CanonicalNumericType<P,PR,PRE>& c2) {
            return mul(std::move(f1),rec(c2)); }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Sub, const CanonicalNumericType<P,PR,PRE>& c1, ScalarFunctionModel<P,ARG,PR,PRE> f2) {
            return add(neg(std::move(f2)),c1); }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Div, const CanonicalNumericType<P,PR,PRE>& c1, ScalarFunctionModel<P,ARG,PR,PRE> f2) {
            return mul(rec(std::move(f2)),c1); }

        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Rec, const ScalarFunctionModel<P,ARG,PR,PRE>& f) {
            return f.apply(Rec()); }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Div, const ScalarFunctionModel<P,ARG,PR,PRE>& f1, const ScalarFunctionModel<P,ARG,PR,PRE>& f2) {
            return mul(f1,rec(f2)); }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Pow, const ScalarFunctionModel<P,ARG,PR,PRE>& f1, Int n2) {
            return generic_pow(f1,n2); }

        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Abs, const ScalarFunctionModel<P,ARG,PR,PRE>& f) {
            return f.apply(Abs()); }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Max, const ScalarFunctionModel<P,ARG,PR,PRE>& f1, const ScalarFunctionModel<P,ARG,PR,PRE>& f2) {
            return hlf(add(add(f1,f2),abs(sub(f1,f2)))); }
        static ScalarFunctionModel<P,ARG,PR,PRE> apply(Min, const ScalarFunctionModel<P,ARG,PR,PRE>& f1, const ScalarFunctionModel<P,ARG,PR,PRE>& f2) {
            return hlf(sub(add(f1,f2),abs(sub(f1,f2)))); }
    */

    typedef ScalarFunctionModel<P, ARG, PR, PRE> FM;
    typedef CanonicalNumericType<P, PR, PRE> X;
    typedef ScalarFunctionModelInterface<P, ARG, PR, PRE> FMI;
    typedef ScalarFunctionModelAlgebraInterface<P, ARG, PR, PRE> FMAI;

    /** @brief Apply the supplied operation to the operands. */
    static FM apply(BinaryElementaryOperator op, const FM& f1, const FM& f2) {
        return FM(&dynamic_cast<FMI&>(*static_cast<FMAI const&>(f1.reference())._apply(op, f2.reference())));
    }
    /** @brief Apply the supplied operation to the operands. */
    static FM apply(BinaryElementaryOperator op, const FM& f1, const X& c2) {
        return FM(&dynamic_cast<FMI&>(*static_cast<FMAI const&>(f1.reference())._apply(op, c2)));
    }
    /** @brief Apply the supplied operation to the operands. */
    static FM apply(BinaryElementaryOperator op, const X& c1, const FM& f2) {
        return FM(&dynamic_cast<FMI&>(*static_cast<FMAI const&>(f2.reference())._rapply(op, c1)));
    }
    /** @brief Apply the supplied operation to the operands. */
    static FM apply(UnaryElementaryOperator op, const FM& f) {
        return FM(&dynamic_cast<FMI&>(*static_cast<FMAI const&>(f.reference())._apply(op)));
    }
    /** @brief Apply the supplied operation to the operands. */
    static FM apply(GradedElementaryOperator op, const FM& f, Int n) {
        return FM(&dynamic_cast<FMI&>(*static_cast<FMAI const&>(f.reference())._apply(op, n)));
    }

    /*
        static FM apply(BinaryElementaryOperator op, const X& c1, const FM& f2) {
            FM s1=factory(f2).create_constant(c1); return op(s1,f2); }

        static FM apply(BinaryElementaryOperator op, const FM& f1, const Number<P>& c2) {
            return op(f1,f1.create_constant(c2)); }
        static FM apply(BinaryElementaryOperator, const Number<P>& c1, const FM& f2) {
            return op(f2.create_constant(c1),f2); }

    */
};

template<class P, class ARG, class PR, class PRE>
inline auto
    FunctionModel<P, RealScalar(ARG), PR, PRE>::operator=(const CanonicalNumericType<P, PR, PRE>& c) -> ScalarFunctionModel<P, ARG, PR, PRE>& {
    (*this) *= nul(c);
    (*this) += c;
    return *this;
}
template<class P, class ARG, class PR, class PRE>
inline auto
    FunctionModel<P, RealScalar(ARG), PR, PRE>::operator=(const Number<P>& c) -> ScalarFunctionModel<P, ARG, PR, PRE>& {
    return (*this) = factory(*this).create(c);
}
template<class P, class ARG, class PR, class PRE>
inline auto
    FunctionModel<P, RealScalar(ARG), PR, PRE>::operator=(const ScalarFunction<P, ARG>& f) -> ScalarFunctionModel<P, ARG, PR, PRE>& {
    return (*this) = factory(*this).create(f);
}
template<class P, class ARG, class PR, class PRE>
inline auto
    FunctionModel<P, RealScalar(ARG), PR, PRE>::operator=(const ScalarFunctionModelInterface<P, ARG, PR, PRE>& f) -> ScalarFunctionModel<P, ARG, PR, PRE>& {
    return (*this) = ScalarFunctionModel<P, ARG, PR, PRE>(f._concrete_clone());
}

/** @brief Proxy referring to a selected component of a composite object. */
template<class V>
struct Element;

template<class M>
class ScaledFunctionPatch;
template<class M>
class VectorScaledFunctionPatch;
/** @brief Proxy referring to a selected component of a composite object. */
template<class M>
struct Element<VectorScaledFunctionPatch<M>> {
    typedef ScaledFunctionPatch<M> Type;
};

typedef ScaledFunctionPatch<ValidatedTaylorModelDP> ValidatedScalarMultivariateTaylorFunctionModelDP;

/** @brief Scalar component view of a vector function model. */
template<class P, class ARG, class PR, class PRE>
class VectorFunctionModelElement
    : public DispatchTranscendentalAlgebraOperations<ScalarFunctionModel<P, ARG, PR, PRE>, CanonicalNumericType<P, PR, PRE>>,
      public ProvideConcreteGenericArithmeticOperators<ScalarFunctionModel<P, ARG, PR, PRE>> {
    VectorFunctionModel<P, ARG, PR, PRE>* _p;
    SizeType _i;

  public:
    typedef typename ScalarFunctionModel<P, ARG, PR, PRE>::GenericType GenericType;
    typedef typename ScalarFunctionModel<P, ARG, PR, PRE>::RangeType RangeType;
    /** @brief Convert to const scalar function model< p, arg, pr, pre >. */
    operator const ScalarFunctionModel<P, ARG, PR, PRE>() const;
    /** @brief Construct a \c VectorFunctionModelElement from p, i. */
    VectorFunctionModelElement(VectorFunctionModel<P, ARG, PR, PRE>* p, SizeType i) : _p(p), _i(i) {}
    /** @brief Assign the supplied value and return this object. */
    VectorFunctionModelElement<P, ARG, PR, PRE>& operator=(const ScalarFunctionModel<P, ARG, PR, PRE>& sf) {
        _p->set(_i, sf);
        return *this;
    }
    /** @brief Assign the supplied value and return this object. */
    VectorFunctionModelElement<P, ARG, PR, PRE>& operator=(const VectorFunctionModelElement<P, ARG, PR, PRE>& sf) {
        return this->operator=(static_cast<const ScalarFunctionModel<P, ARG, PR, PRE>>(sf));
    }
    /** @brief Discard the error bound of the model. */
    Void clobber() {
        ScalarFunctionModel<P, ARG, PR, PRE> sf = _p->get(_i);
        sf.clobber();
        _p->set(_i, sf);
    }
    /** @brief Return the selected function model. */
    const ScalarFunctionModel<P, ARG, PR, PRE> model() const {
        ScalarFunctionModel<P, ARG, PR, PRE> sf = _p->get(_i);
        return sf.model();
    }
    /** @brief Return the error bound. */
    const CanonicalErrorType<P, PRE> error() const {
        ScalarFunctionModel<P, ARG, PR, PRE> sf = _p->get(_i);
        return sf.error();
    }
    /** @brief Set the error. */
    Void set_error(CanonicalErrorType<P, PRE> e) const {
        ScalarFunctionModel<P, ARG, PR, PRE> sf = _p->get(_i);
        sf.set_error(e);
        _p->set(_i, sf);
    }
    /** @brief Set the error. */
    Void set_error(Nat e) const {
        ScalarFunctionModel<P, ARG, PR, PRE> sf = _p->get(_i);
        sf.set_error(e);
        _p->set(_i, sf);
    }
    /** @brief Return the range enclosure. */
    const RangeType range() const {
        ScalarFunctionModel<P, ARG, PR, PRE> sf = _p->get(_i);
        return sf.range();
    }
    /** @brief Test whether the first enclosure is contained in the second. */
    friend Boolean refines(const ScalarFunctionModel<ValidatedTag, ARG, PR, PRE>& f1, const ScalarFunctionModel<ValidatedTag, ARG, PR, PRE>& f2);
    /** @brief Test whether the enclosures cannot represent a common value. */
    friend Boolean inconsistent(const ScalarFunctionModel<ValidatedTag, ARG, PR, PRE>& f1, const ScalarFunctionModel<ValidatedTag, ARG, PR, PRE>& f2);
    /** @brief Return an antiderivative with respect to the selected argument. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> antiderivative(const VectorFunctionModelElement<P, ARG, PR, PRE>& f, SizeType k) {
        return antiderivative(ScalarFunctionModel<P, ARG, PR, PRE>(f), k);
    }
    /** @brief Write the value to the output stream. */
    friend inline OutputStream& operator<<(OutputStream& os, const VectorFunctionModelElement<P, ARG, PR, PRE>& function) {
        return os << static_cast<const ScalarFunctionModel<P, ARG, PR, PRE>>(function);
    }
};

//! \ingroup FunctionModelSubModule
//! \brief Generic vector functions on bounded domains.
template<class P, class ARG, class PR, class PRE>
class FunctionModel<P, RealVector(ARG), PR, PRE>
    : public Handle<FunctionModelInterface<P, RealVector(ARG), PR, PRE>> {
    static_assert(Same<ARG, RealScalar> or Same<ARG, RealVector>, "");
    static_assert(Same<PRE, DoublePrecision> or Same<PRE, MultiplePrecision>, "");
    using RES = RealVector;
    using SIG = RES(ARG);
    using C = typename SignatureTraits<SIG>::BoundedCodomainType;
    using D = typename SignatureTraits<SIG>::BoundedDomainType;

  public:
    typedef FunctionModelInterface<P, SIG, PR, PRE> Interface;

  public:
    typedef VectorFunction<P, ARG> GenericType;
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef D DomainType;
    typedef C CodomainType;

    typedef typename Interface::RangeType RangeType;
    typedef typename Interface::NormType NormType;
    typedef typename Interface::ValueType ValueType;
    typedef typename Interface::ErrorType ErrorType;
    typedef typename Interface::NumericType NumericType;
    typedef typename Interface::GenericNumericType GenericNumericType;
    typedef typename Interface::ArgumentIndexType ArgumentIndexType;

    template<class Y>
    using Argument = typename SignatureTraits<SIG>::template Argument<Y>;
    template<class Y>
    using Result = typename SignatureTraits<SIG>::template Result<Y>;

  public:
    /** @brief Construct a \c FunctionModel from p. */
    inline explicit FunctionModel(Interface* p) : Handle<Interface>(p) {}
    /** @brief Construct a \c FunctionModel from p. */
    inline FunctionModel(UniquePointer<Interface> p) : Handle<Interface>(p->_concrete_clone()) {}

    /** @brief Construct a \c FunctionModel with the specified characteristics. */
    inline FunctionModel() : FunctionModel(nullptr) {}
    /** @brief Construct a \c FunctionModel from vfp. */
    inline FunctionModel(SharedPointer<const Interface> vfp)
        : FunctionModel(vfp->_concrete_clone()) {}
    /** @brief Construct a \c FunctionModel from n, sf. */
    inline FunctionModel(SizeType n, const ScalarFunctionModelInterface<P, ARG, PR, PRE>& sf) : FunctionModel() {
        FunctionModelFactory<P, PR, PRE> factory(sf._factory());
        *this = factory.create_zeros(n, sf.domain());
        for (SizeType i = 0; i != n; ++i) {
            (*this)[i] = sf;
        }
    }
    /** @brief Construct a \c FunctionModel from asf. */
    inline FunctionModel(const Array<ScalarFunctionModel<P, ARG, PR, PRE>>& asf)
        : FunctionModel(asf.size(), asf[0]) {
        for (SizeType i = 0; i != asf.size(); ++i) {
            (*this)[i] = asf[i];
        }
    }
    /** @brief Construct a \c FunctionModel from lsf. */
    inline FunctionModel(const List<ScalarFunctionModel<P, ARG, PR, PRE>>& lsf)
        : FunctionModel(lsf.size(), lsf[0]) {
        for (SizeType i = 0; i != lsf.size(); ++i) {
            (*this)[i] = lsf[i];
        }
    }
    /** @brief Construct a \c FunctionModel from f. */
    inline FunctionModel(const Interface& f) : FunctionModel(f._concrete_clone()) {}
    /** @brief Convert to function< p, sig >. */
    inline operator Function<P, SIG>() const { return Function<P, SIG>(*this->_ptr); }

    /** @brief Return the number of result components. */
    inline SizeType result_size() const { return this->_ptr->result_size(); }
    /** @brief Return the number of scalar arguments. */
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    /** @brief Return the number of elements. */
    inline SizeType size() const { return this->_ptr->result_size(); }
    /** @brief Apply the callable to the supplied arguments. */
    template<class XX>
    inline Vector<XX> operator()(const Vector<XX>& v) const { return this->_ptr->_call(v); }
    /** @brief Evaluate at the supplied argument. */
    template<class XX>
    inline Vector<XX> evaluate(const Vector<XX>& v) const { return this->_ptr->_call(v); }
    /** @brief Return the requested element or stored value. */
    inline const ScalarFunctionModel<P, ARG, PR, PRE> get(SizeType i) const { return ScalarFunctionModel<P, ARG, PR, PRE>(this->_ptr->_concrete_get(i)); }
    /** @brief Replace the selected value. */
    inline Void set(SizeType i, const ScalarFunctionModel<P, ARG, PR, PRE>& sf) { this->pointer()->_concrete_set(i, sf); }
    /** @brief Return the selected element or component. */
    inline const ScalarFunctionModel<P, ARG, PR, PRE> operator[](SizeType i) const { return this->get(i); }
    /** @brief Return the selected element or component. */
    inline VectorFunctionModelElement<P, ARG, PR, PRE> operator[](SizeType i) { return VectorFunctionModelElement<P, ARG, PR, PRE>(this, i); }
    /** @brief Return the selected element or component. */
    inline VectorFunctionModel<P, ARG, PR, PRE> operator[](Range rng) const {
        VectorFunctionModel<P, ARG, PR, PRE> r = factory(*this).create_zeros(rng.size());
        for (SizeType i = 0; i != rng.size(); ++i) {
            r[i] = this->operator[](rng[i]);
        }
        return r;
    }
    /** @brief Return the working precision. */
    inline const PrecisionType precision() const { return this->get(0).precision(); }
    /** @brief Return the domain. */
    inline const DomainType domain() const { return this->_ptr->domain(); }
    /** @brief Return the codomain. */
    inline const CodomainType codomain() const { return this->_ptr->codomain(); }
    /** @brief Return the range enclosure. */
    inline const RangeType range() const { return this->_ptr->range(); }
    /** @brief Return the stored values. */
    inline const Vector<ValueType> values() const { return this->_ptr->_concrete_values(); }
    /** @brief Return the component error bounds. */
    inline const Vector<ErrorType> errors() const { return this->_ptr->_concrete_errors(); }
    /** @brief Return the error bound. */
    inline const ErrorType error() const { return this->_ptr->_concrete_error(); }
    /** @brief Discard the error bound of the model. */
    inline Void clobber() { this->pointer()->clobber(); }
    /** @brief Return the matrix of first partial derivatives. */
    inline const Matrix<NumericType> jacobian(const Vector<NumericType>& x) const
        requires Same<ARG, RealVector>
    {
        Vector<Differential<NumericType>> dfx = this->_ptr->_call(Differential<NumericType>::variables(1u, x));
        return dfx.jacobian();
    }

    /** @brief Restrict to the supplied subset or domain. */
    inline Void restrict(const DomainType& d) { *this = restriction(*this, d); }

  public:
    /** @brief Return a factory for compatible objects. */
    friend FunctionModelCreator<FunctionModelFactory<P, PR, PRE>, ARG> factory(const VectorFunctionModel<P, ARG, PR, PRE>& f) {
        FunctionModelFactory<P, PR, PRE> factory(f._ptr->_factory());
        return FunctionModelCreator<FunctionModelFactory<P, PR, PRE>, ARG>(f.domain(), factory);
    }

  public:
    /** @brief Compose the outer object with the inner object. */
    friend inline ScalarFunctionModel<P, ARG, PR, PRE> compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionModel<P, ARG, PR, PRE>& g) {
        return ScalarFunctionModel<P, ARG, PR, PRE>(g._ptr->_compose(f));
    }
    /** @brief Compose the outer object with the inner object. */
    friend inline ScalarFunctionModel<P, ARG, PR, PRE> compose(const ScalarMultivariateFunctionPatch<P>& f, const VectorFunctionModel<P, ARG, PR, PRE>& g) {
        return ScalarFunctionModel<P, ARG, PR, PRE>(g._ptr->_compose(f));
    }
    /** @brief Compose the outer object with the inner object. */
    friend inline ScalarFunctionModel<P, ARG, PR, PRE> compose(const ScalarFunctionModel<P, ARG, PR, PRE>& f, const VectorFunctionModel<P, ARG, PR, PRE>& g)
        requires Same<ARG, RealVector>
    {
        return ScalarFunctionModel<P, ARG, PR, PRE>(g._ptr->_compose(f));
    }
    /** @brief Compose the outer object with the inner object. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> compose(const VectorMultivariateFunction<P>& f, const VectorFunctionModel<P, ARG, PR, PRE>& g) {
        return VectorFunctionModel<P, ARG, PR, PRE>(g._ptr->_compose(f));
    }
    /** @brief Compose the outer object with the inner object. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> compose(const VectorMultivariateFunctionPatch<P>& f, const VectorFunctionModel<P, ARG, PR, PRE>& g) {
        return VectorFunctionModel<P, ARG, PR, PRE>(g._ptr->_compose(f));
    }
    /** @brief Compose the outer object with the inner object. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> compose(const VectorFunctionModel<P, ARG, PR, PRE>& f, const VectorFunctionModel<P, ARG, PR, PRE>& g)
        requires Same<ARG, RealVector>
    {
        return VectorFunctionModel<P, ARG, PR, PRE>(g._ptr->_compose(f));
    }

    /** @brief Compose without checking domain compatibility. */
    friend inline ScalarFunctionModel<P, ARG, PR, PRE> unchecked_compose(const ScalarFunctionModel<P, RealVector, PR, PRE>& f, const VectorFunctionModel<P, ARG, PR, PRE>& g) {
        return ScalarFunctionModel<P, ARG, PR, PRE>(g._ptr->_unchecked_compose(static_cast<ScalarMultivariateFunction<P>>(f)));
    }
    /** @brief Compose without checking domain compatibility. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> unchecked_compose(const VectorFunctionModel<P, RealVector, PR, PRE>& f, const VectorFunctionModel<P, ARG, PR, PRE>& g) {
        return VectorFunctionModel<P, ARG, PR, PRE>(g._ptr->_unchecked_compose(static_cast<VectorMultivariateFunction<P>>(f)));
    }

    /** @brief Return the sum or unary positive value. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator+(const VectorFunctionModel<P, ARG, PR, PRE>& f) {
        return VectorFunctionModel<P, ARG, PR, PRE>(f._ptr->_concrete_clone());
    }
    /** @brief Return the difference or additive inverse. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator-(const VectorFunctionModel<P, ARG, PR, PRE>& f) {
        VectorFunctionModel<P, ARG, PR, PRE> r = f;
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = -f[i];
        }
        return r;
    }
    /** @brief Return the sum or unary positive value. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator+(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const VectorFunctionModel<P, ARG, PR, PRE>& f2) {
        VectorFunctionModel<P, ARG, PR, PRE> r = f1;
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = f1[i] + f2[i];
        }
        return r;
    }
    /** @brief Return the difference or additive inverse. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator-(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const VectorFunctionModel<P, ARG, PR, PRE>& f2) {
        VectorFunctionModel<P, ARG, PR, PRE> r = f1;
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = f1[i] - f2[i];
        }
        return r;
    }
    /** @brief Return the product. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator*(const ScalarFunctionModel<P, ARG, PR, PRE>& f1, const VectorFunctionModel<P, ARG, PR, PRE>& f2) {
        VectorFunctionModel<P, ARG, PR, PRE> r = f2;
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = f1 * f2[i];
        }
        return r;
    }
    /** @brief Return the product. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator*(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const ScalarFunctionModel<P, ARG, PR, PRE>& f2) {
        VectorFunctionModel<P, ARG, PR, PRE> r = f1;
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = f1[i] * f2;
        }
        return r;
    }
    /** @brief Return the quotient. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator/(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const ScalarFunctionModel<P, ARG, PR, PRE>& f2) {
        VectorFunctionModel<P, ARG, PR, PRE> r = f1;
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = f1[i] / f2;
        }
        return r;
    }
    /** @brief Return the sum or unary positive value. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator+(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const Vector<CanonicalNumericType<P, PR, PRE>>& c2) {
        VectorFunctionModel<P, ARG, PR, PRE> r = f1;
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = f1[i] + c2[i];
        }
        return r;
    }
    /** @brief Return the difference or additive inverse. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator-(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const Vector<CanonicalNumericType<P, PR, PRE>>& c2) {
        VectorFunctionModel<P, ARG, PR, PRE> r = f1;
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = f1[i] - c2[i];
        }
        return r;
    }
    /** @brief Return the sum or unary positive value. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator+(const Vector<CanonicalNumericType<P, PR, PRE>>& c1, const VectorFunctionModel<P, ARG, PR, PRE>& f2);
    /** @brief Return the difference or additive inverse. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator-(const Vector<CanonicalNumericType<P, PR, PRE>>& c1, const VectorFunctionModel<P, ARG, PR, PRE>& f2);
    /** @brief Return the product. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator*(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const CanonicalNumericType<P, PR, PRE>& c2) {
        VectorFunctionModel<P, ARG, PR, PRE> r = f1;
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = f1[i] * c2;
        }
        return r;
    }
    /** @brief Return the product. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator*(const CanonicalNumericType<P, PR, PRE>& c1, const VectorFunctionModel<P, ARG, PR, PRE>& f2) {
        VectorFunctionModel<P, ARG, PR, PRE> r = f2;
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = c1 * f2[i];
        }
        return r;
    }

    /** @brief Return the sum or unary positive value. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator+(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const VectorMultivariateFunction<P>& f2)
        requires Same<ARG, RealVector>
    {
        return f1 + factory(f1).create(f2);
    }
    /** @brief Return the difference or additive inverse. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator-(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const VectorMultivariateFunction<P>& f2)
        requires Same<ARG, RealVector>
    {
        return f1 - factory(f1).create(f2);
    }
    /** @brief Return the product. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator*(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const ScalarMultivariateFunction<P>& f2)
        requires Same<ARG, RealVector>
    {
        return f1 * factory(f1).create(f2);
    }
    /** @brief Return the quotient. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator/(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const ScalarMultivariateFunction<P>& f2)
        requires Same<ARG, RealVector>
    {
        return f1 / factory(f1).create(f2);
    }
    /** @brief Return the sum or unary positive value. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator+(const VectorMultivariateFunction<P>& f1, const VectorFunctionModel<P, ARG, PR, PRE>& f2)
        requires Same<ARG, RealVector>
    {
        return factory(f2).create(f1) + f2;
    }
    /** @brief Return the difference or additive inverse. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator-(const VectorMultivariateFunction<P>& f1, const VectorFunctionModel<P, ARG, PR, PRE>& f2)
        requires Same<ARG, RealVector>
    {
        return factory(f2).create(f1) - f2;
    }
    /** @brief Return the product. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator*(const ScalarMultivariateFunction<P>& f1, const VectorFunctionModel<P, ARG, PR, PRE>& f2)
        requires Same<ARG, RealVector>
    {
        return factory(f2).create(f1) * f2;
    }
    /** @brief Return the quotient. */
    friend inline VectorFunctionModel<P, ARG, PR, PRE> operator/(const ScalarMultivariateFunction<P>& f1, const VectorFunctionModel<P, ARG, PR, PRE>& f2)
        requires Same<ARG, RealVector>
    {
        return factory(f2).create(f1) / f2;
    }

  public:
    /** @brief Return the norm. */
    friend NormType norm(const VectorFunctionModel<P, ARG, PR, PRE>& f) {
        return f._ptr->_concrete_norm();
    }
    /** @brief Embed in the larger argument space. */
    friend VectorFunctionModel<P, ARG, PR, PRE> embed(const DomainType& d1, const VectorFunctionModel<P, ARG, PR, PRE>& f, const DomainType& d2)
        requires Same<ARG, RealVector>
    {
        return VectorFunctionModel<P, ARG, PR, PRE>(f._ptr->_embed(d1, d2));
    }
    /** @brief Embed in the larger argument space. */
    friend VectorFunctionModel<P, ARG, PR, PRE> embed(const DomainType& d, const VectorFunctionModel<P, ARG, PR, PRE>& f)
        requires Same<ARG, RealVector>
    {
        return embed(d, f, DomainType());
    }
    /** @brief Embed in the larger argument space. */
    friend VectorFunctionModel<P, ARG, PR, PRE> embed(const VectorFunctionModel<P, ARG, PR, PRE>& f, const BoxDomainType& d)
        requires Same<ARG, RealVector>
    {
        return embed(DomainType(), f, d);
    }
    /** @brief Embed in the larger argument space. */
    friend VectorFunctionModel<P, ARG, PR, PRE> embed(const VectorFunctionModel<P, ARG, PR, PRE>& f, const IntervalDomainType& d)
        requires Same<ARG, RealVector>
    {
        return embed(f, DomainType(1, d));
    }
    /** @brief Return the object restricted to the supplied domain. */
    friend VectorFunctionModel<P, ARG, PR, PRE> restriction(const VectorFunctionModel<P, ARG, PR, PRE>& f, const DomainType& d) {
        return VectorFunctionModel<P, ARG, PR, PRE>(f._ptr->_restriction(d));
    }

    /** @brief Evaluate at the supplied argument. */
    friend Vector<CanonicalNumericType<P, PR, PRE>> evaluate(const VectorFunctionModel<P, ARG, PR, PRE>& f, const Argument<CanonicalNumericType<P, PR, PRE>>& x) {
        return f._ptr->_call(x);
    }
    /** @brief Evaluate at the supplied argument. */
    friend Vector<Number<P>> evaluate(const VectorFunctionModel<P, ARG, PR, PRE>& f, const Argument<Number<P>>& x) {
        return f._ptr->_call(Argument<CanonicalNumericType<P, PR, PRE>>(x, f.precision()));
    }

    /** @brief Evaluate without checking that the argument belongs to the domain. */
    friend Vector<CanonicalNumericType<P, PR, PRE>> unchecked_evaluate(const VectorFunctionModel<P, ARG, PR, PRE>& f, const Argument<CanonicalNumericType<P, PR, PRE>>& x) {
        return f._ptr->_unchecked_evaluate(x);
    }
    /** @brief Evaluate without checking that the argument belongs to the domain. */
    friend Vector<Number<P>> unchecked_evaluate(const VectorFunctionModel<P, ARG, PR, PRE>& f, const Argument<Number<P>>& x) {
        return f._ptr->_unchecked_evaluate(Argument<CanonicalNumericType<P, PR, PRE>>(x, f.precision()));
    }

    /** @brief Fix the selected argument to the supplied value. */
    friend VectorFunctionModel<P, ARG, PR, PRE> partial_evaluate(const VectorFunctionModel<P, ARG, PR, PRE>& f, SizeType j, const CanonicalNumericType<P, PR, PRE>& c) {
        return VectorFunctionModel<P, ARG, PR, PRE>(f._ptr->_concrete_partial_evaluate(j, c));
    }
    /** @brief Fix the selected argument to the supplied value. */
    friend VectorFunctionModel<P, ARG, PR, PRE> partial_evaluate(const VectorFunctionModel<P, ARG, PR, PRE>& f, SizeType j, const Number<P>& c) {
        return partial_evaluate(f, j, CanonicalNumericType<P, PR, PRE>(c, f.precision()));
    }

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const VectorFunctionModel<P, ARG, PR, PRE>& f) {
        return os << static_cast<GenericType>(f);
    }

    /** @brief Compose without checking domain compatibility. */
    friend ScalarFunctionModel<P, ARG, PR, PRE> unchecked_compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionModel<P, ARG, PR, PRE>& g) {
        const ScalarFunctionModelInterface<P, RealVector, PR, PRE>* fptr = dynamic_cast<const ScalarFunctionModelInterface<P, RealVector, PR, PRE>*>(f.raw_pointer());
        if (fptr) {
            return unchecked_compose(ScalarFunctionModel<P, RealVector, PR, PRE>(*fptr), g);
        } else {
            return compose(f, g);
        }
    }
    /** @brief Compose without checking domain compatibility. */
    friend VectorFunctionModel<P, ARG, PR, PRE> unchecked_compose(const VectorMultivariateFunction<P>& f, const VectorFunctionModel<P, ARG, PR, PRE>& g) {
        const VectorFunctionModelInterface<P, RealVector, PR, PRE>* fptr = dynamic_cast<const VectorFunctionModelInterface<P, RealVector, PR, PRE>*>(f.raw_pointer());
        if (fptr) {
            return unchecked_compose(VectorFunctionModel<P, RealVector, PR, PRE>(*fptr), g);
        } else {
            return compose(f, g);
        }
    }

    /** @brief Join the operands in their natural order. */
    friend VectorFunctionModel<P, ARG, PR, PRE> join(const ScalarFunctionModel<P, ARG, PR, PRE>& f1, const VectorFunctionModel<P, ARG, PR, PRE>& f2) {
        return join(VectorFunctionModel<P, ARG, PR, PRE>(1u, f1), f2);
    }
    /** @brief Join the operands in their natural order. */
    friend VectorFunctionModel<P, ARG, PR, PRE> join(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const ScalarFunctionModel<P, ARG, PR, PRE>& f2) {
        VectorFunctionModel<P, ARG, PR, PRE> r = VectorFunctionModel<P, ARG, PR, PRE>(f1._ptr->_concrete_clone());
        r._ptr->_concrete_adjoin(f2);
        return r;
    }
    /** @brief Join the operands in their natural order. */
    friend VectorFunctionModel<P, ARG, PR, PRE> join(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const VectorFunctionModel<P, ARG, PR, PRE>& f2) {
        return VectorFunctionModel<P, ARG, PR, PRE>(f1._ptr->_concrete_join(*f2._ptr));
    }

    /** @brief Combine the operands over their respective arguments. */
    friend VectorFunctionModel<P, ARG, PR, PRE> combine(const ScalarFunctionModel<P, ARG, PR, PRE>& f1, const ScalarFunctionModel<P, ARG, PR, PRE>& f2) {
        return VectorFunctionModel<P, ARG, PR, PRE>(VectorFunctionModel<P, ARG, PR, PRE>(1, f1)._ptr->_concrete_combine(VectorFunctionModel<P, ARG, PR, PRE>(1, f2)));
    };
    /** @brief Combine the operands over their respective arguments. */
    friend VectorFunctionModel<P, ARG, PR, PRE> combine(const ScalarFunctionModel<P, ARG, PR, PRE>& f1, const VectorFunctionModel<P, ARG, PR, PRE>& f2) {
        return VectorFunctionModel<P, ARG, PR, PRE>(VectorFunctionModel<P, ARG, PR, PRE>(1, f1)._ptr->_concrete_combine(f2));
    };
    /** @brief Combine the operands over their respective arguments. */
    friend VectorFunctionModel<P, ARG, PR, PRE> combine(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const ScalarFunctionModel<P, ARG, PR, PRE>& f2) {
        return VectorFunctionModel<P, ARG, PR, PRE>(f1._ptr->_concrete_combine(VectorFunctionModel<P, ARG, PR, PRE>(1, f2)));
    };
    /** @brief Combine the operands over their respective arguments. */
    friend VectorFunctionModel<P, ARG, PR, PRE> combine(const VectorFunctionModel<P, ARG, PR, PRE>& f1, const VectorFunctionModel<P, ARG, PR, PRE>& f2) {
        return VectorFunctionModel<P, ARG, PR, PRE>(f1._ptr->_concrete_combine(f2));
    }

    /** @brief Return the common refinement of the enclosures. */
    friend inline VectorFunctionModel<ValidatedTag, ARG, PR, PRE> refinement(const VectorFunctionModel<ValidatedTag, ARG, PR, PRE>& f1, const VectorFunctionModel<ValidatedTag, ARG, PR, PRE>& f2) {
        ARIADNE_ASSERT_MSG(f1.size() == f2.size(), "refinement(f1,f2): f1=" << f1 << ", f2=" << f2 << ")");
        VectorFunctionModel<ValidatedTag, ARG, PR, PRE> r = +f1;
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = refinement(f1[i], f2[i]);
        }
        return r;
    }

    /** @brief Return an antiderivative with respect to the selected argument. */
    friend VectorFunctionModel<P, ARG, PR, PRE> antiderivative(const VectorFunctionModel<P, ARG, PR, PRE>& f, ArgumentIndexType j) {
        VectorFunctionModel<P, ARG, PR, PRE> r(f);
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = antiderivative(f[i], j);
        }
        return r;
    }

    /** @brief Return an antiderivative with respect to the selected argument. */
    friend VectorFunctionModel<P, ARG, PR, PRE> antiderivative(const VectorFunctionModel<P, ARG, PR, PRE>& f, ArgumentIndexType j, const CanonicalNumericType<P, PR, PRE>& c) {
        VectorFunctionModel<P, ARG, PR, PRE> r(f);
        for (SizeType i = 0; i != r.size(); ++i) {
            r[i] = antiderivative(f[i], j, c);
        }
        return r;
    }

    /** @brief Return an antiderivative with respect to the selected argument. */
    friend VectorFunctionModel<P, ARG, PR, PRE> antiderivative(const VectorFunctionModel<P, ARG, PR, PRE>& f, ArgumentIndexType j, const Number<P>& c) {
        return antiderivative(f, j, CanonicalNumericType<P, PR, PRE>(c, f[0].value().precision()));
    }
};

// FIXME: Implement for Multiple-Precision versions
/** @brief Evaluate without checking that the argument belongs to the domain. */
template<class P>
inline CanonicalNumeric64Type<P> unchecked_evaluate(const ScalarMultivariateFunction<P>& f, const Vector<CanonicalNumeric64Type<P>>& x) {
    const auto* fmptr = dynamic_cast<const typename ScalarMultivariateFunctionModelDP<P>::Interface*>(f.raw_pointer());
    if (fmptr) {
        return unchecked_evaluate(ScalarMultivariateFunctionModelDP<P>(*fmptr), x);
    }
    const auto* fptr = dynamic_cast<const typename ScalarMultivariateFunctionPatch<P>::Interface*>(f.raw_pointer());
    if (fptr) {
        return unchecked_evaluate(ScalarMultivariateFunctionPatch<P>(*fptr), x);
    }
    return evaluate(f, x);
}

/** @brief Evaluate without checking that the argument belongs to the domain. */
template<class P>
inline Vector<CanonicalNumeric64Type<P>> unchecked_evaluate(const VectorMultivariateFunction<P>& f, const Vector<CanonicalNumeric64Type<P>>& x) {
    const auto* fmptr = dynamic_cast<const typename VectorMultivariateFunctionModelDP<P>::Interface*>(f.raw_pointer());
    if (fmptr) {
        return unchecked_evaluate(VectorMultivariateFunctionModelDP<P>(*fmptr), x);
    }
    const auto* fptr = dynamic_cast<const typename VectorMultivariateFunctionPatch<P>::Interface*>(f.raw_pointer());
    if (fptr) {
        return unchecked_evaluate(VectorMultivariateFunctionPatch<P>(*fptr), x);
    }
    return evaluate(f, x);
}

/** @brief Compose without checking domain compatibility. */
template<class P, class ARG, class PR, class PRE>
inline ScalarFunctionModel<P, ARG, PR, PRE> unchecked_compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionModel<P, ARG, PR, PRE>& g) {
    const auto* fptr = dynamic_cast<const typename ScalarFunctionModel<P, ARG, PR, PRE>::Interface*>(f.raw_pointer());
    if (fptr) {
        return unchecked_compose(ScalarFunctionModel<P, ARG, PR, PRE>(*fptr), g);
    } else {
        return compose(f, g);
    }
}
/** @brief Compose without checking domain compatibility. */
template<class P, class ARG, class PR, class PRE>
inline VectorFunctionModel<P, ARG, PR, PRE> unchecked_compose(const VectorMultivariateFunction<P>& f, const VectorFunctionModel<P, ARG, PR, PRE>& g) {
    const auto* fptr = dynamic_cast<const typename VectorFunctionModel<P, ARG, PR, PRE>::Interface*>(f.raw_pointer());
    if (fptr) {
        return unchecked_compose(VectorFunctionModel<P, ARG, PR, PRE>(*fptr), g);
    } else {
        return compose(f, g);
    }
}

// FIXME: Should be unneeded
/** @brief Compose without checking domain compatibility. */
template<class ARG, class PR, class PRE>
ScalarFunctionModel<ValidatedTag, ARG, PR, PRE> unchecked_compose(const ScalarMultivariateFunctionModel<ValidatedTag, PR, PRE>& f, const VectorFunctionModel<ValidatedTag, ARG, PR, PRE>& g) {
    return VectorFunctionModel<ValidatedTag, ARG, PR, PRE>(g._ptr->_unchecked_compose(f));
}
/** @brief Compose without checking domain compatibility. */
template<class ARG, class PR, class PRE>
VectorFunctionModel<ValidatedTag, ARG, PR, PRE> unchecked_compose(const VectorMultivariateFunctionModel<ValidatedTag, PR, PRE>& f, const VectorFunctionModel<ValidatedTag, ARG, PR, PRE>& g) {
    return VectorFunctionModel<ValidatedTag, ARG, PR, PRE>(g._ptr->_unchecked_compose(f));
}

template<class P, class ARG, class PR, class PRE>
VectorFunctionModelElement<P, ARG, PR, PRE>::operator const ScalarFunctionModel<P, ARG, PR, PRE>() const {
    return ScalarFunctionModel<P, ARG, PR, PRE>(_p->get(_i));
}

template<class P, class... ARGS>
/** @brief Construct a \c FunctionPatch from fm. */
template<class PR, class PRE>
FunctionPatch<P, RealScalar(ARGS...)>::FunctionPatch(FunctionModel<P, RealScalar(ARGS...), PR, PRE> fm)
    : _ptr(fm.pointer()) {}

template<class P, class... ARGS>
/** @brief Construct a \c FunctionPatch from fm. */
template<class PR, class PRE>
FunctionPatch<P, RealVector(ARGS...)>::FunctionPatch(FunctionModel<P, RealVector(ARGS...), PR, PRE> fm)
    : _ptr(fm.pointer()) {}

} // namespace Ariadne

#endif
