/***************************************************************************
 *            function/function_archetypes.hpp
 *
 *  Copyright  2022  Pieter Collins
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

/** @file function/function_archetypes.hpp
 *  @brief Archetypes implementing abstract function concepts.
 */

#ifndef ARIADNE_FUNCTION_ARCHETYPES_HPP
#define ARIADNE_FUNCTION_ARCHETYPES_HPP

#include "function/formula.decl.hpp"

namespace Ariadne {

/** @brief Reference type specifying the requirements of function. */
template<class P, class SIG>
class FunctionArchetype;

/** @brief Reference type specifying the requirements of function. */
template<class SIG>
class FunctionArchetype<ApproximateTag, SIG> {
    using P = ApproximateTag;

  public:
    using DomainType = typename SignatureTraits<SIG>::DomainType;
    using CodomainType = typename SignatureTraits<SIG>::CodomainType;
    using ArgumentSizeType = typename SignatureTraits<SIG>::ArgumentSizeType;
    using ResultSizeType = typename SignatureTraits<SIG>::ResultSizeType;
    using ArgumentIndexType = typename SignatureTraits<SIG>::ArgumentIndexType;
    template<class X>
    using Argument = typename SignatureTraits<SIG>::template Argument<X>;
    template<class X>
    using Result = typename SignatureTraits<SIG>::template Result<X>;

    /** @brief Return the number of scalar arguments. */
    ArgumentSizeType argument_size() const { assert(false); }
    /** @brief Return the number of result components. */
    ResultSizeType result_size() const { assert(false); }

    /** @brief Apply the callable to the supplied arguments. */
    Result<FloatDPApproximation> operator()(const Argument<FloatDPApproximation>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Result<FloatMPApproximation> operator()(const Argument<FloatMPApproximation>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Result<Differential<FloatDPApproximation>> operator()(const Argument<Differential<FloatDPApproximation>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Result<Differential<FloatMPApproximation>> operator()(const Argument<Differential<FloatMPApproximation>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Result<TaylorModel<ApproximateTag, FloatDP>> operator()(const Argument<TaylorModel<ApproximateTag, FloatDP>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Result<TaylorModel<ApproximateTag, FloatMP>> operator()(const Argument<TaylorModel<ApproximateTag, FloatMP>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Result<Formula<ApproximateNumber>> operator()(const Argument<Formula<ApproximateNumber>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Result<ElementaryAlgebra<ApproximateNumber>> operator()(const Argument<ElementaryAlgebra<ApproximateNumber>>&) const { assert(false); }

    /** @brief Return the derivative with respect to the selected argument. */
    friend FunctionArchetype<P, SIG> derivative(FunctionArchetype<P, SIG>, SizeType) { assert(false); }

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream&, const FunctionArchetype<P, SIG>&) { assert(false); }

    /** @brief Evaluate with the supplied arguments. */
    template<class X>
    Result<X> _call(const Argument<X>& x) const { return (*this)(x); }
};

/** @brief Reference type specifying the requirements of function. */
template<class SIG>
class FunctionArchetype<ValidatedTag, SIG>
    : public FunctionArchetype<ApproximateTag, SIG> {
    using P = ValidatedTag;

  public:
    using FunctionArchetype<ApproximateTag, SIG>::operator();
    using FunctionArchetype<ApproximateTag, SIG>::_call;

    /** @brief Apply the callable to the supplied arguments. */
    Scalar<FloatDPBounds> operator()(const Vector<FloatDPBounds>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Scalar<FloatMPBounds> operator()(const Vector<FloatMPBounds>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Scalar<Differential<FloatDPBounds>> operator()(const Vector<Differential<FloatDPBounds>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Scalar<Differential<FloatMPBounds>> operator()(const Vector<Differential<FloatMPBounds>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Scalar<TaylorModel<ValidatedTag, FloatDP>> operator()(const Vector<TaylorModel<ValidatedTag, FloatDP>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Scalar<TaylorModel<ValidatedTag, FloatMP>> operator()(const Vector<TaylorModel<ValidatedTag, FloatMP>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Scalar<TaylorModel<ValidatedTag, FloatDPUpperInterval>> operator()(const Vector<TaylorModel<ValidatedTag, FloatDPUpperInterval>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Scalar<TaylorModel<ValidatedTag, FloatMPUpperInterval>> operator()(const Vector<TaylorModel<ValidatedTag, FloatMPUpperInterval>>&) const { assert(false); }

    /** @brief Apply the callable to the supplied arguments. */
    Scalar<Formula<ValidatedNumber>> operator()(const Vector<Formula<ValidatedNumber>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Scalar<ElementaryAlgebra<ValidatedNumber>> operator()(const Vector<ElementaryAlgebra<ValidatedNumber>>&) const { assert(false); }

    /** @brief Apply the callable to the supplied arguments. */
    Scalar<Function<ValidatedTag, SIG>> operator()(const Vector<Function<ValidatedTag, SIG>>&) const { assert(false); }

    /** @brief Apply the callable to the supplied arguments. */
    inline Scalar<FloatDPBounds> operator()(const Vector<FloatDP>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    inline Scalar<FloatMPBounds> operator()(const Vector<FloatMP>&) const { assert(false); }

    /** @brief Return the derivative with respect to the selected argument. */
    friend FunctionArchetype<P, SIG> derivative(FunctionArchetype<P, SIG>, SizeType) { assert(false); }
};

/** @brief Reference type specifying the requirements of function. */
template<class SIG>
class FunctionArchetype<EffectiveTag, SIG>
    : public FunctionArchetype<ValidatedTag, SIG> {
    using P = EffectiveTag;

  public:
    using FunctionArchetype<ValidatedTag, SIG>::operator();
    using FunctionArchetype<ValidatedTag, SIG>::_call;

    /** @brief Apply the callable to the supplied arguments. */
    Scalar<Real> operator()(const Vector<Real>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Scalar<ElementaryAlgebra<Real>> operator()(const Vector<ElementaryAlgebra<Real>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Scalar<Formula<Real>> operator()(const Vector<Formula<Real>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Scalar<ElementaryAlgebra<EffectiveNumber>> operator()(const Vector<ElementaryAlgebra<EffectiveNumber>>&) const { assert(false); }
    /** @brief Apply the callable to the supplied arguments. */
    Scalar<Formula<EffectiveNumber>> operator()(const Vector<Formula<EffectiveNumber>>&) const { assert(false); }

    /** @brief Return the derivative with respect to the selected argument. */
    friend FunctionArchetype<P, SIG> derivative(FunctionArchetype<P, SIG>, SizeType) { assert(false); }
};

} // namespace Ariadne

#include "function_concepts.hpp"

namespace Ariadne {

using ScalarMultivariate = RealScalar(RealVector);

static_assert(IsApproximateFunction<FunctionArchetype<ApproximateTag, ScalarMultivariate>, ScalarMultivariate>);
static_assert(IsValidatedFunction<FunctionArchetype<ValidatedTag, ScalarMultivariate>, ScalarMultivariate>);
static_assert(IsEffectiveFunction<FunctionArchetype<EffectiveTag, ScalarMultivariate>, ScalarMultivariate>);

static_assert(AFunction<FunctionArchetype<ApproximateTag, ScalarMultivariate>, ApproximateTag, ScalarMultivariate>);
static_assert(AFunction<FunctionArchetype<ValidatedTag, ScalarMultivariate>, ValidatedTag, ScalarMultivariate>);
static_assert(AFunction<FunctionArchetype<EffectiveTag, ScalarMultivariate>, EffectiveTag, ScalarMultivariate>);

} // namespace Ariadne

#endif // ARIADNE_FUNCTION_ARCHETYPES_HPP
