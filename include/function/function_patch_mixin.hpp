/***************************************************************************
 *            function/function_patch_mixin.hpp
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

/*! @file function/function_patch_mixin.hpp
 *  \brief Mixin for functions on bounded sets.
 */

#ifndef ARIADNE_FUNCTION_PATCH_MIXIN_HPP
#define ARIADNE_FUNCTION_PATCH_MIXIN_HPP

#include "../function/domain.hpp"
#include "../function/function.decl.hpp"
#include "../function/function_interface.hpp"
#include "../function/function_mixin.hpp"
#include "../function/function_patch_interface.hpp"
#include "algebra/algebra_interface.hpp"
#include "algebra/algebra_mixin.hpp"
#include "algebra/matrix.hpp"
#include "algebra/operations.hpp"
#include "algebra/vector.hpp"
#include "numeric/number.decl.hpp"
#include "numeric/numeric.hpp"
#include "numeric/operators.hpp"
#include "utility/range.hpp"

namespace Ariadne {

/** @brief Reusable implementation of function patch. */
template<class FP, class P, class SIG>
class FunctionPatchMixin;
template<class FP, class P, class... ARGS>
using ScalarFunctionPatchMixin = FunctionPatchMixin<FP, P, RealScalar(ARGS...)>;
template<class FP, class P, class... ARGS>
using VectorFunctionPatchMixin = FunctionPatchMixin<FP, P, RealVector(ARGS...)>;
template<class FP, class P>
using ScalarUnivariateFunctionPatchMixin = FunctionPatchMixin<FP, P, RealScalar(RealScalar)>;
template<class FP, class P>
using VectorUnivariateFunctionPatchMixin = FunctionPatchMixin<FP, P, RealVector(RealScalar)>;
template<class FP, class P>
using ScalarMultivariateFunctionPatchMixin = FunctionPatchMixin<FP, P, RealScalar(RealVector)>;
template<class FP, class P>
using VectorMultivariateFunctionPatchMixin = FunctionPatchMixin<FP, P, RealVector(RealVector)>;

/** @brief Reusable implementation of function patch. */
template<class FP, class P, class... ARGS>
class FunctionPatchMixin<FP, P, RealScalar(ARGS...)>
    : public FunctionMixin<FP, P, RealScalar(ARGS...), ScalarFunctionPatchInterface<P, ARGS...>> {
    using RES = RealScalar;
    using SIG = RES(ARGS...);
    using D = typename ScalarFunctionPatchInterface<P, ARGS...>::DomainType;
    using C = typename ScalarFunctionPatchInterface<P, ARGS...>::CodomainType;

  public:
    typedef D DomainType;
    typedef C CodomainType;
    using typename FunctionPatchInterface<P, SIG>::ErrorType;
    using typename FunctionPatchInterface<P, SIG>::NormType;
    typedef typename SignatureTraits<SIG>::BoundedRangeType RangeType;
    typedef Number<P> X;

  public:
    /** @brief Return the domain. */
    const DomainType domain() const override { return static_cast<FP const&>(*this).domain(); }
    /** @brief Return the codomain. */
    const CodomainType codomain() const override { return static_cast<FP const&>(*this).codomain(); }

    /** @brief Return a dynamically allocated copy preserving the concrete type using the concrete implementation. */
    ScalarFunctionPatchInterface<P, ARGS...>* _concrete_clone() const override {
        return new FP(static_cast<const FP&>(*this));
    }

    /** @brief Construct a compatible copy of the supplied object. */
    ScalarFunctionPatchInterface<P, ARGS...>* _create_copy() const override {
        return new FP(static_cast<const FP&>(*this));
    }
    /** @brief Construct a zero object with compatible characteristics. */
    ScalarFunctionPatchInterface<P, ARGS...>* _create_zero() const override {
        return new FP(factory(static_cast<FP const&>(*this)).create_zero());
    }
    /** @brief Construct a constant object with the supplied value. */
    ScalarFunctionPatchInterface<P, ARGS...>* _create_constant(const Number<P>& c) const override {
        return new FP(factory(static_cast<FP const&>(*this)).create_constant(c));
    }

    /** @brief Apply the supplied operation to the operands. */
    ScalarFunctionPatchInterface<P, ARGS...>* _apply(UnaryElementaryOperator op) const override {
        return new FP(op(static_cast<FP const&>(*this)));
    }
    /** @brief Apply the supplied operation to the operands. */
    ScalarFunctionPatchInterface<P, ARGS...>* _apply(GradedElementaryOperator op, Int n) const override {
        return new FP(op(static_cast<FP const&>(*this), n));
    }
    /** @brief Apply the supplied operation to the operands. */
    ScalarFunctionPatchInterface<P, ARGS...>* _apply(
        BinaryElementaryOperator op,
        const ElementaryAlgebraInterface<Number<P>>& other) const override {
        const auto* ptr = dynamic_cast<FP const*>(&other);
        ARIADNE_ASSERT(ptr);
        return new FP(op(static_cast<FP const&>(*this), *ptr));
    }
    /** @brief Apply the supplied operation to the operands. */
    ScalarFunctionPatchInterface<P, ARGS...>* _apply(
        BinaryElementaryOperator op, const Number<P>& value) const override {
        return new FP(op(static_cast<FP const&>(*this), value));
    }
    /** @brief Apply the operation with reversed operand order. */
    ScalarFunctionPatchInterface<P, ARGS...>* _rapply(
        BinaryElementaryOperator op, const Number<P>& value) const override {
        return new FP(op(value, static_cast<FP const&>(*this)));
    }

    /** @brief Return the range enclosure. */
    const RangeType _range() const override {
        return RangeType(static_cast<const FP&>(*this).range());
    }
    /** @brief Return the error bound. */
    const ErrorType _error() const override {
        return static_cast<ErrorType>(static_cast<const FP&>(*this).error());
    }
    /** @brief Return the component error bounds. */
    const Scalar<ErrorType> _errors() const override {
        return static_cast<ErrorType>(static_cast<const FP&>(*this).error());
    }

    /** @brief Discard the error bound of the model. */
    inline virtual Void _clobber() override {
        static_cast<FP&>(*this).clobber();
    }
    /** @brief Test whether the first enclosure is contained in the second. */
    inline virtual Bool _refines(const FunctionPatchInterface<P, SIG>& fp) const override {
        ARIADNE_ASSERT(dynamic_cast<const FP*>(&fp));
        return refines(static_cast<const FP&>(*this), dynamic_cast<const FP&>(fp));
    }
    /** @brief Compose the outer object with the inner object. */
    inline virtual ScalarFunctionPatchInterface<P, ARGS...>* _compose(const ScalarFunction<P, RES>& f) const override {
        return heap_move(compose(f, static_cast<const FP&>(*this)));
    }
    /** @brief Compose the outer object with the inner object. */
    inline virtual VectorFunctionPatchInterface<P, ARGS...>* _compose(const VectorFunction<P, RES>& f) const override {
        return heap_move(compose(f, static_cast<const FP&>(*this)));
    }
    /** @brief Compose without checking domain compatibility. */
    inline virtual ScalarFunctionPatchInterface<P, ARGS...>* _unchecked_compose(const ScalarFunction<P, RES>& f) const override {
        return this->_compose(f);
    }
    /** @brief Compose without checking domain compatibility. */
    inline virtual VectorFunctionPatchInterface<P, ARGS...>* _unchecked_compose(const VectorFunction<P, RES>& f) const override {
        return this->_compose(f);
    }

    /** @brief Return the norm. */
    const NormType _norm() const override {
        return static_cast<NormType>(norm(static_cast<const FP&>(*this)));
    }
    /** @brief Return the derivative with respect to the selected argument using the concrete implementation. */
    ScalarFunctionPatchInterface<P, ARGS...>* _concrete_derivative(SizeType j) const override {
        return new FP(derivative(static_cast<const FP&>(*this), j));
    }
    /** @brief Return an antiderivative with respect to the selected argument. */
    ScalarFunctionPatchInterface<P, ARGS...>* _antiderivative(SizeType j) const override {
        return new FP(antiderivative(static_cast<const FP&>(*this), j));
    }
    /** @brief Return an antiderivative with respect to the selected argument. */
    ScalarFunctionPatchInterface<P, ARGS...>* _antiderivative(SizeType j, Number<P> c) const override {
        return new FP(antiderivative(static_cast<const FP&>(*this), j, c));
    }
    /** @brief Return the object restricted to the supplied domain. */
    ScalarFunctionPatchInterface<P, ARGS...>* _restriction(const BoxDomainType& d) const override {
        return new FP(restriction(static_cast<const FP&>(*this), d));
    }
    /** @brief Evaluate without checking that the argument belongs to the domain. */
    Number<P> _unchecked_evaluate(const Vector<Number<P>>& x) const override {
        return unchecked_evaluate(static_cast<const FP&>(*this), x);
    }
    /** @brief Evaluate without checking that the argument belongs to the domain. */
    CanonicalNumericType<P, DP> _unchecked_evaluate(const Vector<CanonicalNumericType<P, DP>>& x) const override {
        auto r = unchecked_evaluate(static_cast<const FP&>(*this), x);
        if constexpr (Convertible<decltype(r), CanonicalNumericType<P, DP>>) {
            return r;
        } else {
            ARIADNE_FAIL_MSG("");
        }
    }
    /** @brief Evaluate without checking that the argument belongs to the domain. */
    CanonicalNumericType<P, MultiplePrecision> _unchecked_evaluate(const Vector<CanonicalNumericType<P, MultiplePrecision>>& x) const override {
        auto r = unchecked_evaluate(static_cast<const FP&>(*this), x);
        if constexpr (Convertible<decltype(r), CanonicalNumericType<P, MultiplePrecision>>) {
            return r;
        } else {
            ARIADNE_FAIL_MSG("");
        }
    }
    /** @brief Fix the selected argument to the supplied value. */
    ScalarFunctionPatchInterface<P, ARGS...>* _partial_evaluate(SizeType j, const Number<P>& c) const override {
        return heap_copy(partial_evaluate(static_cast<const FP&>(*this), j, c));
    }
    /** @brief Embed in the larger argument space. */
    ScalarFunctionPatchInterface<P, ARGS...>* _embed(const BoxDomainType& d1, const BoxDomainType& d2) const override {
        return new FP(embed(d1, static_cast<const FP&>(*this), d2));
    }

    /** @brief Write the value to the output stream. */
    OutputStream& _write(OutputStream& os) const override {
        return os << static_cast<FP const&>(*this);
    }
};

/** @brief Reusable implementation of function patch. */
template<class FP, class P, class... ARGS>
class FunctionPatchMixin<FP, P, RealVector(ARGS...)>
    : public FunctionMixin<FP, P, RealVector(ARGS...), VectorFunctionPatchInterface<P, ARGS...>>,
      public FunctionGetterMixin<FP, P, RealVector(ARGS...)> {
    using RES = RealVector;
    using SIG = RES(ARGS...);
    using C = BoxDomainType;
    using D = typename ScalarFunctionPatchInterface<P, ARGS...>::DomainType;

  public:
    typedef D DomainType;
    typedef C CodomainType;
    using typename FunctionPatchInterface<P, SIG>::ErrorType;
    using typename FunctionPatchInterface<P, SIG>::NormType;
    typedef typename VectorFunctionPatchInterface<P, ARGS...>::RangeType RangeType;
    typedef typename Element<FP>::Type ScalarMultivariateFunctionType;

  public:
    /** @brief Return the domain. */
    const DomainType domain() const override {
        return static_cast<FP const&>(*this).domain();
    }
    /** @brief Return the codomain. */
    const CodomainType codomain() const override {
        return static_cast<FP const&>(*this).codomain();
    }

    /** @brief Return the error bound. */
    const ErrorType _error() const override {
        return static_cast<ErrorType>(static_cast<const FP&>(*this).error());
    }
    /** @brief Return the component error bounds. */
    const Vector<ErrorType> _errors() const override {
        return static_cast<Vector<ErrorType>>(static_cast<const FP&>(*this).errors());
    }

    /** @brief Discard the error bound of the model. */
    inline virtual Void _clobber() override {
        static_cast<FP&>(*this).clobber();
    }
    /** @brief Test whether the first enclosure is contained in the second. */
    inline virtual Bool _refines(const FunctionPatchInterface<P, SIG>& fp) const override {
        ARIADNE_ASSERT(dynamic_cast<const FP*>(&fp));
        return refines(static_cast<const FP&>(*this), dynamic_cast<const FP&>(fp));
    }
    /** @brief Compose the outer object with the inner object. */
    inline virtual ScalarFunctionPatchInterface<P, ARGS...>* _compose(const ScalarFunction<P, RES>& f) const override {
        return heap_move(compose(f, static_cast<const FP&>(*this)));
    }
    /** @brief Compose the outer object with the inner object. */
    inline virtual VectorFunctionPatchInterface<P, ARGS...>* _compose(const VectorFunction<P, RES>& f) const override {
        return heap_move(compose(f, static_cast<const FP&>(*this)));
    }
    // FIXME: Unchecked compose dispatches to compose
    /** @brief Compose without checking domain compatibility. */
    inline virtual ScalarFunctionPatchInterface<P, ARGS...>* _unchecked_compose(const ScalarFunction<P, RES>& f) const override {
        return this->_compose(f);
    }
    /** @brief Compose without checking domain compatibility. */
    inline virtual VectorFunctionPatchInterface<P, ARGS...>* _unchecked_compose(const VectorFunction<P, RES>& f) const override {
        return this->_compose(f);
    }

    /** @brief Return a dynamically allocated copy preserving the concrete type using the concrete implementation. */
    virtual VectorFunctionPatchInterface<P, ARGS...>* _concrete_clone() const override { return new FP(static_cast<const FP&>(*this)); }
    /** @brief Replace the selected value. */
    virtual Void _set(SizeType i, const ScalarFunctionPatchInterface<P, ARGS...>& sf) override {
        if (!dynamic_cast<const typename FP::ScalarMultivariateFunctionType*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorMultivariateFunctionPatch " << *this << " to " << sf << "\n");
        }
        static_cast<FP&>(*this).FP::set(i, dynamic_cast<const ScalarMultivariateFunctionType&>(sf));
    }
    /** @brief Return the derivative with respect to the selected argument using the concrete implementation. */
    virtual VectorFunctionPatchInterface<P, ARGS...>* _concrete_derivative(SizeType) const override {
        ARIADNE_NOT_IMPLEMENTED;
    }
    /** @brief Return an antiderivative with respect to the selected argument. */
    virtual VectorFunctionPatchInterface<P, ARGS...>* _antiderivative(SizeType) const override {
        ARIADNE_NOT_IMPLEMENTED;
    }
    /** @brief Return an antiderivative with respect to the selected argument. */
    virtual VectorFunctionPatchInterface<P, ARGS...>* _antiderivative(SizeType, Number<P>) const override {
        ARIADNE_NOT_IMPLEMENTED;
    }
    /** @brief Return the range enclosure. */
    const RangeType _range() const override {
        return RangeType(static_cast<const FP&>(*this).range());
    }
    /** @brief Return the norm. */
    const NormType _norm() const override {
        return norm(static_cast<const FP&>(*this));
    }
    /** @brief Embed in the larger argument space. */
    VectorFunctionPatchInterface<P, ARGS...>* _embed(const BoxDomainType& d1, const BoxDomainType& d2) const override {
        return heap_copy(embed(d1, static_cast<const FP&>(*this), d2));
    }
    /** @brief Return the object restricted to the supplied domain. */
    VectorFunctionPatchInterface<P, ARGS...>* _restriction(const BoxDomainType& d) const override {
        return new FP(restriction(static_cast<const FP&>(*this), d));
    }
    /** @brief Add the supplied value or component. */
    Void _adjoin(const ScalarFunctionPatchInterface<P, ARGS...>& f) override {
        static_cast<FP&>(*this).FP::adjoin(dynamic_cast<const ScalarMultivariateFunctionType&>(f));
    }
    /** @brief Join the operands in their natural order. */
    VectorFunctionPatchInterface<P, ARGS...>* _join(const VectorFunctionPatchInterface<P, ARGS...>& f) const override {
        return heap_copy(join(static_cast<const FP&>(*this), dynamic_cast<const FP&>(f)));
    }
    /** @brief Combine the operands over their respective arguments. */
    VectorFunctionPatchInterface<P, ARGS...>* _combine(const VectorFunctionPatchInterface<P, ARGS...>& f) const override {
        return heap_copy(combine(static_cast<const FP&>(*this), dynamic_cast<const FP&>(f)));
    }
    /** @brief Evaluate without checking that the argument belongs to the domain. */
    Vector<Number<P>> _unchecked_evaluate(const Vector<Number<P>>& x) const override {
        return unchecked_evaluate(static_cast<const FP&>(*this), x);
    }
    /** @brief Evaluate without checking that the argument belongs to the domain. */
    Vector<CanonicalNumericType<P, DP>> _unchecked_evaluate(const Vector<CanonicalNumericType<P, DP>>& x) const override {
        auto r = unchecked_evaluate(static_cast<const FP&>(*this), x);
        if constexpr (Convertible<decltype(r), Vector<CanonicalNumericType<P, DP>>>) {
            return r;
        } else {
            ARIADNE_FAIL_MSG("");
        }
    }
    /** @brief Evaluate without checking that the argument belongs to the domain. */
    Vector<CanonicalNumericType<P, MultiplePrecision>> _unchecked_evaluate(const Vector<CanonicalNumericType<P, MultiplePrecision>>& x) const override {
        auto r = unchecked_evaluate(static_cast<const FP&>(*this), x);
        if constexpr (Convertible<decltype(r), Vector<CanonicalNumericType<P, MultiplePrecision>>>) {
            return r;
        } else {
            ARIADNE_FAIL_MSG("");
        }
    }

    /** @brief Fix the selected argument to the supplied value. */
    VectorFunctionPatchInterface<P, ARGS...>* _partial_evaluate(SizeType j, const Number<P>& c) const override {
        return heap_copy(partial_evaluate(static_cast<const FP&>(*this), j, c));
    }
};

/** @brief Reusable implementation of function patch factory. */
template<class FCTRY, class P>
class FunctionPatchFactoryMixin
    : public FunctionPatchFactoryInterface<P> {
    typedef BoxDomainType VD;
    typedef IntervalDomainType SD;

  public:
    typedef VD VectorDomainType;
    typedef SD ScalarDomainType;

    /** @brief Return a dynamically allocated copy preserving the concrete type. */
    virtual FunctionPatchFactoryInterface<P>* clone() const override { return new FCTRY(this->upcast()); }
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const override { return os << this->upcast(); }

  private:
    /** @brief Return a dynamically allocated moved value. */
    template<class T>
    static inline T* heap_move(T&& t) { return new T(std::forward<T>(t)); }
    /** @brief Convert to the common interface or value type. */
    inline FCTRY const& upcast() const { return static_cast<FCTRY const&>(*this); }

    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual ScalarUnivariateFunctionPatchInterface<P>* _create(const ScalarDomainType&, const ScalarUnivariateFunctionInterface<P>&) const override {
        ARIADNE_NOT_IMPLEMENTED;
    } // return heap_move(this->upcast().create(domain,function)); };
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual VectorUnivariateFunctionPatchInterface<P>* _create(const ScalarDomainType&, const VectorUnivariateFunctionInterface<P>&) const override {
        ARIADNE_NOT_IMPLEMENTED;
    } // return heap_move(this->upcast().create(domain,function)); };

    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual ScalarMultivariateFunctionPatchInterface<P>* _create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const override {
        return heap_move(this->upcast().create(domain, function));
    };
    /** @brief Construct a compatible object from the supplied value and characteristics. */
    virtual VectorMultivariateFunctionPatchInterface<P>* _create(const VectorDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const override {
        return heap_move(this->upcast().create(domain, function));
    };

    /** @brief Construct a zero object with compatible characteristics. */
    virtual ScalarUnivariateFunctionPatchInterface<P>* _create_zero(const ScalarDomainType&) const override {
        ARIADNE_NOT_IMPLEMENTED;
    } // return heap_move(this->upcast().create_zero(domain)); };
    /** @brief Construct a constant object with the supplied value. */
    virtual ScalarUnivariateFunctionPatchInterface<P>* _create_constant(const ScalarDomainType&, const Number<P>&) const override {
        ARIADNE_NOT_IMPLEMENTED;
    } // return heap_move(this->upcast().create_constant(domain,value)); };
    /** @brief Construct zero components with compatible characteristics. */
    virtual VectorUnivariateFunctionPatchInterface<P>* _create_zeros(SizeType, const ScalarDomainType&) const override {
        ARIADNE_NOT_IMPLEMENTED;
    } // return heap_move(this->upcast().create_zeros(rsize,domain)); };
    /** @brief Construct constant components with the supplied values. */
    virtual VectorUnivariateFunctionPatchInterface<P>* _create_constants(const ScalarDomainType&, const Vector<Number<P>>&) const override {
        ARIADNE_NOT_IMPLEMENTED;
    } // return heap_move(this->upcast().create_constants(domain,values)); };
    /** @brief Construct the identity function. */
    virtual ScalarUnivariateFunctionPatchInterface<P>* _create_identity(const ScalarDomainType&) const override {
        ARIADNE_NOT_IMPLEMENTED;
    } // return heap_move(this->upcast().create_identity(domain)); };

    /** @brief Construct a zero object with compatible characteristics. */
    virtual ScalarMultivariateFunctionPatchInterface<P>* _create_zero(const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_zero(domain));
    };
    /** @brief Construct a constant object with the supplied value. */
    virtual ScalarMultivariateFunctionPatchInterface<P>* _create_constant(const VectorDomainType& domain, const Number<P>& value) const override {
        return heap_move(this->upcast().create_constant(domain, value));
    };
    /** @brief Construct the function returning the selected coordinate. */
    virtual ScalarMultivariateFunctionPatchInterface<P>* _create_coordinate(const VectorDomainType& domain, SizeType index) const override {
        return heap_move(this->upcast().create_coordinate(domain, index));
    };
    /** @brief Construct zero components with compatible characteristics. */
    virtual VectorMultivariateFunctionPatchInterface<P>* _create_zeros(SizeType rsize, const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_zeros(rsize, domain));
    };
    /** @brief Construct constant components with the supplied values. */
    virtual VectorMultivariateFunctionPatchInterface<P>* _create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const override {
        return heap_move(this->upcast().create_constants(domain, values));
    };
    /** @brief Construct the function selecting the supplied coordinates. */
    virtual VectorMultivariateFunctionPatchInterface<P>* _create_projection(const VectorDomainType& domain, Range indices) const override {
        return heap_move(this->upcast().create_projection(domain, indices));
    };
    /** @brief Construct the identity function. */
    virtual VectorMultivariateFunctionPatchInterface<P>* _create_identity(const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_identity(domain));
    };
};

} // namespace Ariadne

#endif
