/***************************************************************************
 *            function/function_wrapper.hpp
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

/** @file function/function_wrapper.hpp
 *  @brief Wrappers adapting concrete functions to abstract interfaces.
 */

#ifndef ARIADNE_FUNCTION_WRAPPER_HPP
#define ARIADNE_FUNCTION_WRAPPER_HPP

#include "function_concepts.hpp"
#include "function_mixin.hpp"

namespace Ariadne {

template<class T, class I>
class Wrapper;

template<class F, class P, class SIG>
    requires AFunction<F, P, SIG>
class FunctionWrapper;

/** @brief Adapter exposing a concrete object through an interface. */
template<class F, class P, class SIG>
struct Wrapper<F, FunctionInterface<P, SIG>> : public FunctionWrapper<F, P, SIG> {
    /** @brief Destroy the object and release its owned resources. */
    virtual ~Wrapper() = default;
    using FunctionWrapper<F, P, SIG>::FunctionWrapper;
};

template<class P, class SIG, class F>
    requires AFunction<F, P, SIG>
FunctionWrapper<F, P, SIG>* make_function_wrapper(F f);

/** @brief Reusable implementation of function wrapper getter. */
template<class F, class P, class SIG>
class FunctionWrapperGetterMixin;

/** @brief Reusable implementation of function wrapper getter. */
template<class F, class P, class... ARGS>
class FunctionWrapperGetterMixin<F, P, RealScalar(ARGS...)> {
};

/** @brief Reusable implementation of function wrapper getter. */
template<class F, class P, class... ARGS>
class FunctionWrapperGetterMixin<F, P, RealVector(ARGS...)>
    : public virtual VectorOfFunctionInterface<P, ARGS...> {
    /** @brief Return the requested element or stored value. */
    virtual ScalarFunctionInterface<P, ARGS...>* _get(SizeType i) const override;
};

/** @brief Adapter exposing function. */
template<class F, class P, class SIG>
    requires AFunction<F, P, SIG>
class FunctionWrapper
    : public FunctionMixin<FunctionWrapper<F, P, SIG>, P, SIG>,
      public FunctionWrapperGetterMixin<FunctionWrapper<F, P, SIG>, P, SIG>,
      public F {
    using typename FunctionInterface<P, SIG>::ArgumentSizeType;
    using typename FunctionInterface<P, SIG>::ResultSizeType;
    using typename FunctionInterface<P, SIG>::ArgumentIndexType;
    template<class X>
    using Argument = typename SignatureTraits<SIG>::template Argument<X>;
    template<class X>
    using Result = typename SignatureTraits<SIG>::template Result<X>;

  public:
    /** @brief Destroy the object and release its owned resources. */
    virtual ~FunctionWrapper() = default;
    /** @brief Construct a \c FunctionWrapper from f. */
    FunctionWrapper(F f) : F(f) {}
    /** @brief Convert to const f &. */
    operator const F&() const { return *this; }
    /** @brief Return the wrapped concrete object. */
    const F& wrapped() const { return *this; }
    /** @brief Return the number of scalar arguments. */
    virtual ArgumentSizeType argument_size() const override { return this->wrapped().argument_size(); }
    /** @brief Return the number of result components. */
    virtual ResultSizeType result_size() const override { return this->wrapped().result_size(); }
    /** @brief Return the derivative with respect to the selected argument. */
    virtual FunctionInterface<P, SIG>* _derivative(ArgumentIndexType k) const override {
        return make_function_wrapper<P, SIG>(derivative(this->wrapped(), k));
    }
    /** @brief Apply the callable to the supplied arguments. */
    template<class T>
    Result<T> operator()(const Argument<T>& x) const { return this->wrapped()(x); }
    /** @brief Write the value to the output stream. */
    friend inline OutputStream& operator<<(OutputStream& os, const FunctionWrapper<F, P, SIG>& f) { return os << f.wrapped(); }
    /** @brief Write the value to the output stream. */
    friend inline OutputStream& operator<<(OutputStream& os, const Representation<FunctionWrapper<F, P, SIG>>& f) { return os << Representation<F>(f.reference().wrapped()); }
    /** @brief Return the selected element or component. */
    template<class I>
    decltype(auto) operator[](I i) const { return this->wrapped()[i]; }
};

template<class P, class SIG>
template<class F>
    requires(not IsFunctionClass<F, SIG>) and AFunction<F, P, SIG>
Function<P, SIG>::Function(const F& f)
    : Handle<const Interface>(make_function_wrapper<P, SIG>(f)) {}

/** @brief Construct function wrapper from the supplied arguments. */
template<class P, class SIG, class F>
    requires AFunction<F, P, SIG>
FunctionWrapper<F, P, SIG>* make_function_wrapper(F f) { return new Wrapper<F, FunctionInterface<P, SIG>>(f); }

/** @brief Extract the requested concrete representation. */
template<class F, class P, class SIG>
const F& extract(const Function<P, SIG>& f) {
    if constexpr (DerivedFrom<F, FunctionInterface<P, SIG>>) {
        return dynamic_cast<const F&>(f.reference());
    } else {
        return dynamic_cast<const FunctionWrapper<F, P, SIG>&>(f.reference());
    }
}

} // namespace Ariadne

#include "function/function_mixin.tpl.hpp"

#endif // ARIADNE_FUNCTION_WRAPPER_HPP
