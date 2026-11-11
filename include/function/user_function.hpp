/***************************************************************************
 *            function/user_function.hpp
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

/*! @file function/user_function.hpp
 *  \brief Wrappers for user functions
 */

#ifndef ARIADNE_USER_FUNCTION_HPP
#define ARIADNE_USER_FUNCTION_HPP

#include <cstdarg>
#include <iosfwd>

#include "algebra/differential.hpp"
#include "algebra/matrix.hpp"
#include "algebra/vector.hpp"
#include "function/formula.hpp"
#include "function/function_interface.hpp"
#include "function/function_mixin.hpp"
#include "function/taylor_model.hpp"
#include "function_mixin.hpp"
#include "numeric/numeric.hpp"
#include "utility/container.hpp"
#include "utility/macros.hpp"
#include "utility/pointer.hpp"

namespace Ariadne {

//! \related ScalarUserFunction
//! \brief Template structure containing scalar function meta-information.
template<SizeType AS, SizeType PS = 0u, DegreeType SM = 255u>
struct ScalarMultivariateFunctionData {
    /** @brief Return the number of scalar arguments. */
    static SizeType argument_size() { return AS; }
    /** @brief Return the number of parameters. */
    static SizeType parameter_size() { return PS; }
    /** @brief Return the supported differentiation order. */
    static DegreeType smoothness() { return SM; }
};

//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//!
//! Given a C++ class F with a (static) template method
//!   <code>template<class R, class A, class P> compute(R& r, const A& a, const P& p);</code>
//! the type <code>ScalarUserFunction<F></code> is an Ariadne function defined by \f$r=f(a)\f$.
//! The constructor for ScalarUserFunction<F> takes a Vector<ApproximateNumericType> or Vector<ValidatedNumericType> argument which is used for \a p.
//!
//! The class F must also define meta-data <c>argument_size(), parameter_size()
//! and smoothness()</c>. These are most easily defined by inheriting from the
//! <tt>ScalarMultivariateFunctionData<AS,PS,SM=SMOOTH></tt> class.
//!
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
template<class F>
class ScalarUserFunction
    : public EffectiveScalarMultivariateFunction {
  private:
    /** @brief Non-owning wrapper selecting representation output for a value. */
    class Representation
        : public ScalarMultivariateFunctionMixin<Representation, EffectiveTag> {
      private:
        Vector<EffectiveNumber> _p;

      public:
        /** @brief Construct a \c Representation from p. */
        Representation(const Vector<EffectiveNumber>& p) : _p(p) {}

        /** @brief Compute the value at the requested accuracy or effort. */
        template<class R, class A>
        inline Void _compute(R& r, const A& a) const { F::compute(r, a, _p); }

        /** @brief Return a dynamically allocated copy preserving the concrete type. */
        virtual Representation* clone() const { return new Representation(*this); }

        /** @brief Return the number of scalar arguments. */
        virtual SizeType argument_size() const { return F::argument_size(); }
        /** @brief Return the number of parameters. */
        virtual SizeType parameter_size() const { return F::parameter_size(); }

        /** @brief Return the derivative with respect to the selected argument. */
        virtual EffectiveScalarMultivariateFunction derivative(SizeType) const { ARIADNE_NOT_IMPLEMENTED; }

        /** @brief Return the vector of first partial derivatives. */
        virtual Covector<FloatDPApproximation> gradient(const Vector<FloatDPApproximation>& x) const {
            return this->evaluate(Differential<FloatDPApproximation>::variables(1u, x)).gradient();
        }
        /** @brief Return the vector of first partial derivatives. */
        virtual Covector<FloatDPBounds> gradient(const Vector<FloatDPBounds>& x) const {
            return this->evaluate(Differential<FloatDPBounds>::variables(1u, x)).gradient();
        }

        /** @brief Write the representation of the object. */
        virtual OutputStream& repr(OutputStream& os) const {
            return os << "USER";
        }
        /** @brief Write the value to the output stream. */
        virtual OutputStream& _write(OutputStream& os) const {
            return os << "ScalarUserFunction( argument_size=" << this->argument_size() << " )";
        }
    };

  public:
    /** @brief Construct a \c ScalarUserFunction with the specified characteristics. */
    ScalarUserFunction() : EffectiveScalarMultivariateFunction(new Representation(Vector<EffectiveNumber>(this->parameter_size()))) {}
    /** @brief Construct a \c ScalarUserFunction from p. */
    ScalarUserFunction(const Vector<Real>& p) : EffectiveScalarMultivariateFunction(new Representation(Vector<EffectiveNumber>(p))) {}
    /** @brief Construct a \c ScalarUserFunction from p. */
    ScalarUserFunction(const Vector<EffectiveNumber>& p) : EffectiveScalarMultivariateFunction(new Representation(p)) {}

    /** @brief Return the number of parameters. */
    SizeType parameter_size() const { return F().parameter_size(); }
    // const Vector<ValidatedNumericType> parameters() const { return _p; }
};

//! \brief A convenience class defining the meta-data of an %Ariadne function.
template<SizeType RS, SizeType AS, SizeType PS = 0u, DegreeType SM = 255u>
class VectorMultivariateFunctionData {
  public:
    //!
    /** @brief Return the number of result components. */
    static SizeType result_size() { return RS; }
    //!
    /** @brief Return the number of scalar arguments. */
    static SizeType argument_size() { return AS; }
    //!
    /** @brief Return the number of parameters. */
    static SizeType parameter_size() { return PS; }
    //!
    /** @brief Return the supported differentiation order. */
    static DegreeType smoothness() { return SM; }
};

//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//!
//! Given a C++ class F with a (static) template method
//!   <code>template<class R, class A, class P> compute(R& r, const A& a, const P& p);</code>
//! the type <code>VectorUserFunction<F></code> is an Ariadne function defined by \f$r=f(a)\f$.
//! The constructor for VectorUserFunction<F> takes a Vector<EffectiveNumericType> argument which is used for \a p.
//!
//! The class F must also define meta-data <c>result_size(), argument_size(), parameter_size()
//! and smoothness()</c> as static functions. These are most easily defined by inheriting from the
//! <tt>VectorMultivariateFunctionData<RS,AS,PS,SM=SMOOTH></tt> class.
//!
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
template<class F>
class VectorUserFunction
    : public EffectiveVectorMultivariateFunction {
    typedef DoublePrecision PR;

  private:
    /** @brief Non-owning wrapper selecting representation output for a value. */
    class Representation
        : public VectorMultivariateFunctionMixin<Representation, EffectiveTag> {
      public:
        /** @brief Construct a \c Representation from p. */
        Representation(const Vector<EffectiveNumber>& p) : _p(p) {}

        /** @brief Return a dynamically allocated copy preserving the concrete type. */
        virtual Representation* clone() const { return new Representation(*this); }

        /** @brief Return the number of result components. */
        virtual SizeType result_size() const { return F::result_size(); }
        /** @brief Return the number of scalar arguments. */
        virtual SizeType argument_size() const { return F::argument_size(); }
        /** @brief Return the number of parameters. */
        virtual SizeType parameter_size() const { return F::parameter_size(); }

        /** @brief Compute the value at the requested accuracy or effort. */
        template<class R, class A>
        inline Void _compute(R& r, const A& a) const { F::compute(r, a, _p); }

        /** @brief Return the matrix of first partial derivatives. */
        virtual Matrix<FloatDPApproximation> jacobian(const Vector<FloatDPApproximation>& x) const {
            return Ariadne::jacobian(this->evaluate(Differential<FloatDPApproximation>::variables(1u, x)));
        }
        /** @brief Return the matrix of first partial derivatives. */
        virtual Matrix<FloatDPBounds> jacobian(const Vector<FloatDPBounds>& x) const {
            return Ariadne::jacobian(this->evaluate(Differential<FloatDPBounds>::variables(1u, x)));
        }

        /** @brief Return the requested element or stored value. */
        virtual EffectiveScalarMultivariateFunction::Interface* _get(SizeType) const { ARIADNE_NOT_IMPLEMENTED; }
        /** @brief Return the selected element or component. */
        virtual EffectiveScalarMultivariateFunction operator[](SizeType) const { ARIADNE_NOT_IMPLEMENTED; }

        // TODO: Find a better way for writing functions which can handle transformations which may not have a
        // _write() method or operator<<.
        /** @brief Write the value to the output stream. */
        virtual OutputStream& _write(OutputStream& os) const {
            return os << "VectorUserFunction( result_size=" << this->result_size() << ", argument_size=" << this->argument_size() << " )";
        }

        Vector<EffectiveNumber> _p;
    };

  public:
    /** @brief Construct a \c VectorUserFunction from p. */
    VectorUserFunction(const Vector<EffectiveNumber>& p) : EffectiveVectorMultivariateFunction(new Representation(p)) {}
    /** @brief Return the stored parameters. */
    const Vector<EffectiveNumber>& parameters() const { return dynamic_cast<const Representation*>(this->raw_pointer())->_p; }
};

} // namespace Ariadne

#endif /* ARIADNE_USER_FUNCTION_HPP */
