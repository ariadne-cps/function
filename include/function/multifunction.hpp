/***************************************************************************
 *            function/multifunction.hpp
 *
 *  Copyright 2008-21  Pieter Collins
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Function. If not, see <https://www.gnu.org/licenses/>.
 */

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! @file function/multifunction.hpp
 *  \brief Multivalued functions
 */

#ifndef ARIADNE_MULTIFUNCTION_HPP
#define ARIADNE_MULTIFUNCTION_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "algebra/covector.hpp"
#include "algebra/differential.hpp"
#include "algebra/vector.hpp"
#include "function/declarations.hpp"
#include "function/function.hpp"
#include "numeric/numeric.hpp"
#include "utility/container.hpp"
#include "utility/macros.hpp"
#include "utility/metaprogramming.hpp"
#include "utility/pointer.hpp"
// #include "function/function_model.hpp"

#include "geometry/box.hpp"
#include "geometry/set.hpp"
#include "interval/interval.hpp"

namespace Ariadne {

/** @brief Polymorphic representation of a point in a generic space. */
template<class P>
class GenericPoint {
    /** @brief Polymorphic interface for . */
    class Interface {
        /** @brief Destroy the object and release its owned resources. */
        virtual ~Interface() = 0;
        /** @brief Return the lower and upper bounds. */
        virtual Vector<Interval<Number<P>>> _bounds() const;
        /** @brief Apply the supplied operation to the operands. */
        GenericPoint<P> _apply(const VectorFunction<P>& f) const;
    };
    SharedPointer<const Interface> _ptr;

  public:
    /** @brief Construct a \c GenericPoint from ptr. */
    GenericPoint(SharedPointer<const Interface> ptr) : _ptr(ptr) {}
    /** @brief Return the lower and upper bounds. */
    Vector<Interval<Number<P>>> bounds() { return this->_ptr->_bounds(); }
    /** @brief Return the selected element or component. */
    Interval<Number<P>> operator[](SizeType i) { return this->bounds()[i]; }
    /** @brief Apply the supplied operation to the operands. */
    friend GenericPoint<P> apply(const VectorFunction<P>& f, const GenericPoint<P>& pt) { return pt._ptr->_apply(f); }
};

/** @brief Generic point represented as the image of another point. */
template<class P>
class ImageGenericPoint : public GenericPoint<P>::Interface {
    VectorFunction<P> _h;
    /** @brief Construct a \c ImageGenericPoint from h. */
    ImageGenericPoint(const VectorFunction<P>& h) : _h(h) {}
    /** @brief Apply the supplied operation to the operands. */
    friend ImageGenericPoint<P> apply(const VectorFunction<P>& f, const GenericPoint<P>& pt) {
        return ImageGenericPoint<P>(f(pt._h));
    }
};

template<class P, class SIG, template<class, class> class SET = LocatedSet>
class Multifunction;

template<class P>
using VectorMultivariateMultifunction = Multifunction<P, RealVector(RealVector)>;
using ValidatedVectorMultivariateMultifunction = VectorMultivariateMultifunction<ValidatedTag>;

template<class P, class SIG, template<class, class> class SET = LocatedSet>
class MultifunctionInterface;

/** @brief Polymorphic interface for multifunction. */
template<class P, class SIG, template<class, class> class SET>
class MultifunctionInterface {
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    using RES = typename SignatureTraits<SIG>::ResultKind;

  public:
    template<class Y>
    using Argument = typename SignatureTraits<SIG>::template Argument<Y>;
    /** @brief Destroy the object and release its owned resources. */
    virtual ~MultifunctionInterface() = default;
    /** @brief Evaluate with the supplied arguments. */
    virtual SET<P, RES> _call(const Argument<Number<P>>& x) const = 0;
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

/** @brief Compile-time predicate testing whether a type is multifunction. */
template<class F, class P, class SIG, template<class, class> class SET>
struct IsMultifunction {
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    using RES = typename SignatureTraits<SIG>::ResultKind;
    template<class Y>
    using Argument = typename SignatureTraits<SIG>::template Argument<Y>;

    /** @brief Evaluate the compile-time or runtime predicate for the supplied arguments. */
    template<class FF, class = decltype(declval<SET<P, RES>>() = declval<FF>()(declval<Argument<Number<P>>>()))>
    static std::true_type test(int);
    /** @brief Evaluate the compile-time or runtime predicate for the supplied arguments. */
    template<class FF>
    static std::false_type test(...);
    static const bool value = decltype(test<F>(1))::value;
};

template<class F, class P, class SIG, template<class, class> class SET>
concept AMultifunction = IsMultifunction<F, P, SIG, SET>::value;

//! \ingroup Function
//! \brief Functions \f$\X\mvto\Y\f$ returning a \ref LocatedSet at each point.
template<class P, class SIG, template<class, class> class SET>
class Multifunction
    : public Handle<const MultifunctionInterface<P, SIG, SET>> {
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    using RES = typename SignatureTraits<SIG>::ResultKind;
    using ResultSetType = SET<P, RES>;

  public:
    typedef MultifunctionInterface<P, SIG, SET> Interface;
    template<class Y>
    using Argument = typename SignatureTraits<SIG>::template Argument<Y>;

    using Handle<const Interface>::Handle;
    //! \brief <p/>
    /** @brief Construct a \c Multifunction from mf. */
    template<AMultifunction<P, SIG, SET> MF>
    explicit Multifunction(MF const& mf);
    //! \brief <p/>
    /** @brief Apply the callable to the supplied arguments. */
    SET<P, RES> operator()(const Argument<Number<P>>& x) const { return this->reference()._call(x); }
    //! \brief <p/>
    /** @brief Apply the supplied operation to the operands. */
    friend SET<P, RES> apply(const Multifunction<P, SIG, SET>& f, const SET<P, ARG>& x);
    //! \brief <p/>
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Multifunction<P, SIG, SET>& mf) { return mf.reference()._write(os); }
};

/** @brief Adapter exposing multifunction. */
template<class MF, class P, class SIG, template<class, class> class SET>
class MultifunctionWrapper
    : public virtual MultifunctionInterface<P, SIG, SET> {
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    using RES = typename SignatureTraits<SIG>::ResultKind;

  public:
    template<class Y>
    using Argument = typename SignatureTraits<SIG>::template Argument<Y>;

    /** @brief Construct a \c MultifunctionWrapper from mf. */
    MultifunctionWrapper(MF mf) : _mf(mf) {}
    /** @brief Return the underlying base object. */
    MF const& base() const { return this->_mf; }
    /** @brief Evaluate with the supplied arguments. */
    virtual SET<P, RES> _call(const Argument<Number<P>>& x) const final override {
        return static_cast<SET<P, RES>>(this->_mf(x));
    }
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const final override {
        return os << this->base();
    }

  private:
    MF _mf;
};

template<class P, class SIG, template<class, class> class SET>
template<AMultifunction<P, SIG, SET> MF>
Multifunction<P, SIG, SET>::Multifunction(MF const& mf)
    : Multifunction<P, SIG, SET>(std::make_shared<const MultifunctionWrapper<MF, P, SIG, SET>>(mf)) {
}

/** @brief Function representation parameterized by its computational paradigm and signature. */
template<class P, class RES, class ARG>
class Function<P, CompactSet<P, RES>(ARG)>
    : public Multifunction<P, RES(ARG), CompactSet> {
    using Multifunction<P, RES(ARG), CompactSet>::Multifunction;
};
/** @brief Function representation parameterized by its computational paradigm and signature. */
template<class P, class RES, class ARG>
class Function<P, LocatedSet<P, RES>(ARG)>
    : public Multifunction<P, RES(ARG), LocatedSet> {
    using Multifunction<P, RES(ARG), LocatedSet>::Multifunction;
};

} // namespace Ariadne

#include "geometry/set_wrapper.hpp"

namespace Ariadne {

template<class P, class SIG, template<class, class> class SET = LocatedSet>
class FunctionSet;

/** @brief Compile-time predicate testing whether a type is function set. */
template<class F, class P, class SIG, template<class, class> class SET>
struct IsFunctionSet {
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    using RES = typename SignatureTraits<SIG>::ResultKind;
    template<class Y>
    using Argument = typename SignatureTraits<SIG>::template Argument<Y>;

    /** @brief Evaluate the compile-time or runtime predicate for the supplied arguments. */
    template<class FF, class = decltype(declval<SET<P, RES>>() = declval<FF>()(declval<Argument<Number<P>>>()))>
    static std::true_type test(int);
    /** @brief Evaluate the compile-time or runtime predicate for the supplied arguments. */
    template<class FF>
    static std::false_type test(...);
    static const bool value = decltype(test<F>(1))::value;
};

template<class F, class P, class SIG, template<class, class> class SET>
concept AFunctionSet = IsFunctionSet<F, P, SIG, SET>::value;

/** @brief Polymorphic interface for function set. */
template<class P, class SIG, template<class, class> class SET>
class FunctionSetInterface {
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    using RES = typename SignatureTraits<SIG>::ResultKind;

  public:
    template<class Y>
    using Argument = typename SignatureTraits<SIG>::template Argument<Y>;
    /** @brief Destroy the object and release its owned resources. */
    virtual ~FunctionSetInterface() = default;
    /** @brief Evaluate with the supplied arguments. */
    virtual SET<P, RES> _call(const Argument<Number<P>>& x) const = 0;
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

/** @brief Representation of a set of functions. */
template<class P, class SIG, template<class, class> class SET>
class FunctionSet
    : public Handle<const FunctionSetInterface<P, SIG, SET>> {
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    using RES = typename SignatureTraits<SIG>::ResultKind;

  public:
    typedef FunctionSetInterface<P, SIG, SET> Interface;
    template<class Y>
    using Argument = typename SignatureTraits<SIG>::template Argument<Y>;

    using Handle<const Interface>::Handle;

    //! \brief <p/>
    /** @brief Construct a \c FunctionSet from fs. */
    template<AFunctionSet<P, SIG, SET> FS>
    explicit FunctionSet(FS const& fs);
    //! \brief <p/>
    /** @brief Apply the callable to the supplied arguments. */
    SET<P, RES> operator()(const Argument<Number<P>>& x) const { return this->reference()._call(x); }
    //! \brief <p/>
    /** @brief Apply the supplied operation to the operands. */
    friend SET<P, RES> apply(const FunctionSet<P, SIG, SET>& fs, const SET<P, ARG>& x);
    //! \brief <p/>
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const FunctionSet<P, SIG, SET>& fs) { return fs.reference()._write(os); }
};

/** @brief Adapter exposing function set. */
template<class FS, class P, class SIG, template<class, class> class SET>
class FunctionSetWrapper
    : public virtual FunctionSetInterface<P, SIG, SET> {
    using ARG = typename SignatureTraits<SIG>::ArgumentKind;
    using RES = typename SignatureTraits<SIG>::ResultKind;

  public:
    template<class X>
    using Argument = typename SignatureTraits<SIG>::template Argument<X>;

    /** @brief Construct a \c FunctionSetWrapper from mf. */
    FunctionSetWrapper(FS mf) : _fs(mf) {}
    /** @brief Return the underlying base object. */
    FS const& base() const { return this->_fs; }
    /** @brief Evaluate with the supplied arguments. */
    virtual SET<P, RES> _call(const Argument<Number<P>>& x) const final override {
        return static_cast<SET<P, RES>>(this->_fs(x));
    }
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const final override {
        return os << this->base();
    }

  private:
    FS _fs;
};

template<class P, class SIG, template<class, class> class SET>
template<AFunctionSet<P, SIG, SET> FS>
FunctionSet<P, SIG, SET>::FunctionSet(FS const& fs)
    : Handle<const Interface>(std::make_shared<const FunctionSetWrapper<FS, P, SIG, SET>>(fs)) {}

/** @brief Representation of a compact set with paradigm-specific predicates. */
template<class RES, class ARG>
class CompactSet<ValidatedTag, RES(ARG)>
    : public FunctionSet<ValidatedTag, RES(ARG), CompactSet> {
    using Base = FunctionSet<ValidatedTag, RES(ARG), CompactSet>;

  public:
    using Base::Base;
    /** @brief Construct a \c CompactSet from fs. */
    template<AFunctionSet<ValidatedTag, RES(ARG), CompactSet> FS>
    explicit CompactSet(FS const& fs) : Base(fs) {}
};
/** @brief Representation of a set supporting location predicates. */
template<class RES, class ARG>
class LocatedSet<ValidatedTag, RES(ARG)>
    : public FunctionSet<ValidatedTag, RES(ARG), LocatedSet> {
    using Base = FunctionSet<ValidatedTag, RES(ARG), LocatedSet>;

  public:
    using Base::Base;
    /** @brief Construct a \c LocatedSet from fs. */
    template<AFunctionSet<ValidatedTag, RES(ARG), LocatedSet> FS>
    explicit LocatedSet(FS const& fs) : Base(fs) {}
};

template<class A1, class A2>
using JoinType = decltype(join(declval<A1>(), declval<A2>()));

/** @brief Function whose values are sets of functions. */
template<class P, class RES, class PAR, class ARG, template<class, class> class SET>
class FunctionSetFunction {
    typedef JoinType<ARG, PAR> ARGS;

  public:
    template<class Y>
    using Argument = typename DomainTraits<ARG>::template Type<Y>;
    typedef SET<P, RES(PAR)> FunctionSetType;
    /** @brief Polymorphic interface for . */
    class Interface {
      private:

      public:
        /** @brief Evaluate with the supplied arguments. */
        virtual SET<P, RES(PAR)> _call(const Argument<Number<P>>&) const = 0;
        /** @brief Fix an argument group and return a function of the remaining arguments. */
        virtual Function<P, SET<P, RES>(ARGS)> _curry() const = 0;
        /** @brief Write the value to the output stream. */
        virtual OutputStream& _write(OutputStream&) const = 0;
    };
    /** @brief Apply the callable to the supplied arguments. */
    SET<P, RES(PAR)> operator()(const Argument<Number<P>>& x) const { return this->_handle._ptr->_call(x); }
    /** @brief Fix an argument group and return a function of the remaining arguments. */
    Function<P, SET<P, RES>(ARGS)> curry() const { return this->_handle._ptr->_curry(); }

  private:
    Handle<Interface> _handle;
};

/** @brief Function representation parameterized by its computational paradigm and signature. */
template<class P, class RES, class PAR, class ARG>
class Function<P, CompactSet<P, RES(PAR)>(ARG)>
    : public FunctionSetFunction<P, RES, PAR, ARG, CompactSet> {};

template<class P, template<class, class> class SET = LocatedSet>
using Multiflow = Function<P, SET<P, RealVector(Real)>(RealVector)>;
template<class P>
using CompactMultiflow = Function<P, CompactSet<P, RealVector(Real)>(RealVector)>;
template<class P>
using LocatedMultiflow = Multiflow<P, LocatedSet>;

using ValidatedCompactMultiflow = Multiflow<ValidatedTag, CompactSet>;

} // namespace Ariadne

#endif

