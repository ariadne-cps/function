/***************************************************************************
 *            function/domain.hpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
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

/*! @file function/domain.hpp
 *  \brief Interval and box domains for functions.
 */

#ifndef ARIADNE_DOMAIN_HPP
#define ARIADNE_DOMAIN_HPP

#include "function/function_traits.hpp"
#include "geometry/box.hpp"
#include "interval/interval.hpp"

namespace Ariadne {

//! \ingroup FunctionModule
//! \brief The type used as a bounded domain for a univariate function.
using IntervalDomainType = Interval<FloatDP>;
//! \ingroup FunctionModule
//! \brief The type used as a bounded domain for a multivariate function.
using BoxDomainType = Box<Interval<FloatDP>>;

//! \ingroup FunctionModule \brief The domain of an entire univariate function.
class RealDomain {
  public:
    typedef SizeOne DimensionType; //!< .
    typedef IndexZero IndexType;   //!< .

    constexpr RealDomain() {}                                                                                                   //!< .
    constexpr RealDomain(SizeOne) {}                                                                                            //!< .
    constexpr SizeOne dimension() const { return SizeOne(); }                                                                   //!< .
    operator IntervalDomainType() const { return IntervalDomainType(-inf, +inf); }                                              //!< .
    friend RealDomain intersection(const RealDomain &, const RealDomain &) { return RealDomain(); }                             //!< .
    friend Bool operator==(const RealDomain &, const RealDomain &) { return true; }                                             //!< .
    friend Bool operator==(const RealDomain &dom1, const IntervalDomainType &dom2) { return IntervalDomainType(dom1) == dom2; } //!< .
    friend Bool operator==(const IntervalDomainType &dom1, const RealDomain &dom2) { return dom1 == IntervalDomainType(dom2); } //!< .
    friend OutputStream &operator<<(OutputStream &os, const RealDomain &) { return os << "R"; }                                 //!< .
};

//! \ingroup FunctionModule \brief The domain of an entire multivariate function.
class EuclideanDomain {
    SizeType _dim;

  public:
    typedef SizeType DimensionType; //!< .
    typedef SizeType IndexType;     //!< .

    constexpr EuclideanDomain(SizeType dim) : _dim(dim) {}                                                        //!< .
    constexpr EuclideanDomain(SizeType dim, RealDomain) : _dim(dim) {}                                            //!< .
    constexpr SizeType dimension() const { return this->_dim; }                                                   //!< .
    constexpr RealDomain operator[](SizeType) { return RealDomain(); }                                            //!< .
    operator BoxDomainType() const { return BoxDomainType(this->dimension(), IntervalDomainType(RealDomain())); } //!< .
    friend EuclideanDomain intersection(const EuclideanDomain &dom1, const EuclideanDomain &dom2) {
        ARIADNE_PRECONDITION(dom1 == dom2);
        return dom1;
    }                                                                                                                                                         //!< .
    friend EuclideanDomain product(const EuclideanDomain &dom1, const EuclideanDomain &dom2) { return EuclideanDomain(dom1.dimension() + dom2.dimension()); } //!< .
    friend EuclideanDomain product(const EuclideanDomain &dom1, const RealDomain &dom2) { return EuclideanDomain(dom1.dimension() + dom2.dimension()); }      //!< .
    friend Bool operator==(const EuclideanDomain &dom1, const EuclideanDomain &dom2) { return dom1.dimension() == dom2.dimension(); }                         //!< .
    friend Bool operator==(const EuclideanDomain &dom1, const BoxDomainType &dom2) { return BoxDomainType(dom1) == dom2; }                                    //!< .
    friend Bool operator==(const BoxDomainType &dom1, const EuclideanDomain &dom2) { return dom1 == BoxDomainType(dom2); }                                    //!< .
    friend OutputStream &operator<<(OutputStream &os, const EuclideanDomain &dom) { return os << "R" << dom.dimension(); }                                    //!< .
};

//! \ingroup FunctionModule \brief The signed unit interval \f$[-1:+1]\f$.
class UnitInterval {
  public:
    constexpr UnitInterval() {}                                                                             //!< .
    constexpr SizeOne dimension() const { return SizeOne(); }                                               //!< .
    operator IntervalDomainType() const { return IntervalDomainType(-1, +1); }                              //!< .
    friend UnitInterval intersection(const UnitInterval &, const UnitInterval &) { return UnitInterval(); } //!< .
    friend Bool operator==(const UnitInterval &, const UnitInterval &) { return true; }                     //!< .
    friend OutputStream &operator<<(OutputStream &os, const UnitInterval &) { return os << "[-1:+1]"; }     //!< .
};

//! \ingroup FunctionModule \brief The signed unit interval \f$[-1:+1]^n\f$.
class UnitBox {
    SizeType _dim;

  public:
    constexpr UnitBox(SizeType dim) : _dim(dim) {}                                                                  //!< .
    constexpr UnitBox(SizeType dim, UnitInterval) : _dim(dim) {}                                                    //!< .
    constexpr SizeType dimension() const { return this->_dim; }                                                     //!< .
    constexpr UnitInterval operator[](SizeType) { return UnitInterval(); }                                          //!< .
    operator BoxDomainType() const { return BoxDomainType(this->dimension(), IntervalDomainType(UnitInterval())); } //!< .
    friend UnitBox intersection(const UnitBox &dom1, const UnitBox &dom2) {
        ARIADNE_PRECONDITION(dom1 == dom2);
        return dom1;
    }                                                                                                                              //!< .
    friend UnitBox product(const UnitBox &dom1, const UnitBox &dom2) { return UnitBox(dom1.dimension() + dom2.dimension()); }      //!< .
    friend UnitBox product(const UnitBox &dom1, const UnitInterval &dom2) { return UnitBox(dom1.dimension() + dom2.dimension()); } //!< .
    friend Bool operator==(const UnitBox &dom1, const UnitBox &dom2) { return dom1.dimension() == dom2.dimension(); }              //!< .
    friend OutputStream &operator<<(OutputStream &os, const UnitBox &dom) { return os << "[-1:+1]^" << dom.dimension(); }          //!< .
};

template<class S>
struct ElementTraits;
template<class UB>
struct ElementTraits<Interval<UB>> : DomainTraits<RealScalar> {};
template<class IVL>
struct ElementTraits<Box<IVL>> : DomainTraits<RealVector> {};
template<>
struct ElementTraits<RealDomain> : DomainTraits<RealScalar> {};
template<>
struct ElementTraits<EuclideanDomain> : DomainTraits<RealVector> {};
template<>
struct ElementTraits<UnitInterval> : DomainTraits<RealScalar> {};
template<>
struct ElementTraits<UnitBox> : DomainTraits<RealVector> {};

template<class... TS>
using CartesianProductType = decltype(product(declval<TS>()...));

} // namespace Ariadne

#endif /* ARIADNE_DOMAIN_HPP */
