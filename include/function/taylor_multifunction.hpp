/***************************************************************************
 *            function/taylor_multifunction.hpp
 *
 *  Copyright 2012-21  Pieter Collins
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

/*! @file function/taylor_multifunction.hpp
 *  \brief Multivalued functions based on Taylor polynomials
 */

#ifndef ARIADNE_TAYLOR_MULTIFUNCTION_HPP
#define ARIADNE_TAYLOR_MULTIFUNCTION_HPP

#include "algebra/algebra.hpp"
#include "algebra/differential.hpp"
#include "algebra/evaluate.hpp"
#include "algebra/matrix.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/polynomial.hpp"
#include "algebra/vector.hpp"
#include "geometry/box.hpp"
#include "geometry/interval_set.hpp"
#include "geometry/set_wrapper.hpp"
#include "interval/interval.hpp"
#include "multifunction.hpp"
#include "numeric/numeric.hpp"
#include "scaled_function_patch.hpp"
#include "taylor_model.hpp"

namespace Ariadne {

inline ValidatedLowerKleenean operator<(const UpperBound<FloatMP>& x1, const LowerBound<FloatDP>& x2) {
    if (x1.raw() < x2.raw()) {
        return true;
    } else {
        return ValidatedLowerKleenean(ValidatedKleenean(LogicalValue::INDETERMINATE));
    }
}
inline ValidatedLowerKleenean operator>(const LowerBound<FloatMP>& x1, const UpperBound<FloatDP>& x2) {
    if (x1.raw() > x2.raw()) {
        return true;
    } else {
        return ValidatedLowerKleenean(ValidatedKleenean(LogicalValue::INDETERMINATE));
    }
}

template<class M>
class ScaledFunctionPatchMixin;

template<class F>
class ScaledFunctionPatchFactory<ValidatedIntervalTaylorModel<F>> {
    using M = ValidatedIntervalTaylorModel<F>;

  public:
    using Paradigm = typename M::Paradigm;
    using PrecisionType = typename M::PrecisionType;
    using ErrorPrecisionType = typename M::ErrorPrecisionType;
    using PropertiesType = typename M::PropertiesType;

    explicit ScaledFunctionPatchFactory(const PropertiesType& properties) : _properties(properties) {}
    PropertiesType properties() const { return _properties; }

  private:
    PropertiesType _properties;
};

template<class F>
class ScaledFunctionPatchMixin<ValidatedIntervalTaylorModel<F>>
    : public MultifunctionInterface<ValidatedTag, Real(RealVector), CompactSet> {
    using M = ValidatedIntervalTaylorModel<F>;
    using P = ValidatedTag;
    using RES = Real;
    const ScaledFunctionPatch<M>& _upcast() const { return static_cast<const ScaledFunctionPatch<M>&>(*this); }

  protected:
    virtual CompactSet<P, RES> _call(const Argument<Number<P>>& x) const override;
    // The concrete scaled patch supplies its formatter.
    virtual OutputStream& _write(OutputStream&) const override = 0;
};

template<class F>
class VectorScaledFunctionPatchMixin<ValidatedIntervalTaylorModel<F>>
    : public MultifunctionInterface<ValidatedTag, RealVector(RealVector), CompactSet> {
    using M = ValidatedIntervalTaylorModel<F>;
    using P = ValidatedTag;
    using RES = RealVector;
    const VectorScaledFunctionPatch<M>& _upcast() const { return static_cast<const VectorScaledFunctionPatch<M>&>(*this); }

  protected:
    virtual CompactSet<P, RES> _call(const Argument<Number<P>>& x) const override;
    // The concrete scaled patch supplies its formatter.
    virtual OutputStream& _write(OutputStream&) const override = 0;
};

template<class F>
class ValidatedIntervalTaylorFunctionModel
    : public ScaledFunctionPatch<ValidatedIntervalTaylorModel<F>> {
    using P = ValidatedTag;
    using M = ValidatedIntervalTaylorModel<F>;
    using RES = Real;

  public:
    ValidatedIntervalTaylorFunctionModel(const ScaledFunctionPatch<M>& f) : ScaledFunctionPatch<M>(f) {}
    using ScaledFunctionPatch<M>::ScaledFunctionPatch;
    Interval<UpperBound<F>> operator()(const Vector<Bounds<F>>& x) const {
        return evaluate(this->model(), Ariadne::unscale(x, this->domain()));
    }
    CompactSet<P, RES> operator()(const Vector<Number<P>>& x) const { return this->_call(x); }
};

template<class F>
class ValidatedVectorIntervalTaylorFunctionModel
    : public VectorScaledFunctionPatch<ValidatedIntervalTaylorModel<F>> {
    using P = ValidatedTag;
    using M = ValidatedIntervalTaylorModel<F>;
    using RES = RealVector;

  public:
    ValidatedVectorIntervalTaylorFunctionModel(const VectorScaledFunctionPatch<M>& f) : VectorScaledFunctionPatch<M>(f) {}
    using VectorScaledFunctionPatch<M>::VectorScaledFunctionPatch;
    Box<Interval<UpperBound<F>>> operator()(const Vector<Bounds<F>>& x) const {
        return evaluate(this->models(), Ariadne::unscale(x, this->domain()));
    }
    CompactSet<P, RES> operator()(const Vector<Number<P>>& x) const { return this->_call(x); }
};

template<class F>
class ValidatedIntervalTaylorFunctionSet
    //    : public FunctionSet<ValidatedTag,Real(RealVector),CompactSet>::Interface
    : public CompactSet<ValidatedTag, Real(RealVector)>::Interface {
    using P = ValidatedTag;
    using M = ValidatedIntervalTaylorModel<F>;
    using RES = Real;
    ValidatedIntervalTaylorFunctionModel<F> _function_model;

  public:
    ValidatedIntervalTaylorFunctionSet(const ScaledFunctionPatch<M>& f) : _function_model(f) {}
    CompactSet<P, RES> operator()(const Vector<Number<P>>& x) const { return this->_function_model(x); }
    friend OutputStream& operator<<(OutputStream& os, const ValidatedIntervalTaylorFunctionSet<F>& fs) {
        return os << "FunctionSet(" << fs._function_model << ")";
    }

  private:
    virtual CompactSet<P, RES> _call(const Vector<Number<P>>& x) const { return this->operator()(x); }
    virtual OutputStream& _write(OutputStream& os) const { return os << *this; }
};

template<class F>
class ValidatedVectorIntervalTaylorFunctionSet
    //    : public FunctionSet<ValidatedTag,Real(RealVector),CompactSet>::Interface
    : public CompactSet<ValidatedTag, RealVector(RealVector)>::Interface {
    using P = ValidatedTag;
    using M = ValidatedIntervalTaylorModel<F>;
    using RES = RealVector;
    ValidatedVectorIntervalTaylorFunctionModel<F> _function_model;

  public:
    ValidatedVectorIntervalTaylorFunctionSet(const VectorScaledFunctionPatch<M>& f) : _function_model(f) {}
    CompactSet<P, RES> operator()(const Vector<Number<P>>& x) const { return this->_function_model(x); }
    friend OutputStream& operator<<(OutputStream& os, const ValidatedVectorIntervalTaylorFunctionSet<F>& fs) {
        return os << "FunctionSet(" << fs._function_model << ")";
    }

  private:
    virtual CompactSet<P, RES> _call(const Vector<Number<P>>& x) const { return this->operator()(x); }
    virtual OutputStream& _write(OutputStream& os) const { return os << *this; }
};

template<class F>
class ValidatedIntervalTaylorCurriedTrajectorySet
    : public CompactSet<ValidatedTag, RealVector(Real)>::Interface {
    using P = ValidatedTag;
    using M = ValidatedIntervalTaylorModel<F>;
    using RES = RealVector;
    ValidatedVectorIntervalTaylorFunctionModel<F> _flow_model;
    Vector<Number<P>> _initial_state;

  public:
    ValidatedIntervalTaylorCurriedTrajectorySet(const ValidatedVectorIntervalTaylorFunctionModel<F>& phi, Vector<Number<P>> x0)
        : _flow_model(phi), _initial_state(x0) {}
    CompactSet<P, RES> operator()(const Scalar<Number<P>>& t) const { return this->_flow_model(join(this->_initial_state, t)); }
    friend OutputStream& operator<<(OutputStream& os, const ValidatedIntervalTaylorCurriedTrajectorySet<F>& trs) {
        return os << "TrajectorySet(" << trs._flow_model << ", " << trs._initial_state << ")";
    }

  private:
    virtual CompactSet<P, RES> _call(const Scalar<Number<P>>& t) const override {
        return this->operator()(t);
    }
    virtual OutputStream& _write(OutputStream& os) const override { return os << *this; }
};

template<class F>
class ValidatedIntervalTaylorMultiflow
    : public ValidatedCompactMultiflow::Interface {
    using P = ValidatedTag;
    using M = ValidatedIntervalTaylorModel<F>;
    using RES = Real;
    typedef ValidatedVectorIntervalTaylorFunctionModel<F> FlowModelType;
    typedef ValidatedIntervalTaylorCurriedTrajectorySet<F> TrajectorySetType;
    typedef ValidatedVectorIntervalTaylorFunctionModel<F> ReachMapType;
    FlowModelType _flow_model;

  public:
    ValidatedIntervalTaylorMultiflow(const VectorScaledFunctionPatch<M>& f) : _flow_model(f) {}
    TrajectorySetType operator()(const Vector<Number<P>>& x0) const {
        return TrajectorySetType(this->_flow_model, x0);
    }
    ReachMapType curry() const {
        return ReachMapType(this->_flow_model);
    }
    friend OutputStream& operator<<(OutputStream& os, const ValidatedIntervalTaylorMultiflow<F>& phi) { return phi._write(os); }

  private:

  public:
    virtual CompactSet<P, RealVector(Real)> _call(const Vector<Number<P>>& x0) const override {
        return CompactSet<P, RealVector(Real)>(std::make_shared<TrajectorySetType>(this->_flow_model, x0));
    }
    virtual Function<P, CompactSet<P, RealVector>(RealVector)> _curry() const override {
        return Function<P, CompactSet<P, RealVector>(RealVector)>(std::make_shared<ReachMapType>(this->curry()));
    }
    virtual OutputStream& _write(OutputStream& os) const override { return os << "Multiflow(" << this->_flow_model << ")"; }
};

} // namespace Ariadne

#endif

