/***************************************************************************
 *            box_function.cpp
 *
 *  Copyright  2026  Luca Geretti
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

#include "algebra/algebra.hpp"
#include "function/formula.hpp"
#include "function/function.hpp"
#include "function/taylor_model.hpp"
#include "geometry/box.hpp"

namespace Ariadne {
typedef Interval<FloatDPUpperBound> UpperIntervalType;

namespace {

Vector<FloatDPBounds> singleton_bounds(const Vector<UpperIntervalType>& x) {
    Vector<FloatDPBounds> r(x.size(), FloatDPBounds(dp));
    for (SizeType i = 0; i != x.size(); ++i) {
        r[i] = cast_singleton(x[i]);
    }
    return r;
}

UpperBoxType upper_box(const Vector<FloatDPBounds>& x) {
    UpperBoxType r(x.size());
    for (SizeType i = 0; i != x.size(); ++i) {
        r[i] = make_interval(x[i]);
    }
    return r;
}

} // namespace

UpperIntervalType apply(const ValidatedScalarUnivariateFunction& f, const UpperIntervalType& ivl) {
    if (definitely(ivl.is_empty())) {
        return UpperIntervalType::empty_interval();
    }
    return make_interval(f(cast_singleton(ivl)));
}
UpperBoxType apply(const ValidatedVectorUnivariateFunction& f, const UpperIntervalType& ivl) {
    if (definitely(ivl.is_empty())) {
        return UpperBoxType(f.result_size(), UpperIntervalType::empty_interval());
    }
    return upper_box(f(cast_singleton(ivl)));
}
UpperIntervalType apply(const ValidatedScalarMultivariateFunction& f, const UpperBoxType& bx) {
    if (definitely(bx.is_empty())) {
        return UpperIntervalType::empty_interval();
    }
    return make_interval(f(singleton_bounds(cast_vector(bx))));
}
UpperBoxType apply(const ValidatedVectorMultivariateFunction& f, const UpperBoxType& bx) {
    if (definitely(bx.is_empty())) {
        return UpperBoxType(f.result_size(), UpperIntervalType::empty_interval());
    }
    return upper_box(f(singleton_bounds(cast_vector(bx))));
}

UpperIntervalType image(const UpperIntervalType& ivl, const ValidatedScalarUnivariateFunction& f) {
    return apply(f, ivl);
}
UpperBoxType image(const UpperIntervalType& ivl, const ValidatedVectorUnivariateFunction& f) {
    return apply(f, ivl);
}
UpperIntervalType image(const UpperBoxType& bx, const ValidatedScalarMultivariateFunction& f) {
    return apply(f, bx);
}
UpperBoxType image(const UpperBoxType& bx, const ValidatedVectorMultivariateFunction& f) {
    return apply(f, bx);
}

} // namespace Ariadne
