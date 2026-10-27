#include "geometry/box.hpp"
#include "function/function.hpp"
#include "function/taylor_model.hpp"
#include "function/formula.hpp"
#include "algebra/algebra.hpp"

namespace Ariadne {
typedef Interval<FloatDPUpperBound> UpperIntervalType;

namespace {

Vector<FloatDPBounds> singleton_bounds(Vector<UpperIntervalType> const& x) {
    Vector<FloatDPBounds> r(x.size(),FloatDPBounds(dp));
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=cast_singleton(x[i]); }
    return r;
}

UpperBoxType upper_box(Vector<FloatDPBounds> const& x) {
    UpperBoxType r(x.size());
    for(SizeType i=0; i!=x.size(); ++i) { r[i]=make_interval(x[i]); }
    return r;
}

} // namespace

UpperIntervalType apply(ValidatedScalarUnivariateFunction const& f, UpperIntervalType const& ivl) {
    if (definitely(ivl.is_empty())) { return UpperIntervalType::empty_interval(); }
    return make_interval(f(cast_singleton(ivl))); }
UpperBoxType apply(ValidatedVectorUnivariateFunction const& f, UpperIntervalType const& ivl) {
    if (definitely(ivl.is_empty())) { return UpperBoxType(f.result_size(),UpperIntervalType::empty_interval()); }
    return upper_box(f(cast_singleton(ivl))); }
UpperIntervalType apply(ValidatedScalarMultivariateFunction const& f, UpperBoxType const& bx) {
    if (definitely(bx.is_empty())) { return UpperIntervalType::empty_interval(); }
    return make_interval(f(singleton_bounds(cast_vector(bx)))); }
UpperBoxType apply(ValidatedVectorMultivariateFunction const& f, UpperBoxType const& bx) {
    if (definitely(bx.is_empty())) { return UpperBoxType(f.result_size(),UpperIntervalType::empty_interval()); }
    return upper_box(f(singleton_bounds(cast_vector(bx)))); }

UpperIntervalType image(UpperIntervalType const& ivl, ValidatedScalarUnivariateFunction const& f) {
    return apply(f,ivl); }
UpperBoxType image(UpperIntervalType const& ivl, ValidatedVectorUnivariateFunction const& f) {
    return apply(f,ivl); }
UpperIntervalType image(UpperBoxType const& bx, ValidatedScalarMultivariateFunction const& f) {
    return apply(f,bx); }
UpperBoxType image(UpperBoxType const& bx, ValidatedVectorMultivariateFunction const& f) {
    return apply(f,bx); }


} // namespace Ariadne
