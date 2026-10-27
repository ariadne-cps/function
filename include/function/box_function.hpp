// Validated images of intervals and boxes under Function objects.
// The definitions are in src/function/box_function.cpp.
#ifndef ARIADNE_FUNCTION_BOX_FUNCTION_HPP
#define ARIADNE_FUNCTION_BOX_FUNCTION_HPP

#include "geometry/box.decl.hpp"
#include "function/function.decl.hpp"

namespace Ariadne {

UpperIntervalType apply(ValidatedScalarUnivariateFunction const& f, UpperIntervalType const& ivl);
UpperBoxType apply(ValidatedVectorUnivariateFunction const& f, UpperIntervalType const& ivl);
UpperIntervalType apply(ValidatedScalarMultivariateFunction const& f, UpperBoxType const& bx);
UpperBoxType apply(ValidatedVectorMultivariateFunction const& f, UpperBoxType const& bx);

UpperIntervalType image(UpperIntervalType const& ivl, ValidatedScalarUnivariateFunction const& f);
UpperBoxType image(UpperIntervalType const& ivl, ValidatedVectorUnivariateFunction const& f);
UpperIntervalType image(UpperBoxType const& bx, ValidatedScalarMultivariateFunction const& f);
UpperBoxType image(UpperBoxType const& bx, ValidatedVectorMultivariateFunction const& f);

} // namespace Ariadne
#endif // ARIADNE_FUNCTION_BOX_FUNCTION_HPP
