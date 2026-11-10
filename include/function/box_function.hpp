/***************************************************************************
 *            function/box_function.hpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

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
