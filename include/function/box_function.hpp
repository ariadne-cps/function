/***************************************************************************
 *            function/box_function.hpp
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Function.  If not, see <https://www.gnu.org/licenses/>.
 */

// Validated images of intervals and boxes under Function objects.
// The definitions are in src/function/box_function.cpp.

/** @file function/box_function.hpp
 *  @brief Functions defined on geometric boxes.
 */

#ifndef ARIADNE_FUNCTION_BOX_FUNCTION_HPP
#define ARIADNE_FUNCTION_BOX_FUNCTION_HPP

#include "function/function.decl.hpp"
#include "geometry/box.decl.hpp"

namespace Ariadne {

UpperIntervalType apply(const ValidatedScalarUnivariateFunction& f, const UpperIntervalType& ivl);
UpperBoxType apply(const ValidatedVectorUnivariateFunction& f, const UpperIntervalType& ivl);
UpperIntervalType apply(const ValidatedScalarMultivariateFunction& f, const UpperBoxType& bx);
UpperBoxType apply(const ValidatedVectorMultivariateFunction& f, const UpperBoxType& bx);

UpperIntervalType image(const UpperIntervalType& ivl, const ValidatedScalarUnivariateFunction& f);
UpperBoxType image(const UpperIntervalType& ivl, const ValidatedVectorUnivariateFunction& f);
UpperIntervalType image(const UpperBoxType& bx, const ValidatedScalarMultivariateFunction& f);
UpperBoxType image(const UpperBoxType& bx, const ValidatedVectorMultivariateFunction& f);

} // namespace Ariadne
#endif // ARIADNE_FUNCTION_BOX_FUNCTION_HPP
