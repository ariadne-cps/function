/***************************************************************************
 *            function/taylor_function.decl.hpp
 *
 *  Forward declarations for Taylor function factories.
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

/** @file function/taylor_function.decl.hpp
 *  @brief Declarations of Taylor function types.
 */

#ifndef ARIADNE_TAYLOR_FUNCTION_DECL_HPP
#define ARIADNE_TAYLOR_FUNCTION_DECL_HPP

#include "function/declarations.hpp"

namespace Ariadne {

template<class F>
class Sweeper;
template<class P>
class FunctionPatchFactoryInterface;

/** @brief Construct taylor function patch factory from the supplied arguments. */
FunctionPatchFactoryInterface<ValidatedTag>* make_taylor_function_patch_factory();
/** @brief Construct taylor function patch factory from the supplied arguments. */
FunctionPatchFactoryInterface<ValidatedTag>* make_taylor_function_patch_factory(const Sweeper<FloatDP>& sweeper);

} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_DECL_HPP
