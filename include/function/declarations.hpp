/***************************************************************************
 *            function/declarations.hpp
 *
 *  Copyright  2011-26  Pieter Collins
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


#ifndef ARIADNE_FUNCTION_DECLARATIONS_HPP
#define ARIADNE_FUNCTION_DECLARATIONS_HPP

#include <iosfwd>

#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"

#include "foundation/paradigm.hpp"
#include "foundation/logical.decl.hpp"

#include "function/function.decl.hpp"

namespace Ariadne {


template<class P, class F> class AffineModel;
template<class P, class F> class TaylorModel;

} // namespace Ariadne

#endif /* ARIADNE_FUNCTION_DECLARATIONS_HPP */
