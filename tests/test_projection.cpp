/***************************************************************************
 * This file is part of Ariadne, distributed under the GNU GPL version 3
 * or (at your option) any later version.
 ****************************************************************************/
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "function/projection.hpp"
#include "macro/testing.hpp"
#include <sstream>
using namespace Ariadne;

Void test_projection() {
    const Projection projection(3u,Array<SizeType>{2u,0u,2u});
    ARIADNE_TEST_EQUALS(projection.argument_size(),3u);
    ARIADNE_TEST_EQUALS(projection.result_size(),3u);
    ARIADNE_TEST_EQUALS(projection.indices(),Array<SizeType>({2u,0u,2u}));
    ARIADNE_TEST_EQUALS(projection[1u],0u);
    ARIADNE_TEST_EQUALS(projection.index(2u),2u);
    const Vector<FloatDPApproximation> point({1,2,3},dp);
    ARIADNE_TEST_EQUALS(projection(point),Vector<FloatDPApproximation>({3,1,3},dp));
    const Projection empty(0u,Array<SizeType>{});
    ARIADNE_TEST_EQUALS(empty.argument_size(),0u);
    ARIADNE_TEST_EQUALS(empty.result_size(),0u);
    ARIADNE_TEST_EQUALS(empty(Vector<FloatDPApproximation>(0u,dp)).size(),0u);
    ARIADNE_TEST_FAIL(Projection(3u,Array<SizeType>({0u,3u})));
    std::ostringstream text;text<<projection;
    ARIADNE_TEST_EQUALS(text.str(),String("Projection([2,0,2])"));
}
Int main() {
    ARIADNE_TEST_CALL(test_projection());
    return ARIADNE_TEST_FAILURES;
}
