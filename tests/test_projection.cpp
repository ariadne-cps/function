/***************************************************************************
 * test_projection.cpp
 * or (at your option) any later version.
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

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "function/projection.hpp"
#include "ariadne-testing.hpp"
#include <sstream>
using namespace Ariadne;

class TestProjection : public Ariadne::Testing::ContextFixture {
public:
    using ContextFixture::ContextFixture;
Void test_projection() {
    const Projection projection(3u,Array<SizeType>{2u,0u,2u});
    ARIADNE_TESTING_EXPECT_EQUAL(projection.argument_size(), 3u);
    ARIADNE_TESTING_EXPECT_EQUAL(projection.result_size(), 3u);
    ARIADNE_TESTING_EXPECT_EQUAL(projection.indices(), Array<SizeType>({2u,0u,2u}));
    ARIADNE_TESTING_EXPECT_EQUAL(projection[1u], 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(projection.index(2u), 2u);
    const Vector<FloatDPApproximation> point({1,2,3},dp);
    ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (projection(point)); return snapshot.array(); }()), ([&] { auto snapshot = (Vector<FloatDPApproximation>({3,1,3},dp)); return snapshot.array(); }()), [](auto const& left, auto const& right) { return same(left,right); });
    const Projection empty(0u,Array<SizeType>{});
    ARIADNE_TESTING_EXPECT_EQUAL(empty.argument_size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(empty.result_size(), 0u);
    ARIADNE_TESTING_EXPECT_EQUAL(empty(Vector<FloatDPApproximation>(0u,dp)).size(), 0u);
    ARIADNE_TESTING_EXPECT_EXCEPTION(Projection(3u,Array<SizeType>({0u,3u})));
    std::ostringstream text;text<<projection;
    ARIADNE_TESTING_EXPECT_EQUAL(text.str(), String("Projection([2,0,2])"));
}
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_projection");
    ARIADNE_TESTING_ADD_TEST(TestProjection, test_projection);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
