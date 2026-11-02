/***************************************************************************
 * This file is part of Ariadne, distributed under the GNU GPL version 3
 * or (at your option) any later version.
 ****************************************************************************/
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
    ARIADNE_TESTING_EXPECT_EQUAL(projection(point), Vector<FloatDPApproximation>({3,1,3},dp));
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
