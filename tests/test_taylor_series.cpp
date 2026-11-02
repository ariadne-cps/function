/***************************************************************************
 *            test_taylor_series.cpp
 *
 *  This file is part of Ariadne, distributed under the GNU GPL version 3
 *  or (at your option) any later version.
 ****************************************************************************/
#include "numeric/numeric.hpp"
#include "algebra/series.hpp"
#include "function/taylor_series.hpp"
#include "ariadne-testing.hpp"
#include <sstream>
using namespace Ariadne;

class TestTaylorSeries : public Ariadne::Testing::ContextFixture {
public:
    using ContextFixture::ContextFixture;
Void test_polynomial_series() {
    IntervalDomainType domain(-1,1);
    TaylorSeries<FloatDPBounds> empty_series(domain,3u);
    ARIADNE_TESTING_EXPECT_EQUAL(empty_series.degree(), 3u);
    ARIADNE_TESTING_EXPECT_EQUAL(empty_series.expansion().size(), 4u);
    ARIADNE_TESTING_EXPECT_EQUAL(empty_series.error().raw(), FloatDP(0,dp));
    TaylorSeries<FloatDPBounds> square(Sqr(),domain,FloatDP(0,dp),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(square.degree(), 2u);
    ARIADNE_TESTING_EXPECT_EQUAL(square[0u], FloatDP(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(square[1u], FloatDP(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(square[2u], FloatDP(1,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(square.error().raw(), FloatDP(0,dp));
    square.sweep(FloatDP(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(square[2u], FloatDP(1,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(square.error().raw(), FloatDP(0,dp));
    square.sweep(FloatDP(1,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(square[2u], FloatDP(0,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(square.error().raw(), FloatDP(1,dp));
    TaylorSeries<FloatDPBounds> constant(Sqr(),IntervalDomainType(2,2),FloatDP(2,dp),0u);
    ARIADNE_TESTING_EXPECT_EQUAL(constant[0u], FloatDP(4,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(constant.error().raw(), FloatDP(0,dp));
    std::ostringstream text;text<<square;
    ARIADNE_TESTING_PROVE(text.str().find("TS(")==0u);
}
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_taylor_series");
    ARIADNE_TESTING_ADD_TEST(TestTaylorSeries, test_polynomial_series);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
