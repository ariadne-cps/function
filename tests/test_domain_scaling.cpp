/***************************************************************************
 *            test_domain_scaling.cpp
 *
 *  Copyright  2026  Luca Geretti
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

#include <sstream>

#include "function/domain.hpp"
#include "function/scaling.hpp"
#include "numeric/numeric.hpp"

#include "ariadne-testing.hpp"
using namespace Ariadne;

class TestDomainScaling : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    Void test_domains() {
        RealDomain real, explicit_real{SizeOne()};
        ARIADNE_TESTING_EXPECT_EQUAL(real.dimension(), 1u);
        ARIADNE_TESTING_PROVE(real == explicit_real);
        ARIADNE_TESTING_PROVE(intersection(real, explicit_real) == real);
        IntervalDomainType entire = real;
        ARIADNE_TESTING_EXPECT_EQUAL(entire.lower_bound(), FloatDP(-infty, dp));
        ARIADNE_TESTING_EXPECT_EQUAL(entire.upper_bound(), FloatDP(+infty, dp));
        ARIADNE_TESTING_PROVE(real == entire);
        ARIADNE_TESTING_PROVE(entire == real);
        ARIADNE_TESTING_PROVE(!(real == IntervalDomainType(-1, 1)));
        ARIADNE_TESTING_PROVE(!(IntervalDomainType(-1, 1) == real));
        EuclideanDomain plane(2u), another_plane(2u, real), line(1u);
        ARIADNE_TESTING_EXPECT_EQUAL(plane.dimension(), 2u);
        ARIADNE_TESTING_PROVE(plane[0u] == real);
        ARIADNE_TESTING_PROVE(plane == another_plane);
        ARIADNE_TESTING_PROVE(!(plane == line));
        ARIADNE_TESTING_PROVE(intersection(plane, another_plane) == plane);
        ARIADNE_TESTING_EXPECT_EXCEPTION(intersection(plane, line));
        ARIADNE_TESTING_EXPECT_EQUAL(product(plane, line).dimension(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(product(plane, real).dimension(), 3u);
        BoxDomainType entire_plane = plane;
        ARIADNE_TESTING_EXPECT_EQUAL(entire_plane.dimension(), 2u);
        ARIADNE_TESTING_PROVE(plane == entire_plane);
        ARIADNE_TESTING_PROVE(entire_plane == plane);
        ARIADNE_TESTING_PROVE(!(plane == BoxDomainType::unit_box(2u)));
        ARIADNE_TESTING_PROVE(!(BoxDomainType::unit_box(2u) == plane));
        UnitInterval unit;
        ARIADNE_TESTING_EXPECT_EQUAL(unit.dimension(), 1u);
        ARIADNE_TESTING_PROVE(unit == UnitInterval());
        ARIADNE_TESTING_PROVE(intersection(unit, unit) == unit);
        ARIADNE_TESTING_EXPECT_EQUAL(IntervalDomainType(unit), IntervalDomainType(-1, 1));
        // Both conversion APIs must describe the same signed unit interval.
        ARIADNE_TESTING_EXPECT_EQUAL(unit.operator IntervalDomainType(), IntervalDomainType(unit));
        UnitBox square(2u), another_square(2u, unit), segment(1u);
        ARIADNE_TESTING_EXPECT_EQUAL(square.dimension(), 2u);
        ARIADNE_TESTING_PROVE(square[0u] == unit);
        ARIADNE_TESTING_PROVE(square == another_square);
        ARIADNE_TESTING_PROVE(!(square == segment));
        ARIADNE_TESTING_PROVE(intersection(square, another_square) == square);
        ARIADNE_TESTING_EXPECT_EXCEPTION(intersection(square, segment));
        ARIADNE_TESTING_EXPECT_EQUAL(product(square, segment).dimension(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(product(square, unit).dimension(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(BoxDomainType(square), BoxDomainType::unit_box(2u));
        ARIADNE_TESTING_EXPECT_EQUAL(BoxDomainType(UnitBox(0u)).dimension(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(BoxDomainType(EuclideanDomain(0u)).dimension(), 0u);
        std::ostringstream text;
        text << real << " " << plane << " " << unit << " " << square;
        ARIADNE_TESTING_EXPECT_EQUAL(text.str(), String("R R2 [-1:+1] [-1:+1]^2"));
    }

    Void test_scaling() {
        IntervalDomainType interval(2, 6);
        ARIADNE_TESTING_EXPECT_EQUAL(med(interval), Dyadic(4));
        ARIADNE_TESTING_EXPECT_EQUAL(rad(interval), Dyadic(2));
        ARIADNE_TESTING_EXPECT_SAME(med_apprx(interval), FloatDPApproximation(4, dp));
        ARIADNE_TESTING_EXPECT_SAME(rad_apprx(interval), FloatDPApproximation(2, dp));
        ARIADNE_TESTING_PROVE(definitely(med_val(interval) == FloatDPBounds(4, dp)));
        ARIADNE_TESTING_PROVE(definitely(rad_val(interval) == FloatDPBounds(2, dp)));
        ARIADNE_TESTING_EXPECT_EQUAL(convert_interval(interval, dp), interval);
        MultiplePrecision precision(128u);
        auto mp_interval = convert_interval(interval, precision);
        ARIADNE_TESTING_EXPECT_EQUAL(mp_interval.lower_bound(), FloatMP(2, precision));
        ARIADNE_TESTING_EXPECT_EQUAL(mp_interval.upper_bound(), FloatMP(6, precision));
        Scaling scaling(interval);
        Unscaling unscaling(interval);
        ARIADNE_TESTING_PROVE(scaling.domain() == UnitInterval());
        ARIADNE_TESTING_EXPECT_EQUAL(scaling.codomain(), interval);
        ARIADNE_TESTING_EXPECT_EQUAL(unscaling.domain(), interval);
        ARIADNE_TESTING_PROVE(unscaling.codomain() == UnitInterval());
        for (Int k = -1; k <= 1; ++k) {
            FloatDPApproximation x(k, dp);
            ARIADNE_TESTING_EXPECT_SAME(scaling(x), FloatDPApproximation(4 + 2 * k, dp));
            ARIADNE_TESTING_EXPECT_SAME(unscaling(scaling(x)), x);
        }
        ARIADNE_TESTING_EXPECT_SAME(unscale(FloatDPApproximation(3, dp), IntervalDomainType(3, 3)), FloatDPApproximation(0, dp));
        BoxDomainType box{
            { 2,  6},
            {-3, -1}
        };
        VectorScaling vector_scaling(box);
        VectorUnscaling vector_unscaling(box);
        ARIADNE_TESTING_EXPECT_EQUAL(vector_scaling.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(vector_unscaling.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(vector_scaling.codomain(), box);
        ARIADNE_TESTING_EXPECT_EQUAL(vector_unscaling.domain(), box);
        ARIADNE_TESTING_EXPECT_EQUAL(vector_scaling[0u].codomain(), box[0u]);
        ARIADNE_TESTING_EXPECT_EQUAL(vector_unscaling[1u].domain(), box[1u]);
        Vector<FloatDPApproximation> point({-1, 1}, dp), expected({2, -1}, dp);
        ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (vector_scaling(point)); return snapshot.array(); }()), ([&] { auto snapshot = (expected); return snapshot.array(); }()), [](const auto& left, const auto& right) { return same(left, right); });
        ARIADNE_TESTING_EXPECT_RANGES_BY(([&] { auto snapshot = (vector_unscaling(expected)); return snapshot.array(); }()), ([&] { auto snapshot = (point); return snapshot.array(); }()), [](const auto& left, const auto& right) { return same(left, right); });
        Vector<FloatDPApproximation> wrong_dimension({0}, dp);
        ARIADNE_TESTING_EXPECT_EXCEPTION(vector_scaling(wrong_dimension));
        ARIADNE_TESTING_EXPECT_EXCEPTION(vector_unscaling(wrong_dimension));
        Vector<FloatDPApproximation> empty_point(0u, dp);
        ARIADNE_TESTING_EXPECT_EQUAL(VectorScaling(BoxDomainType())(empty_point).size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(VectorUnscaling(BoxDomainType())(empty_point).size(), 0u);
        UnscalingException error("invalid domain ", interval);
        ARIADNE_TESTING_EXPECT_EQUAL(error.domain, interval);
        ARIADNE_TESTING_PROVE(String(error.what()).find("invalid domain") != String::npos);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_domain_scaling");
    ARIADNE_TESTING_ADD_TEST(TestDomainScaling, test_domains);
    ARIADNE_TESTING_ADD_TEST(TestDomainScaling, test_scaling);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
