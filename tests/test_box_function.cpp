/***************************************************************************
 *            test_box_function.cpp
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Function. If not, see <https://www.gnu.org/licenses/>.
 */

// Coverage of validated Function images over intervals and boxes,
// including empty-domain dimensionality.
#include "geometry/box.hpp"
#include "function/function.hpp"
#include "function/box_function.hpp"
#include "function/domain.hpp"
#include "ariadne-testing.hpp"
using namespace Ariadne;

class TestBoxFunction : public Ariadne::Testing::ContextFixture {
public:
    using ContextFixture::ContextFixture;
Void test_validated_function_images() {
        // Exercise scalar and vector images over both inhabited and empty
        // domains.  Empty images must retain their result dimension.
        ValidatedScalarUnivariateFunction scalar_univariate=
            ValidatedScalarUnivariateFunction::coordinate();
        ValidatedVectorUnivariateFunction vector_univariate=
            ValidatedVectorUnivariateFunction::zeros(2u);
        ValidatedScalarMultivariateFunction scalar_multivariate=
            ValidatedScalarMultivariateFunction::coordinate(2u,0u);
        ValidatedVectorMultivariateFunction vector_multivariate=
            ValidatedVectorMultivariateFunction::identity(2u);

        UpperIntervalType segment(0.0_x,1.0_x);
        UpperIntervalType empty_segment=UpperIntervalType::empty_interval();
        UpperBoxType rectangle({{0.0_x,1.0_x},{2.0_x,3.0_x}});
        UpperBoxType empty_rectangle(2u,empty_segment);

        auto scalar_segment_image=image(segment,scalar_univariate);
        ARIADNE_TESTING_INSPECT("scalar_segment_image", scalar_segment_image);
        ARIADNE_TESTING_PROVE(definitely(subset(scalar_segment_image,ExactIntervalType(-1.0_x,2.0_x))));
        auto vector_segment_image=image(segment,vector_univariate);
        ARIADNE_TESTING_INSPECT("vector_segment_image", vector_segment_image);
        ARIADNE_TESTING_EXPECT_EQUAL(vector_segment_image.dimension(), 2u);
        ARIADNE_TESTING_PROVE(definitely(subset(vector_segment_image,ExactBoxType({{-1.0_x,1.0_x},{-1.0_x,1.0_x}}))));

        auto scalar_rectangle_image=image(rectangle,scalar_multivariate);
        ARIADNE_TESTING_INSPECT("scalar_rectangle_image", scalar_rectangle_image);
        ARIADNE_TESTING_PROVE(definitely(subset(scalar_rectangle_image,ExactIntervalType(-1.0_x,2.0_x))));
        auto vector_rectangle_image=image(rectangle,vector_multivariate);
        ARIADNE_TESTING_INSPECT("vector_rectangle_image", vector_rectangle_image);
        ARIADNE_TESTING_EXPECT_EQUAL(vector_rectangle_image.dimension(), 2u);
        ARIADNE_TESTING_PROVE(definitely(subset(vector_rectangle_image,ExactBoxType({{-1.0_x,2.0_x},{1.0_x,4.0_x}}))));

        auto empty_scalar_segment=image(empty_segment,scalar_univariate);
        ARIADNE_TESTING_PROVE(definitely(empty_scalar_segment.is_empty()));
        auto empty_vector_segment=image(empty_segment,vector_univariate);
        ARIADNE_TESTING_EXPECT_EQUAL(empty_vector_segment.dimension(), 2u);
        ARIADNE_TESTING_PROVE(definitely(empty_vector_segment.is_empty()));

        auto empty_scalar_rectangle=image(empty_rectangle,scalar_multivariate);
        ARIADNE_TESTING_PROVE(definitely(empty_scalar_rectangle.is_empty()));
        auto empty_vector_rectangle=image(empty_rectangle,vector_multivariate);
        ARIADNE_TESTING_EXPECT_EQUAL(empty_vector_rectangle.dimension(), 2u);
        ARIADNE_TESTING_PROVE(definitely(empty_vector_rectangle.is_empty()));
    }


Void test_euclidean_domain_product() {
    ARIADNE_TESTING_EXPECT_EQUAL(product(EuclideanDomain(1u),EuclideanDomain(2u),EuclideanDomain(3u)).dimension(), 6u);
}
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_box_function");
    ARIADNE_TESTING_ADD_TEST(TestBoxFunction, test_validated_function_images);
    ARIADNE_TESTING_ADD_TEST(TestBoxFunction, test_euclidean_domain_product);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
