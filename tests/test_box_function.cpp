// Coverage of validated Function images over intervals and boxes,
// including empty-domain dimensionality.
#include "geometry/box.hpp"
#include "function/function.hpp"
#include "function/box_function.hpp"
#include "function/domain.hpp"
#include "macro/testing.hpp"
using namespace Ariadne;

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
        ARIADNE_TEST_PRINT(scalar_segment_image);
        ARIADNE_TEST_ASSERT(definitely(subset(scalar_segment_image,ExactIntervalType(-1.0_x,2.0_x))));
        auto vector_segment_image=image(segment,vector_univariate);
        ARIADNE_TEST_PRINT(vector_segment_image);
        ARIADNE_TEST_EQUALS(vector_segment_image.dimension(),2u);
        ARIADNE_TEST_ASSERT(definitely(subset(vector_segment_image,ExactBoxType({{-1.0_x,1.0_x},{-1.0_x,1.0_x}}))));

        auto scalar_rectangle_image=image(rectangle,scalar_multivariate);
        ARIADNE_TEST_PRINT(scalar_rectangle_image);
        ARIADNE_TEST_ASSERT(definitely(subset(scalar_rectangle_image,ExactIntervalType(-1.0_x,2.0_x))));
        auto vector_rectangle_image=image(rectangle,vector_multivariate);
        ARIADNE_TEST_PRINT(vector_rectangle_image);
        ARIADNE_TEST_EQUALS(vector_rectangle_image.dimension(),2u);
        ARIADNE_TEST_ASSERT(definitely(subset(vector_rectangle_image,ExactBoxType({{-1.0_x,2.0_x},{1.0_x,4.0_x}}))));

        auto empty_scalar_segment=image(empty_segment,scalar_univariate);
        ARIADNE_TEST_ASSERT(definitely(empty_scalar_segment.is_empty()));
        auto empty_vector_segment=image(empty_segment,vector_univariate);
        ARIADNE_TEST_EQUALS(empty_vector_segment.dimension(),2u);
        ARIADNE_TEST_ASSERT(definitely(empty_vector_segment.is_empty()));

        auto empty_scalar_rectangle=image(empty_rectangle,scalar_multivariate);
        ARIADNE_TEST_ASSERT(definitely(empty_scalar_rectangle.is_empty()));
        auto empty_vector_rectangle=image(empty_rectangle,vector_multivariate);
        ARIADNE_TEST_EQUALS(empty_vector_rectangle.dimension(),2u);
        ARIADNE_TEST_ASSERT(definitely(empty_vector_rectangle.is_empty()));
    }


Void test_euclidean_domain_product() {
    ARIADNE_TEST_EQUALS(product(EuclideanDomain(1u),EuclideanDomain(2u),EuclideanDomain(3u)).dimension(),6u);
}
int main() {
    ARIADNE_TEST_CALL(test_validated_function_images());
    ARIADNE_TEST_CALL(test_euclidean_domain_product());
    return ARIADNE_TEST_FAILURES;
}
