#include <gtest/gtest.h>

#include "eaps/twobody/conic.hpp"

using eaps::twobody::ConicOrbitTypeFromEccentricity;
using eaps::twobody::ConicOrbitType;

TEST(ConicOrbitTypeFromEccentricityTest, Basic) {
    EXPECT_EQ(ConicOrbitTypeFromEccentricity(0.0), ConicOrbitType::Circular);
    EXPECT_EQ(ConicOrbitTypeFromEccentricity(0.4), ConicOrbitType::Elliptical);
    EXPECT_EQ(ConicOrbitTypeFromEccentricity(1.0), ConicOrbitType::Parabolic);
    EXPECT_EQ(ConicOrbitTypeFromEccentricity(12.0), ConicOrbitType::Hyperbolic);
}

TEST(ConicOrbitTypeFromEccentricityTest, ThrowsOnBadEccentricity) {
    EXPECT_THROW(ConicOrbitTypeFromEccentricity(-1.0), std::invalid_argument);
}