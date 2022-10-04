#include <gtest/gtest.h>

#include "eaps/twobody/keplerian_elements.hpp"

#include <limits>

using eaps::twobody::KeplerianElements;

TEST(KeplerianElementsTest, BasicConstruction) {
    KeplerianElements elements{1.1, 0.2, 3.3, 4.4, 5.5, 6.6};
    EXPECT_EQ(elements.SemiMajorAxis(), 1.1);
    EXPECT_EQ(elements.Eccentricity(), 0.2);
    EXPECT_EQ(elements.Inclination(), 3.3);
    EXPECT_EQ(elements.ArgumentOfPeriapsis(), 4.4);
    EXPECT_EQ(elements.RightAscension(), 5.5);
    EXPECT_EQ(elements.TrueAnomaly(), 6.6);
}

TEST(KeplerianElementsTest, InvalidConstruction) {
    double infinity = std::numeric_limits<double>::infinity();

    EXPECT_THROW(KeplerianElements(10.3, 1.5), std::invalid_argument);
    EXPECT_THROW(KeplerianElements(0.0, 1.5), std::invalid_argument);
    EXPECT_THROW(KeplerianElements(0.0, 0.5), std::invalid_argument);
    EXPECT_THROW(KeplerianElements(-10.3, 0.5), std::invalid_argument);
    EXPECT_THROW(KeplerianElements(0.0, 0.0), std::invalid_argument);
    EXPECT_THROW(KeplerianElements(0.0, 1.0), std::invalid_argument);
    EXPECT_THROW(KeplerianElements(1.1, 1.0), std::invalid_argument);
    EXPECT_THROW(KeplerianElements(-1.1, 1.0), std::invalid_argument);
    EXPECT_THROW(KeplerianElements( infinity, 0.3), std::invalid_argument);
    EXPECT_THROW(KeplerianElements(12.3, 1.0), std::invalid_argument);
    EXPECT_NO_THROW(KeplerianElements(infinity, 1.0));
}