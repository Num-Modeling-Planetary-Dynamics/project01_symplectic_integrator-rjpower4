#include <gtest/gtest.h>

#include "eaps/celestial_body.hpp"

TEST(CelestialBodyTest, PositiveGmsShouldNotThrow)
{
    EXPECT_NO_THROW(eaps::CelestialBody("earth", 399, 1.3));
}

TEST(CelestialBodyTest, NegativeGmConstructionShouldThrow)
{
    EXPECT_THROW(eaps::CelestialBody("earth", 399, -1.3), std::invalid_argument);
}

TEST(CelestialBodyTest, ZeroGmConstructionShouldThrow)
{
    EXPECT_THROW(eaps::CelestialBody("earth", 399, 0.0), std::invalid_argument);
}

TEST(CelestialBodyTest, BasicConstructionShouldAssignValues)
{
    auto cb = eaps::CelestialBody("earth", 399, 1.1);
    EXPECT_EQ(cb.Name(), "earth");
    EXPECT_EQ(cb.NaifId(), 399);
    EXPECT_EQ(cb.Gm(), 1.1);
}