#include <gtest/gtest.h>

#include "eaps/vector3.hpp"

using eaps::Vector3;

TEST(Vector3Test, SimpleConstruction)
{
    Vector3 vector{2, 3, 4};

    EXPECT_EQ(vector.X(), 2.0);
    EXPECT_EQ(vector.Y(), 3.0);
    EXPECT_EQ(vector.Z(), 4.0);
}
