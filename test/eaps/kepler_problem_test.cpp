#include <gtest/gtest.h>

#include <array>

#include "eaps/kepler_problem.hpp"

TEST(KeplerProblemTest, ValidEccentricitiesShouldNotThrow)
{
    EXPECT_NO_THROW(eaps::KeplerProblem(0.0, 0.1));
    EXPECT_NO_THROW(eaps::KeplerProblem(0.1, 0.1));
}

TEST(KeplerProblemTest, NegativeEccentricityShouldThrow)
{
    EXPECT_THROW(eaps::KeplerProblem(-0.3, 0.2), std::invalid_argument);
}

TEST(KeplerProblemTest, ParabolicEccentricityShouldThrow)
{
    EXPECT_THROW(eaps::KeplerProblem(1.0, 0.2), std::invalid_argument);
}

TEST(KeplerProblemTest, HyperbolicEccentricityShouldThrow)
{
    EXPECT_THROW(eaps::KeplerProblem(12.0, 0.2), std::invalid_argument);
}

TEST(KeplerProblemTest, TrivialKeplerCaseShouldReturnMeanAnomaly)
{
    auto problem = eaps::KeplerProblem{0.0, 0.3};
    EXPECT_DOUBLE_EQ(SolveKepler(problem).value(), 0.3);
}

TEST(KeplerProblemTest, KeplerSolutionShouldBeRight)
{
    std::array<double, 10> eccs = {
        9.7059278176061570e-01, 9.5716694824294557e-01, 4.8537564872284122e-01, 8.0028046888880011e-01,
        1.4188633862721534e-01, 4.2176128262627499e-01, 9.1573552518906709e-01, 7.9220732955955442e-01,
        9.5949242639290300e-01, 6.5574069915658684e-01,
    };
    std::array<double, 10> mas = {
        2.2438309411206783e-01, 5.3352367785303008e+00, 5.8684526513151845e+00, 4.2646187524686239e+00,
        4.7610216551101256e+00, 4.6692390050105752e+00, 2.4644350462159261e+00, 4.1184890487446983e+00,
        1.0755976816423436e+00, 4.4362184064364305e+00,
    };
    std::array<double, 10> eas = {
        1.0807460149582571e+00, 4.4189767041879744e+00, 5.5400461507465248e+00, 3.7847015877886920e+00,
        4.6197437955391694e+00, 4.2853529961860355e+00, 2.7845172347897078e+00, 3.6992517915192042e+00,
        1.9624400944612785e+00, 3.9582668223006068e+00,
    };

    for (int i = 0; i < eccs.size(); i++)
    {
        auto problem = eaps::KeplerProblem(eccs[i], mas[i]);
        EXPECT_NEAR(eaps::SolveKepler(problem).value(), eas[i], 1e-11);
    }
}