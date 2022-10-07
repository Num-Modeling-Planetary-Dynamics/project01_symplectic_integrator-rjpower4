#include "eaps/kepler_problem.hpp"

#include <cmath>
#include <stdexcept>

namespace eaps
{

using KP = KeplerProblem;

KP::KeplerProblem(double ecc, double ma) : ecc_{ecc}, ma_{ma}
{
    if (ecc < 0.0)
    {
        throw std::invalid_argument("negative eccentricity");
    }

    if (ecc >= 1.0)
    {
        throw std::invalid_argument("eccentricity >= 1 for kepler problem");
    }
}

double KP::Eccentricity() const
{
    return ecc_;
}

double KP::MeanAnomaly() const
{
    return ma_;
}

// ---------------------------------------------------------------------------------------
// Kepler Error and Partials Evaluation
// ---------------------------------------------------------------------------------------

static double KeplerError(const KP &problem, const double ea)
{
    return ea - problem.Eccentricity() * std::sin(ea) - problem.MeanAnomaly();
}

static double KeplerErrorD1(const KP &problem, const double ea)
{
    return 1.0 - problem.Eccentricity() * std::cos(ea);
}

static double KeplerErrorD2(const KP &problem, const double ea)
{
    return problem.Eccentricity() * std::sin(ea);
}

static double KeplerErrorD3(const KP &problem, const double ea)
{
    return problem.Eccentricity() * std::cos(ea);
}

// ---------------------------------------------------------------------------------------
// Update Methods
// ---------------------------------------------------------------------------------------

double DanbyUpdate(const KeplerProblem &problem, double ea)
{
    double err = KeplerError(problem, ea);
    double de1 = KeplerErrorD1(problem, ea);
    double de2 = KeplerErrorD2(problem, ea);
    double de3 = KeplerErrorD3(problem, ea);

    double d1 = -err / de1;
    double d2 = -err / (de1 + 0.5 * d1 * de2);
    double d3 = -err / (de1 + 0.5 * d2 * de2 + (1.0 / 6.0) * d2 * d2 * de3);

    return ea + d3;
}

double SimpleUpdate(const KeplerProblem &problem, double ea)
{
    return problem.MeanAnomaly() + problem.Eccentricity() * std::sin(ea);
}

double NewtonRaphsonUpdate(const KeplerProblem &problem, double ea)
{
    double err = KeplerError(problem, ea);
    double der = KeplerErrorD1(problem, ea);

    return ea - err / der;
}

// ---------------------------------------------------------------------------------------
// Solve
// ---------------------------------------------------------------------------------------
std::optional<double> SolveKepler(const KeplerProblem &problem, double tolerance, int max_iter)
{
    double guess = problem.MeanAnomaly();
    double error = KeplerError(problem, guess);

    int iter_count = 0;
    while (std::abs(error) > tolerance && iter_count < max_iter)
    {
        guess = DanbyUpdate(problem, guess);
        error = KeplerError(problem, guess);
        iter_count += 1;
    }

    if (iter_count == max_iter)
    {
        return std::nullopt;
    }

    return guess;
}

} // namespace eaps
