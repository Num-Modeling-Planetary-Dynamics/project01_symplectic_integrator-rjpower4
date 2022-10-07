#ifndef EAPS_KEPLER_PROBLEM_HPP_
#define EAPS_KEPLER_PROBLEM_HPP_

#include "eaps/config.hpp"

#include <optional>

namespace eaps
{

EAPS_API class KeplerProblem
{
  public:
    KeplerProblem(double ecc, double ma);

    double Eccentricity() const;
    double MeanAnomaly() const;

  private:
    double ecc_; ///< Eccentricity
    double ma_;  ///< Mean Anomaly
};

std::optional<double> SolveKepler(const KeplerProblem &problem, double tolerance = 1e-12, int max_iter = 10);

} // namespace eaps

#endif // EAPS_KEPLER_PROBLEM_HPP_