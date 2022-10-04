#ifndef EAPS_TWOBODY_KEPLER_HPP_
#define EAPS_TWOBODY_KEPLER_HPP_

#include "eaps/config.hpp"

#include <optional>
#include <stdexcept>
#include <array>

namespace eaps::twobody
{

double KeplerError(double ecc, double mean_anom, double ecc_anom);
double KeplerErrorPrime(double ecc, double mean_anom, double ecc_anom);
double KeplerErrorDoublePrime(double ecc, double mean_anom, double ecc_anom);
double KeplerErrorTriplePrime(double ecc, double mean_anom, double ecc_anom);

template <class T> class KeplerSolver
{
  public:
    KeplerSolver(double tolerance) : KeplerSolver(tolerance, 100) {}

    KeplerSolver(double tolerance, int max_iter) : tolerance_{tolerance}, max_iter_{max_iter}, iteration_method_{}
    {

        if (tolerance <= 0.0)
            throw std::invalid_argument("non-positive tolerance");
        if (max_iter <= 0)
            throw std::invalid_argument("non-positive max iter");
    }

    std::optional<double> Solve(double ecc, double mean_anom);

  private:
    double tolerance_ = 1e-12;
    int max_iter_ = 100;
    T iteration_method_;
};


template <class T> std::optional<double> KeplerSolver<T>::Solve(double ecc, double mean_anom)
{
    double guess = mean_anom;

}

} // namespace eaps::twobody

#endif // EAPS_TWOBODY_KEPLER_HPP_
