#ifndef EAPS_KEPLERIAN_ELEMENTS_HPP_
#define EAPS_KEPLERIAN_ELEMENTS_HPP_

#include "eaps/config.hpp"

namespace eaps::twobody
{

EAPS_API class KeplerianElements
{
  public:
    explicit KeplerianElements(double sma, double ecc = 0.0, double inc = 0.0, double aop = 0.0, double raan = 0.0,
                               double ta = 0.0);

    [[nodiscard]] double SemiMajorAxis() const noexcept;
    [[nodiscard]] double Eccentricity() const noexcept;
    [[nodiscard]] double Inclination() const noexcept;
    [[nodiscard]] double ArgumentOfPeriapsis() const noexcept;
    [[nodiscard]] double RightAscension() const noexcept;
    [[nodiscard]] double TrueAnomaly() const noexcept;

    [[nodiscard]] double SemiMinorAxis() const noexcept;
    [[nodiscard]] double SemiLatusRectum() const noexcept;
    [[nodiscard]] double PeriapsisRadius() const noexcept;
    [[nodiscard]] double ApoapsisRadius() const noexcept;

    [[nodiscard]] double RadiusAt(double ta) const noexcept;
    [[nodiscard]] double Radius() const noexcept;

  private:
    double sma_;  ///< Semi-major axis
    double ecc_;  ///< Eccentricity
    double inc_;  ///< Inclination
    double aop_;  ///< Argument of Periapsis
    double raan_; ///< Right Ascension
    double ta_;   ///< True Anomaly
};

} // namespace eaps::twobody

#endif // EAPS_KEPLERIAN_ELEMENTS_HPP_