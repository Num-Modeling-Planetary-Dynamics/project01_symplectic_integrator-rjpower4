#include "eaps/twobody/keplerian_elements.hpp"

#include <stdexcept>
#include <cmath>

#include "eaps/twobody/conic.hpp"

namespace eaps::twobody
{

KeplerianElements::KeplerianElements(double sma, double ecc, double inc, double aop, double raan, double ta)
    : sma_{sma}, ecc_{ecc}, inc_{inc}, aop_{aop}, raan_{raan}, ta_{ta}
{
    if (ecc < 0.0) {
        throw std::invalid_argument("negative eccentricity");
    }

    if (sma == 0) {
        throw std::invalid_argument("zero semi-major axis");
    }

    if (sma < 0.0 && ecc < 1.0) {
        throw std::invalid_argument("negative semi-major axis, elliptical eccentricity");
    }

    if (sma > 0.0 && ecc > 1.0) {
        throw std::invalid_argument("positive semi-major axis, hyperbolic eccentricity");
    }

    if (std::isinf(sma) && ecc != 1.0) {
        throw std::invalid_argument("infinite semi-major axis with non-unit eccentricity");
    }

    if (!std::isinf(sma) && ecc == 1.0) {
        throw std::invalid_argument("finite semi-major axis with unit eccentricity");
    }
}

double KeplerianElements::SemiMajorAxis() const noexcept
{
    return sma_;
}

double KeplerianElements::Eccentricity() const noexcept
{
    return ecc_;
}

double KeplerianElements::Inclination() const noexcept
{
    return inc_;
}

double KeplerianElements::ArgumentOfPeriapsis() const noexcept
{
    return aop_;
}

double KeplerianElements::RightAscension() const noexcept
{
    return raan_;
}

double KeplerianElements::TrueAnomaly() const noexcept
{
    return ta_;
}

double KeplerianElements::SemiMinorAxis() const noexcept
{
    return eaps::twobody::SemiMinorAxis(SemiMajorAxis(), Eccentricity());
}

double KeplerianElements::SemiLatusRectum() const noexcept
{
    return eaps::twobody::SemiLatusRectum(SemiMajorAxis(), Eccentricity());
}

double KeplerianElements::PeriapsisRadius() const noexcept
{
    return eaps::twobody::PeriapsisRadius(SemiMajorAxis(), Eccentricity());
}

double KeplerianElements::ApoapsisRadius() const noexcept
{
    return eaps::twobody::ApoapsisRadius(SemiMajorAxis(), Eccentricity());
}

double KeplerianElements::RadiusAt(double ta) const noexcept
{
    return eaps::twobody::ConicRadius(SemiMajorAxis(), Eccentricity(), ta);
}

double KeplerianElements::Radius() const noexcept
{
    return RadiusAt(TrueAnomaly());
}

} // namespace eaps::twobody