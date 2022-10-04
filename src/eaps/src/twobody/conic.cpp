#include "eaps/twobody/conic.hpp"

#include <cmath>
#include <stdexcept>

namespace eaps::twobody
{

ConicOrbitType ConicOrbitTypeFromEccentricity(double ecc)
{
    if (ecc < 0)
    {
        throw std::invalid_argument("negative eccentricity");
    }

    if (ecc == 0.0)
    {
        return ConicOrbitType::Circular;
    }

    if (ecc < 1.0)
    {
        return ConicOrbitType::Elliptical;
    }

    if (ecc == 1.0)
    {
        return ConicOrbitType::Parabolic;
    }
    
    return ConicOrbitType::Hyperbolic;
}

double SemiMinorAxis(double sma, double ecc)
{
    return sma * std::sqrt(1.0 - ecc * ecc);
}

double SemiLatusRectum(double sma, double ecc)
{
    return sma * (1.0 - ecc * ecc);
}

double PeriapsisRadius(double sma, double ecc)
{
    return sma * (1.0 - ecc);
}

double ApoapsisRadius(double sma, double ecc)
{
    return sma * (1.0 + ecc);
}

double ConicRadius(double sma, double ecc, double ta)
{
    double p = SemiLatusRectum(sma, ecc);
    return p / (1 + ecc * std::cos(ta));
}

} // namespace eaps::twobody
