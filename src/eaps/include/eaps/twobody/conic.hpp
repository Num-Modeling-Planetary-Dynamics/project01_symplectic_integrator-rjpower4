#ifndef EAPS_TWOBODY_CONIC_HPP_
#define EAPS_TWOBODY_CONIC_HPP_

namespace eaps::twobody
{

enum class ConicOrbitType
{
    Circular,
    Elliptical,
    Parabolic,
    Hyperbolic,
};

ConicOrbitType ConicOrbitTypeFromEccentricity(double ecc);

double SemiMinorAxis(double sma, double ecc);
double SemiLatusRectum(double sma, double ecc);
double PeriapsisRadius(double sma, double ecc);
double ApoapsisRadius(double sma, double ecc);
double ConicRadius(double sma, double ecc, double ta);

} // namespace eaps::twobody

#endif // EAPS_TWOBODY_CONIC_HPP_