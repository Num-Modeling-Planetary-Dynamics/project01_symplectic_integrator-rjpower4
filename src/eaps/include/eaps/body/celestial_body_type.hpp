#ifndef EAPS_CELESTIAL_BODY_TYPE_HPP_
#define EAPS_CELESTIAL_BODY_TYPE_HPP_

#include "eaps/config.hpp"

namespace eaps
{

enum class CelestialBodyType
{
    Barycenter,
    Star,
    Planet,
    DwarfPlanet,
    Moon,
    Asteroid,
    Comet,
    Unknown,
};

} // namespace eaps

#endif // EAPS_CELESTIAL_BODY_TYPE_HPP_