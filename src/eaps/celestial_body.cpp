#include "eaps/celestial_body.hpp"

#include <stdexcept>

namespace eaps
{

CelestialBody::CelestialBody(std::string_view name, eaps::NaifId naif_id, double gm)
    : name_{name}, naif_id_{naif_id}, gm_{gm}
{
    if (gm <= 0)
    {
        throw std::invalid_argument("non-positive gm");
    }
}

const std::string &CelestialBody::Name() const
{
    return name_;
}

NaifId CelestialBody::NaifId() const
{
    return naif_id_;
}

double CelestialBody::Gm() const
{
    return gm_;
}

} // namespace eaps
