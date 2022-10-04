#ifndef EAPS_CELESTIAL_BODY_HPP_
#define EAPS_CELESTIAL_BODY_HPP_

#include "eaps/config.hpp"

#include "eaps/body/celestial_body_type.hpp"
#include "eaps/twobody/keplerian_elements.hpp"

#include <optional>
#include <string>
#include <string_view>

namespace eaps
{

EAPS_API struct CelestialBody
{
    std::string name;
    NaifId naif_id;
    double gm;
    NaifId parent_naif_id;
};

} // namespace eaps

#endif // EAPS_CELESTIAL_BODY_HPP_