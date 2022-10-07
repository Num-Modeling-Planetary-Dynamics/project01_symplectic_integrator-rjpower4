#ifndef EAPS_CELESTIAL_BODY_HPP_
#define EAPS_CELESTIAL_BODY_HPP_

#include "eaps/config.hpp"

#include <string>
#include <string_view>

namespace eaps
{

EAPS_API class CelestialBody
{
  public:
    CelestialBody(std::string_view name, NaifId naif_id, double gm);

    /// @brief Return the name of the body
    const std::string &Name() const;

    /// @brief Return the NAIF (Spice) ID of the body
    NaifId NaifId() const;

    /// @brief Return the gravitational parameter of the body
    double Gm() const;

  private:
    std::string name_;
    eaps::NaifId naif_id_;
    double gm_;
};

} // namespace eaps

#endif // EAPS_CELESTIAL_BODY_HPP_