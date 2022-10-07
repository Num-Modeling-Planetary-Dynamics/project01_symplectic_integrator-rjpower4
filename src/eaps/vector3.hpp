#ifndef EAPS_VECTOR_3_HPP_
#define EAPS_VECTOR_3_HPP_

#include "eaps/config.hpp"

#include <array>

namespace eaps
{

EAPS_API class Vector3
{
  public:
    Vector3(double x, double y, double z) : data_{x, y, z} {}

    double X() const;
    double Y() const;
    double Z() const;

  private:
    std::array<double, 3> data_;
};

} // namespace eaps

#endif // EAPS_VECTOR_3_HPP_