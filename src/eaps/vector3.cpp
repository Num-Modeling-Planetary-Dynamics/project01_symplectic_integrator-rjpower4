#include "eaps/vector3.hpp"

namespace eaps
{

double Vector3::X() const
{
    return data_[0];
}

double Vector3::Y() const
{
    return data_[1];
}

double Vector3::Z() const
{
    return data_[2];
}

} // namespace eaps