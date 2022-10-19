# ========================================================================================
# File: dimensions.jl
# Brief: definition of the dimension set used to store a scaling of standard SI dimensions
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

"""
    DimensionSet

Representation of a scaling of standard kilometer-kilogram-second units.
"""
struct DimensionSet
    mass::Float64
    length::Float64
    time::Float64
end

"""
    DimensionsSet(v=1.0; mass=v, length=v, time=v)

Construct a new dimension set ensuring validity.

This is the preferred constructor for dimension sets.
"""
function DimensionSet(v=1.0; mass=v, length=v, time=v)
    mass > 0.0 || error("non-positive masss")
    length > 0.0 || error("non-positive length")
    time > 0.0 || error("non-positive time")
    return DimensionSet(mass, length, time)
end