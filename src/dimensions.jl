# ========================================================================================
# File: dimensions.jl
# Brief: definition of the dimension set used to store a scaling of standard SI dimensions
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

# ----------------------------------------------------------------------------------------
# Custom Exception for invalid dimensions
# ----------------------------------------------------------------------------------------
struct InvalidDimensionError <: Exception
    dimension::String
    value::Float64
end

function Base.showerror(io::Base.IO, e::InvalidDimensionError)
    print(io, "invalid ", e.dimension, " dimension: ", e.value)
end

# ----------------------------------------------------------------------------------------
# Dimension Set type and External Constructors
# ----------------------------------------------------------------------------------------
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
function DimensionSet(v = 1.0; mass = v, length = v, time = v)
    mass > 0.0 || throw(InvalidDimensionError("mass", mass))
    length > 0.0 || throw(InvalidDimensionError("length", length))
    time > 0.0 || throw(InvalidDimensionError("time", time))
    return DimensionSet(mass, length, time)
end

# ----------------------------------------------------------------------------------------
# Accessors
# ----------------------------------------------------------------------------------------
mass_unit(ds::DimensionSet) = ds.mass
length_unit(ds::DimensionSet) = ds.length
time_unit(ds::DimensionSet) = ds.time

# ----------------------------------------------------------------------------------------
# Derived unit values
# ----------------------------------------------------------------------------------------
velocity_unit(ds::DimensionSet) = length_unit(ds) / time_unit(ds)
acceleration_unit(ds::DimensionSet) = velocity_unit(ds) / time_unit(ds)
force_unit(ds::DimensionSet) = acceleration_unit(ds) * mass_unit(ds)
