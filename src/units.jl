# ========================================================================================
# File   : src/units.jl
# Brief  : Definition of unit set and associated function
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================

# ----------------------------------------------------------------------------------------
# Unit Set Type and Constructors
# ----------------------------------------------------------------------------------------
"""
    UnitSet

Representation of a set of units defined by their relation to the Kg-Km-Sec system.
"""
struct UnitSet
    mass::Float64
    length::Float64
    time::Float64

    # Define inner constructor to enforce positive unit invariant
    function UnitSet(m, l, t)
        m > zero(m) || throw(ArgumentError("non-positive mass"))
        l > zero(l) || throw(ArgumentError("non-positive length"))
        t > zero(t) || throw(ArgumentError("non-positive time"))
        return new(m, l, t)
    end
end

"""
    UnitSet(v=1.0; mass=v, length=v, time=v)

Construct a UnitSet with a constant scaling and/or setting units individually.
"""
function UnitSet(v=1.0; mass=v, length=v, time=v)
    return UnitSet(mass, length, time)
end

"""
    unit_set(x)

Return the unit set associated with `x`.
"""
unit_set(s::UnitSet) = s

# ----------------------------------------------------------------------------------------
# Accessors
#
#   Here we also defining an interface such that any type, T,  that provides a method, 
#   `unit_set`, such that `unit_set(a::T)::UnitSet` can be passed to the `mass_unit`, 
#   `length_unit`, and `time_unit` methods thus inheriting all derived unit computations.
# ----------------------------------------------------------------------------------------
"""
    mass_unit(::UnitSet)
    mass_unit(x)

Retrieve the mass unit associated with the unit set.
"""
mass_unit(ds::UnitSet) = ds.mass
mass_unit(x) = mass_unit(unit_set(x))

"""
    length_unit(::UnitSet)
    length_unit(x)

Retrieve the length unit associated with the unit set.
"""
length_unit(ds::UnitSet) = ds.length
length_unit(x) = length_unit(unit_set(x))

"""
    time_unit(::UnitSet)
    time_unit(x)

Retrieve the time unit associated with the unit set.
"""
time_unit(ds::UnitSet) = ds.time
time_unit(x) = time_unit(unit_set(x))

# ----------------------------------------------------------------------------------------
# Derived units
# ----------------------------------------------------------------------------------------
"""
    velocity_unit(ds)

Compute the velocity scaling from the length and time scalings.
"""
velocity_unit(ds) = length_unit(ds) / time_unit(ds)

"""
    acceleration_unit(ds)

Compute the acceleration scaling.
"""
acceleration_unit(ds) = velocity_unit(ds) / time_unit(ds)

"""
    force_unit(ds)

Compute the force scaling.
"""
force_unit(ds) = mass_unit(ds) * acceleration_unit(ds)

"""
    gm_unit(ds)

Compute the gm scaling.
"""
gm_unit(ds) = length_unit(ds)^3 / time_unit(ds)^2

"""
    energy_unit(ds)

Compute the energy scaling.
"""
energy_unit(ds) = length_unit(ds) * force_unit(ds)

# ----------------------------------------------------------------------------------------
# Correctly scaled physical constants
# ----------------------------------------------------------------------------------------
gravitational_constant(ds) = GRAVITATIONAL_CONSTANT / gm_unit(ds) * mass_unit(ds)
astronomical_unit(us) = ASTRONOMICAL_UNIT / length_unit(us)
standard_gravity(us) = STANDARD_GRAVITY / acceleration_unit(us)
seconds_per_day(us) = SECONDS_PER_DAY / time_unit(us)
seconds_per_year(us) = SECONDS_PER_YEAR / time_unit(us)
speed_of_light(us) = SPEED_OF_LIGHT / velocity_unit(us)