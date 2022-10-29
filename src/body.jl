# ========================================================================================
# File   : src/body.jl
# Brief  : Definition of body type
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================

struct Body
    name::String
    naif_id::Int32
    gm::Float64

    # Ensure that gm > 0
    function Body(name, naif_id, gm)
        if gm <= 0.0
            throw(ArgumentError("non-positive gm"))
        end
        return new(name, naif_id, gm)
    end
end

# ----------------------------------------------------------------------------------------
# Accessors
# ----------------------------------------------------------------------------------------
name(b::Body) = b.name
naif_id(b::Body) = b.naif_id
gravity_parameter(b::Body) = b.gm

# ----------------------------------------------------------------------------------------
# Pre-Defined Bodies
# ----------------------------------------------------------------------------------------
const SUN = Body("sun", 10, 1.32712440041939300537e+11)
const SOLAR_SYSTEM_BARYCENTER = Body("solar_system_barycenter", 0, NaN)
const EARTH_BARYCENTER = Body("earth_barycenter", 3, 4.03503235502259805799e+05)
const MARS_BARYCENTER = Body("mars_barycenter", 4, 4.28283752140000215149e+04)
const JUPITER_BARYCENTER = Body("jupiter_barycenter", 5, 1.26712764800000205636e+08)
const SATURN_BARYCENTER = Body("saturn_barycenter", 6, 3.79405852000000029802e+07)
const URANUS_BARYCENTER = Body("uranus_barycenter", 7, 5.79454860000000800937e+06)
const NEPTUNE_BARYCENTER = Body("neptune_barycenter", 8, 6.83652710058002267033e+06)
const PLUTO_BARYCENTER = Body("pluto_barycenter", 9, 9.77000000000000682121e+02)