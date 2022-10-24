# ========================================================================================
# File: body.jl
# Brief: definition of body data type
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

"""
    Body
"""
struct Body
    name::String
    naif_id::Int32
    gm::Float64

    function Body(name, naif_id, gm)
        isempty(name) && throw(ArgumentError("empty body name"))
        gm <= 0.0 && throw(ArgumentError("non-positive gravitational parameter"))
        return new(name, naif_id, gm)
    end
end

Body(; name, naif_id, gm) = Body(name, naif_id, gm)

# ----------------------------------------------------------------------------------------
# Accessors
# ----------------------------------------------------------------------------------------
name(b::Body) = b.name
naif_id(b::Body) = b.naif_id
gravity_parameter(b::Body) = b.gm

# ----------------------------------------------------------------------------------------
# Derived
# ----------------------------------------------------------------------------------------
mass(b::Body) = gravity_parameter(b) / GRAVITATIONAL_CONSTANT

# ----------------------------------------------------------------------------------------
# Pre-Defined Bodies
# ----------------------------------------------------------------------------------------
const SUN = Body("sun", 10, 1.32712440041939300537e+11)

# Barycenters
const EARTH_BARYCENTER = Body("earth_barycenter", 3, 4.03503235502259805799e+05)
const MARS_BARYCENTER = Body("mars_barycenter", 4, 4.28283752140000215149e+04)
const JUPITER_BARYCENTER = Body("jupiter_barycenter", 5, 1.26712764800000205636e+08)
const SATURN_BARYCENTER = Body("saturn_barycenter", 6, 3.79405852000000029802e+07)
const URANUS_BARYCENTER = Body("uranus_barycenter", 7, 5.79454860000000800937e+06)
const NEPTUNE_BARYCENTER = Body("neptune_barycenter", 8, 6.83652710058002267033e+06)
const PLUTO_BARYCENTER = Body("pluto_barycenter", 9, 9.77000000000000682121e+02)

# Planets
const MERCURY = Body("mercury", 199, 2.20317800000000206637e+04)
const VENUS = Body("venus", 299, 3.24858592000000004191e+05)
const EARTH = Body("earth", 399, 3.98600435436095925979e+05)
const MARS = Body("mars", 499, 4.28283736206990870414e+04)
const JUPITER = Body("jupiter", 599, 1.26686534921800792217e+08)
const SATURN = Body("saturn", 699, 3.79312074986522421241e+07)
const URANUS = Body("uranus", 799, 5.79395132227900903672e+06)
const NEPTUNE = Body("neptune", 899, 6.83509950243967212737e+06)

# Dwarf Planets :(
const PLUTO = Body("pluto", 999, 8.69613817760874894702e+02)

# Moons
const MOON = Body("moon", 301, 4.90280006616379614570e+03)
const PHOBOS = Body("phobos", 401, 7.08754606689445234083e-04)
const DEIMOS = Body("deimos", 402, 9.61556964812031314061e-05)
const IO = Body("io", 501, 5.95991603341040354280e+03)
const EUROPA = Body("europa", 502, 3.20273877492289193469e+03)
const GANYMEDE = Body("ganymede", 503, 9.88783445333414420020e+03)
const CALLISTO = Body("callisto", 504, 7.17928936139727011323e+03)
const MIMAS = Body("mimas", 601, 2.50352288466179473403e+00)
const ENCELADUS = Body("enceladus", 602, 7.21129208547998867829e+00)
const TETHYS = Body("tethys", 603, 4.12111720770130247615e+01)
const DIONE = Body("dione", 604, 7.31163532292319331418e+01)
const RHEA = Body("rhea", 605, 1.53942204554534214367e+02)
const TITAN = Body("titan", 606, 8.97813884530737595924e+03)
const HYPERION = Body("hyperion", 607, 3.71879171419166820733e-01)
const ARIEL = Body("ariel", 701, 8.34634443177047700146e+01)
const UMBRIEL = Body("umbriel", 702, 8.50933809448938802689e+01)
const TITANIA = Body("titania", 703, 2.26943700374124802011e+02)
const OBERON = Body("oberon", 704, 2.05323430253562293046e+02)
const TRITON = Body("triton", 801, 1.42759814072503399984e+03)
const CHARON = Body("charon", 901, 1.05879988860188106514e+02)
