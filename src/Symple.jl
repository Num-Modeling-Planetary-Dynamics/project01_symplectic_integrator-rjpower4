module Symple

# ----------------------------------------------------------------------------------------
# External Packages
# ----------------------------------------------------------------------------------------
using StaticArrays: FieldVector, SVector, SMatrix
using SPICE: furnsh, spkez, spkezr, str2et
using LinearAlgebra: norm, dot, cross

# ----------------------------------------------------------------------------------------
# Exports
# ----------------------------------------------------------------------------------------
export Simulation, SimulationState
export LagrangianState, HamiltonianState
export SUN, NEPTUNE_BARYCENTER, PLUTO_BARYCENTER
export naif_id, spice_state, load_kernels
export kepler_danby

# ----------------------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------------------
"""
GRAVITATIONAL_CONSTANT

Universal Constant of Gravitation

This is also known as the Newtonian Constant of Gravitation and the Cavendish
Gravitational Constant. In the classical (Newtonian) formulation of gravity,
this value serves as the constant of proportionality between the force
exterted between two bodies and the sum of their products divided by the
square of the distance between them.

Note that the value of this constant is only really known to about four
significant figures. Therefore, if its use can be avoided it should be.

* Dimensions : [M^(-1) L^3 T^(-2)]
* Units      : [km^3 / (kg . s^2)]
* Sources    : [Nist](https://physics.nist.gov/cuu/Constants/)
               [ArXiv](https://arxiv.org/abs/1507.07956)
"""
const GRAVITATIONAL_CONSTANT = 6.6743e-20

"""
    ASTRONOMICAL_UNIT

Astronomical unit

Unit of length that is about equal to the distance from the Earth to the Sun.
According to WikiPedia, the distance from the Earth to the Sun varies by about
3% over the course of a year so this value could be (in theory) about 3% off
the distance for any given epoch. However, the use case of this unit is that
it is a convenient scaling factor for distances of solar system magnitudes as
it is officially defined.

* Dimensions : [L]
* Units      : [km]
* Sources    : [IAU](https://www.iau.org/static/resolutions/IAU2012_English.pdf)
"""
const ASTRONOMICAL_UNIT = 149597870.700

"""
    SECONDS_PER_DAY

Seconds in a Julian Day

Number of seconds in a Julian day. It's what it says on the tin...

* Dimensions : [T]
* Units      : [s]
* Sources    : SPICE
"""
const SECONDS_PER_DAY = 86400.0

"""
    SECONDS_PER_YEAR

Seconds in a Julian Year

Number of seconds in a Julian year. Did you also get confused about the
SECONDS_PER_DAY one...?

* Dimensions : [T]
* Units      : [s]
* Sources    : SPICE
"""
const SECONDS_PER_YEAR = 31557600.0

# ----------------------------------------------------------------------------------------
# Dimension Set
# ----------------------------------------------------------------------------------------
"""
    InvalidDimensionError

Custom Exception for invalid dimensions
"""
struct InvalidDimensionError <: Exception
    dimension::String
    value::Float64
end

function Base.showerror(io::Base.IO, e::InvalidDimensionError)
    return print(io, "invalid ", e.dimension, " dimension: ", e.value)
end

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

# Accessors
mass_unit(ds::DimensionSet) = ds.mass
length_unit(ds::DimensionSet) = ds.length
time_unit(ds::DimensionSet) = ds.time

# Interface
mass_unit(ds_provider) = mass_unit(dimension_set(ds_provider))
length_unit(ds_provider) = length_unit(dimension_set(ds_provider))
time_unit(ds_provider) = time_unit(dimension_set(ds_provider))

# Derived unit values
velocity_unit(ds) = length_unit(ds) / time_unit(ds)
acceleration_unit(ds) = velocity_unit(ds) / time_unit(ds)
force_unit(ds) = acceleration_unit(ds) * mass_unit(ds)
gm_unit(ds) = time_unit(ds)^2 / length_unit(ds)^3

# ----------------------------------------------------------------------------------------
# Body
# ----------------------------------------------------------------------------------------
"""
    Body

Abstraction of a celestial body/barycenter
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
naif_id(b::Body) = b.naif_id
gravity_parameter(b::Body) = b.gm

# ----------------------------------------------------------------------------------------
# Pre-defined bodies
# ----------------------------------------------------------------------------------------
const SUN = Body("sun", 10, 1.32712440041939300537e+11)
const SOLAR_SYSTEM_BARYCENTER = Body("solar_system_barycenter", 0, NaN)

# Barycenters
const EARTH_BARYCENTER = Body("earth_barycenter", 3, 4.03503235502259805799e+05)
const MARS_BARYCENTER = Body("mars_barycenter", 4, 4.28283752140000215149e+04)
const JUPITER_BARYCENTER = Body("jupiter_barycenter", 5, 1.26712764800000205636e+08)
const SATURN_BARYCENTER = Body("saturn_barycenter", 6, 3.79405852000000029802e+07)
const URANUS_BARYCENTER = Body("uranus_barycenter", 7, 5.79454860000000800937e+06)
const NEPTUNE_BARYCENTER = Body("neptune_barycenter", 8, 6.83652710058002267033e+06)
const PLUTO_BARYCENTER = Body("pluto_barycenter", 9, 9.77000000000000682121e+02)

# ----------------------------------------------------------------------------------------
# Base state
# ----------------------------------------------------------------------------------------
"""
    AbstractState{T}

Abstract state type for 6-element states.
"""
abstract type AbstractState{T} <: FieldVector{6, T} end

"""
    coordinates(s::State)

Return the coordinates (position components) of the state.
"""
coordinates(s::AbstractState) = SVector{3}(s[1], s[2], s[3])

# ----------------------------------------------------------------------------------------
# Lagrangian State
# ----------------------------------------------------------------------------------------
"""
    LagrangianState

State in the Lagrangian basis, i.e., using coordinates and velocities.
"""
struct LagrangianState <: AbstractState{Float64}
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
end

"""
    LagrangianState(r::AbstractVector, v::AbstractVector)

Construct a Lagrangian state from the position and velocity vectors
"""
function LagrangianState(r::AbstractVector, v::AbstractVector)
    return LagrangianState(r[1], r[2], r[3], v[1], v[2], v[3])
end

"""
    velocities(ls::LagrangianState)

Return the velocity measure numbers of the state.
"""
velocities(ls::LagrangianState) = SVector{3}(ls[4], ls[5], ls[6])

"""
    radial_velocity(ls::LagrangianState)

Return the signed scalar radial velocity of the state.
"""
function radial_velocity(ls::LagrangianState)
    pos = coordinates(ls)
    vel = velocities(ls)
    return dot(pos, vel) / norm(pos)
end

"""
    specific_angular_momentum(ls::LagrangianState)

Return the specific angular momentum.
"""
function specific_angular_momentum(ls::LagrangianState)
    pos = coordinates(ls)
    vel = coordinates(ls)
    return cross(pos, vel)
end

# ----------------------------------------------------------------------------------------
# Hamiltonian State
# ----------------------------------------------------------------------------------------
"""
    HamiltonianState

State in the Hamiltonian basis, i.e., using coordinates and momenta.
"""
struct HamiltonianState <: AbstractState{Float64}
    x::Float64
    y::Float64
    z::Float64
    px::Float64
    py::Float64
    pz::Float64
end

"""
    HamiltonianState(r::AbstractVector, p::AbstractVector)

Construct a Hamiltonian state from the position and momenta vectors
"""
function HamiltonianState(r::AbstractVector, p::AbstractVector)
    return HamiltonianState(r[1], r[2], r[3], p[1], p[2], p[3])
end

"""
    momenta(::HamiltonianState)

Return the momenta components of the Hamiltonian state.
"""
function momenta(hs::HamiltonianState)
    return SVector{3}(hs[4], hs[5], hs[6])
end

# ----------------------------------------------------------------------------------------
# Conversion between Lagrangian and Hamiltonian states
# ----------------------------------------------------------------------------------------
"""
    LagrangianState(m, s::HamiltonianState)

Convert the Hamiltonian state to a Lagrangian state.
"""
function LagrangianState(m, s::HamiltonianState)
    return LagrangianState(coordinates(s), momenta(s) / m)
end

"""
    HamiltonianState(m, s::HamiltonianState)

Convert the Lagrangian state to a Hamiltonian state.
"""
function HamiltonianState(m, s::LagrangianState)
    return HamiltonianState(coordinates(s), velocities(s) * m)
end

# ----------------------------------------------------------------------------------------
# Two-Body Utilities
# ----------------------------------------------------------------------------------------
function circular_speed(gm, r)
    return sqrt(gm / r)
end

function specific_energy(gm, a::AbstractFloat)
    return -gm / (2a)
end

function specific_energy(gm, r::AbstractFloat, v::AbstractFloat)
    return v^2 / 2 - gm / r
end

function specific_energy(gm, pos::AbstractVector, vel::AbstractVector)
    return specific_energy(gm, norm(pos), norm(vel))
end

function specific_energy(b::Body, pos::AbstractVector, vel::AbstractVector)
    return specific_energy(gravity_parameter(b), norm(pos), norm(vel))
end

function specific_energy(gm, x::LagrangianState)
    return specific_energy(gm, coordinates(x), velocities(x))
end

function specific_energy(b::Body, x::LagrangianState)
    return specific_energy(gravity_parameter(b), x)
end

function semi_major_axis(gm, spec_energy)
    return -gm / (2spec_energy)
end

function semi_minor_axis(sma, ecc)
    return sma * sqrt(1 - ecc^2)
end

function semi_latus_rectum(sma, ecc)
    return sma * (1 - ecc^2)
end

function radius(sma, ecc, ta)
    p = semi_latus_rectum(sma, ecc)
    return p / (1 + ecc * cos(ta))
end

function mean_motion(gm, sma)
    return sqrt(gm / sma^3)
end

function eccentric_anomaly(ecc, ta)
    ta_bound = mod(ta, 2π)
    alpha = sqrt((1 + ecc) / (1 - ecc))
    tane2 = tan(ta_bound / 2) / alpha
    return mod(2atan(tane2), 2π)
end

function mean_anomaly(ecc, eccentric_anomaly)
    return eccentric_anomaly - ecc * sin(eccentric_anomaly)
end

function period(gm, sma)
    return 2π / mean_motion(gm, sma)
end

function inclination(angular_momentum::AbstractVector)
    return acos(angular_momentum[3] / norm(angular_momentum))
end

function perifocal_to_inertial_rotation(inc, aop, raan)
    sa, ca = sincos(aop)
    si, ci = sincos(inc)
    sr, cr = sincos(raan)

    r11 = cr * ca - sr * sa * ci
    r12 = -cr * sa - sr * ca * ci
    r13 = sr * si

    r21 = sr * ca + cr * sa * ci
    r22 = -sr * sa + cr * ca * ci
    r23 = -cr * si

    r31 = sa * si
    r32 = ca * si
    r33 = ci

    return SMatrix{3, 3}(r11, r21, r31, r12, r22, r32, r13, r23, r33)
end

# ----------------------------------------------------------------------------------------
# Keplerian Elements
# ----------------------------------------------------------------------------------------
"""
    KeplerianElements

Collection of Keplerian elements
"""
struct KeplerianElements <: FieldVector{6, Float64}
    sma::Float64
    ecc::Float64
    inc::Float64
    aop::Float64
    raan::Float64
    ta::Float64
end

function KeplerianElements(sma; ecc = 0.0, inc = 0.0, aop = 0.0, raan = 0.0, ta = 0.0)
    sma != 0.0 || throw(ArgumentError("invalid semi-major axis: $sma"))
    ecc >= 0.0 || throw(ArgumentError("invalid eccentricity: $ecc"))

    if sma .* (ecc - 1.0) > 0.0
        throw(ArgumentError("invalid semi-major axis and eccentricity"))
    elseif ecc == 1.0 && !isinf(sma)
        throw(ArgumentError("invalid semi-major axis and eccentricity"))
    end

    return KeplerianElements(sma, ecc, inc, aop, raan, ta)
end

semi_major_axis(ke::KeplerianElements) = ke.sma
eccentricity(ke::KeplerianElements) = ke.ecc
inclination(ke::KeplerianElements) = ke.inc
argument_of_periapsis(ke::KeplerianElements) = ke.aop
right_ascension(ke::KeplerianElements) = ke.raan
true_anomaly(ke::KeplerianElements) = ke.ta

function specific_energy(gm, ke::KeplerianElements)
    return specific_energy(gm, semi_major_axis(ke))
end

function semi_minor_axis(ke::KeplerianElements)
    return semi_minor_axis(semi_major_axis(ke), eccentricity(ke))
end

function semi_latus_rectum(ke::KeplerianElements)
    return semi_latus_rectum(semi_major_axis(ke), eccentricity(ke))
end

function mean_motion(gm, ke::KeplerianElements)
    return mean_motion(gm, semi_major_axis(ke))
end

function period(gm, ke::KeplerianElements)
    return period(gm, semi_major_axis(ke))
end

function period(b::Body, ke::KeplerianElements)
    return period(gravity_parameter(b), ke)
end

function eccentric_anomaly(ke::KeplerianElements; ta = true_anomaly(ke))
    return eccentric_anomaly(eccentricity(ke), ta)
end

function radius(ke::KeplerianElements; ta = true_anomaly(ke))
    return radius(semi_major_axis(ke), eccentricity(ke), ta)
end

function perifocal_to_inertial_rotation(ke::KeplerianElements)
    inc = inclination(ke)
    aop = argument_of_periapsis(ke)
    raan = right_ascension(ke)
    return perifocal_to_inertial_rotation(inc, aop, raan)
end

"""
    KeplerianElements(gm, pos::AbstractVector, vel::AbstractVector)

Convert from position and velocity to Keplerian elements.
"""
function KeplerianElements(gm, posv::AbstractVector, velv::AbstractVector)
    pos = SVector{3}(posv)
    vel = SVector{3}(velv)
    unit = SVector{3}(one(eltype(pos)), zero(eltype(pos)), zero(eltype(pos)))

    normalize(v) = v ./ norm(v)
    angle_between(v1, v2) = acos(dot(v1, v2))
    wrap_if(val, condition) = condition ? 2pi - val : val

    r = norm(pos)
    v = norm(vel)

    rhat = pos ./ r

    sma = 1 / (2 / r - v^2 / gm)
    h = cross(pos, vel)
    hmag = norm(h)
    ecc = sqrt(1 - hmag^2 / (gm * sma))
    inc = inclination(h)

    ecc_vector = cross(vel, h) / gm - rhat
    # ecc = sqrt(1 - norm(h)^2 / (gm * sma))
    e_hat = ecc != 0.0 ? ecc_vector ./ ecc : unit

    nodal_axis = SVector{3}(-h[2], h[1], 0.0)
    nodal_axis_hat = nodal_axis ./ norm(nodal_axis)
    # nodal_axis_hat = h[2] == 0 ? unit : normalize(nodal_axis)

    raan = wrap_if(acos(nodal_axis_hat[1]), nodal_axis_hat[2] < 0.0)
    aop = wrap_if(angle_between(nodal_axis_hat, e_hat), ecc_vector[3] < 0.0)
    ta = acos(clamp(dot(e_hat, rhat), -1.0, 1.0))

    return KeplerianElements(sma, ecc, inc, aop, raan, ta)
end

function KeplerianElements(gm, x::LagrangianState)
    return KeplerianElements(gm, coordinates(x), velocities(x))
end

function KeplerianElements(b::Body, x::LagrangianState)
    return KeplerianElements(gravity_parameter(b), x)
end

"""
    LagrangianState(gm, ke::KeplerianElements)

Convert Keplerian elements into the 6-dimensional state vector.
"""
function LagrangianState(gm, ke::KeplerianElements)
    z = zero(eltype(ke))
    rot_mat = perifocal_to_inertial_rotation(inclination(ke),
                                             argument_of_periapsis(ke) + true_anomaly(ke),
                                             right_ascension(ke))
    ta = true_anomaly(ke)
    p = semi_latus_rectum(ke)
    h = sqrt(gm * p)

    r_peri = SVector{3}(radius(ke), z, z)
    r_inert = rot_mat * r_peri

    v_transverse = h / norm(r_peri)
    v_radial = h / p * eccentricity(ke) * sin(ta)
    v_peri = SVector{3}(v_radial, v_transverse, z)
    v_inert = rot_mat * v_peri

    return LagrangianState(r_inert, v_inert)
end

# ----------------------------------------------------------------------------------------
# Solving Kepler's problem
# ----------------------------------------------------------------------------------------
"""
    kepler_initial_guess(eccentricity, mean_anomaly)

Return an initial guess to the eccentric anomaly given the mean anomaly and eccentricity.

## References
* Danby, J. M. A., and Burkardt, T. M. “The Solution of Kepler’s Equation, I.” Celestial
  Mechanics, Vol. 31, No. 2, 1983, pp. 95-107. https://doi.org/10.1007/BF01686811.
"""
function _kepler_initial_guess(eccentricity, mean_anomaly)
    return mean_anomaly
end

"""
    _kepler_error(ecc, ma, ea)

Compute error in Kepler problem along with first 3 derivatives w.r.t. eccentric anomaly.
"""
function _kepler_error(ecc, ma, ea)
    se, ce = sincos(ea)
    ese = ecc * se
    ece = ecc * ce
    return (ea - ese - ma, 1 - ece, ese, ece)
end

"""
    kepler_danby(eccentricity, mean_anomaly; tolerance=1e-12, max_iter=10)

Solve Kepler's problem returning the eccentric anomaly associated with given mean anomaly.

Uses Danby's method for solving the problem.

## References
* Danby, J. M. A. “The Solution of Kepler's Equation, III.” Celestial Mechanics, Vol. 40,
  Nos. 3-4, 1987, pp. 303-312. https://doi.org/10.1007/BF01235847.
* Danby, J. M. A., and Burkardt, T. M. “The Solution of Kepler's Equation, I.” Celestial
  Mechanics, Vol. 31, No. 2, 1983, pp. 95-107. https://doi.org/10.1007/BF01686811.
* Burkardt, T. M., and Danby, J. M. A. “The Solution of Kepler's Equation, II.” Celestial
  Mechanics, Vol. 31, No. 3, 1983, pp. 317-328. https://doi.org/10.1007/BF01844230.
"""
function kepler_danby(eccentricity, mean_anomaly; tolerance = 1e-12, max_iter = 10)
    ea = _kepler_initial_guess(eccentricity, mean_anomaly)

    (f, f1, f2, f3) = _kepler_error(eccentricity, mean_anomaly, ea)

    iter_count = 0
    while abs(f) > tolerance && iter_count < max_iter
        d1 = -f / f1
        d2 = -f / (f1 + (1 // 2) * d1 * f2)
        d3 = -f / (f1 + (1 // 2) * d2 * f2 + (1 // 6) * d2^2 * f3)

        ea += d3
        iter_count += 1

        (f, f1, f2, f3) = _kepler_error(eccentricity, mean_anomaly, ea)
    end

    return (ea, abs(f) <= tolerance)
end

# ----------------------------------------------------------------------------------------
# Kepler Propagation
# ----------------------------------------------------------------------------------------
"""
    angular_difference(a, b)

Return the correctly signed counter-clockwise angular difference from a to b.
"""
function angular_difference(a, b)
    ay, ax = sincos(a)
    by, bx = sincos(b)
    theta = atan(ax * by - ay * bx, ax * bx + ay * by)
    return theta > 0 ? theta : 2π + theta
end

"""
    kepler_propagate(gm, dt, x0::LagrangianState)

Propagate the state in the Keplerian two-body model.

## Arguments

- `gm`: Gravitational parameter of the central body
- `dt`: Time to propagate for (may be negative)
- `x0::LagrangianState`: Initial state of the orbiting body
"""
function kepler_propagate(gm, dt, x0::LagrangianState)
    els = KeplerianElements(gm, x0)
    a = semi_major_axis(els)
    n = mean_motion(gm, els)
    dt = mod(dt, 2π / n)

    ea_initial = eccentric_anomaly(els)
    ma_initial = ea_initial - eccentricity(els) * sin(ea_initial)

    ma_final = ma_initial + n * dt
    ea_final, did_converge = kepler_danby(eccentricity(els), ma_final)
    if !did_converge
        error("unable to converge eccentric anomaly with ma=$(ma_final), gm=$gm x0=$(x0)")
    end
    err = ma_final - ea_final + eccentricity(els) * sin(ea_final)

    de = angular_difference(ea_initial, ea_final)
    sde, cde = sincos(de)

    q0 = coordinates(x0)
    dq0 = velocities(x0)

    r0 = norm(q0)

    f = a / r0 * (cde - 1) + 1
    g = dt + (sde - (ea_final - ea_initial)) / n

    new_coordinates = f * q0 + g * dq0
    r = norm(new_coordinates)

    df = -a^2 / (r * r0) * n * sde
    dg = a / r * (cde - 1) + 1

    new_velocities = df * q0 + dg * dq0

    dEn = specific_energy(gm, x0) -
          specific_energy(gm, LagrangianState(new_coordinates, new_velocities))
    println("$(err) -- $(dEn)")

    return LagrangianState(new_coordinates, new_velocities)
end

function kepler_propagate2(gm, dt, x0::LagrangianState)
    rv0 = coordinates(x0)
    vv0 = velocities(x0)

    rmag0 = norm(rv0)
    vmag2 = dot(vv0, vv0)
    h = cross(rv0, vv0)
    hmag2 = dot(h, h)
    a = 1 / (2 / rmag0 - vmag2 / gm)
    ecc = sqrt(1 - hmag2 / (gm * a))

    n = sqrt(gm / a^3)

    E0 = if ecc > 1e-12
        acos(-(rmag0 - a) / (a * ecc))
    else
        0.0
    end

    if dot(rv0, vv0) < -1e-12
        E0 = 2π - E0
    end

    M0 = E0 - ecc * sin(E0)
    M = M0 + n * dt

    E, _ = kepler_danby(ecc, M)
    dE = E - E0

    f = a / rmag0 * (cos(dE) - 1.0) + 1.0
    g = dt + (sin(dE) - dE) / n

    rv = f * rv0 + g * vv0
    rmag = norm(rv)

    fdot = -a^2 / (rmag * rmag0) * n * sin(dE)
    gdot = a / rmag * (cos(dE) - 1.0) + 1.0

    vv = fdot * rv0 + gdot * vv0

    return LagrangianState(rv, vv)
end

# ----------------------------------------------------------------------------------------
# Spice Utilities
# ----------------------------------------------------------------------------------------
"""
    SPICE_KERNEL_EXTENSIONS

Tuple of extensions assoicated with spice kernels
"""
const SPICE_KERNEL_EXTENSIONS = (".bsp", ".tpc", ".tls")

"""
    is_kernel_path(path; require_exists=true)

Determine if the path given points to a valid spice kernel based on path only.

If `require_exists` is `false`, then only extension of file is checked. If it is true,
then also return `false` if the file does not exist.
"""
function is_kernel_path(path::String; require_exists = true)
    _, ext = splitext(path)
    return (isfile(path) || !require_exists) && (ext in SPICE_KERNEL_EXTENSIONS)
end

"""
    list_kernels(path::String)

Return a vector of paths to all kernels in directory.
"""
function list_kernels(path::String)
    isdir(path) || ArgumentError("path not a directory: " + path)
    (r, _, files) = first(walkdir(normpath(path)))
    root = abspath(r)
    kernels = String[]
    for p in files
        pth = joinpath(root, p)
        if is_kernel_path(pth)
            push!(kernels, pth)
        end
    end
    return kernels
end

"""
    load_kernel(pth::String)

Load the kernel at the given path throwing error if not found.
"""
function load_kernel(pth::String)
    isfile(path) || ArgumentError("kernel file not found: " + pth)
    furnsh(pth)
    return nothing
end

"""
    load_kernels(paths::AbstractVector)

Load all kernels with path in `paths`. Return error and don't load any if any not found.
"""
function load_kernels(pths::AbstractVector)
    if any(x -> !isfile(x), pths)
        error("kernel(s) not found")
    end
    foreach(furnsh, pths)
    return nothing
end

"""
    ephemeris_time(s::String)

Compute the ephemeris time from the time string.
"""
ephemeris_time(s::String) = str2et(s)

"""
    spice_state(id::Integer, et; frame="J2000", observer=0)

Retrieve the state of the body with id `id` at the given ephemeris time.
"""
function spice_state(id::Integer, et; frame = "J2000", observer = 0)
    out, _ = spkez(id, et, frame, "NONE", observer)
    return out
end

"""
    spice_state(name::String, et; frame="J2000", observer="SOLAR_SYSTEM_BARYCENTER")

Retrieve the state of the body with name `name` at the given ephemeris time.
"""
function spice_state(name::String, et; frame = "J2000",
                     observer = "SOLAR_SYSTEM_BARYCENTER")
    out, _ = spkezr(name, et, frame, "NONE", observer)
    return out
end

# ----------------------------------------------------------------------------------------
# Simulation
# ----------------------------------------------------------------------------------------
struct Simulation
    primary::Body
    secondaries::Vector{Body}
    dimensions::DimensionSet
end

function Simulation(primary, secondaries; dimensions = DimensionSet())
    return Simulation(primary, secondaries, dimensions)
end

primary(s::Simulation) = s.primary
secondaries(s::Simulation) = s.secondaries
dimension_set(s::Simulation) = s.dimensions

function gravity_parameters(s::Simulation)
    gmu = gm_unit(s)
    gm_p = (gravity_parameter(primary(s))) ./ gmu
    gm_s = (gravity_parameter.(secondaries(s))) ./ gmu
    return (gm_p, gm_s)
end

function total_gravity_parameter(s::Simulation)
    gmu = gm_unit(s)
    gm_total = gravity_parameter(primary(s))
    for b in secondaries(s)
        gm_total += gravity_parameter(b)
    end
    return gm_total / gmu
end

function gravity_parameter(s::Simulation, i::Integer)
    gm_raw = if i == 0
        gravity_parameter(primary(s))
    else
        gravity_parameter(secondaries(s)[i])
    end
    return gm_raw / gm_unit(s)
end

# ----------------------------------------------------------------------------------------
# Simulation State 
# ----------------------------------------------------------------------------------------
mutable struct SimulationState
    primary_state::LagrangianState
    secondary_states::Vector{LagrangianState}
end

primary(s::SimulationState) = s.primary_state
secondaries(s::SimulationState) = s.secondary_states

# ----------------------------------------------------------------------------------------
# Heliocentric <-> Barycentric
# ----------------------------------------------------------------------------------------
function heliocentric_to_barycentric(s::Simulation, x::SimulationState)
    gm_0, gm_i = gravity_parameters(s)
    gm_total = gm_0 + sum(gm_i)

    x_0 = primary(x)
    x_i = secondaries(x)

    q_bary_0 = gm_0 * coordinates(x_0)
    v_bary_0 = gm_0 * velocities(x_0)

    for (gm, x) in zip(gm_i, x_i)
        q_bary_0 += coordinates(x)
        v_bary_0 += gm * velocities(x)
    end
    x_bary_0 = LagrangianState(q_bary_0, v_bary_0) / gm_total

    x_bary_i = Vector{LagrangianState}(undef, length(x_i))
    for (k, x) in enumerate(x_i)
        x_bary_i[k] = LagrangianState(coordinates(x) - coordinates(x_bary_0),
                                      velocities(x) - velocities(x_bary_0))
    end

    return SimulationState(x_bary_0, x_bary_i)
end

function barycentric_to_heliocentric(s::Simulation, x::LagrangianState) end

# ----------------------------------------------------------------------------------------
# Dynamics
# ----------------------------------------------------------------------------------------
function linear_drift!(ss::SimulationState, s::Simulation, dt)
    z = SVector{3}(0.0, 0.0, 0.0)
    gm_0 = gravity_parameter(s, 0)
    gm_total = total_gravity_parameter(s)

    # dq_0 = gm_0 * velocities(primary(ss)) / gm_total
    # ss.primary_state += dt * LagrangianState(dq_0, z)

    for (i, x) in enumerate(secondaries(ss))
        q_dot = gravity_parameter(s, i) * velocities(x) / gm_0
        ss.secondary_states[i] += dt * LagrangianState(q_dot, z)
    end

    return nothing
end

function interaction_kick!(ss::SimulationState, s::Simulation, dt)
    secs = secondaries(s)
    n_secondaries = length(secs)
    for i in 1:(n_secondaries - 1)
        q_i = coordinates(secondaries(ss)[i])
        dx_i = LagrangianState(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        gm_i = gravity_parameter(s, i)

        for j in 2:n_secondaries
            gm_j = gravity_parameter(s, j)
            q_j = coordinates(secondaries(ss)[j])
            q_ij = q_j - q_i
            r = norm(q_ij)
            dx_i += LagrangianState(SVector{3}(0.0, 0.0, 0.0), gm_j / r^3 * q_ij)
        end
        ss.secondary_states[i] += dt * dx_i
    end
end

function kepler_drift!(ss::SimulationState, s::Simulation, dt)
    states = secondaries(ss)
    gm_0 = gravity_parameter(s, 0)

    for (i, x) in enumerate(states)
        x_new = kepler_propagate2(gm_0, dt, x)
        ss.secondary_states[i] = x_new
    end
end

# ----------------------------------------------------------------------------------------
# Stepping
# ----------------------------------------------------------------------------------------
function do_step!(ss::SimulationState, s::Simulation, dt)
    linear_drift!(ss, s, dt / 2)
    interaction_kick!(ss, s, dt / 2)
    kepler_drift!(ss, s, dt)
    interaction_kick!(ss, s, dt / 2)
    linear_drift!(ss, s, dt / 2)
    return nothing
end

# ----------------------------------------------------------------------------------------
# System energy
# ----------------------------------------------------------------------------------------
function energy(ss::SimulationState, s::Simulation)
    gm_total = total_gravity_parameter(s)
    gm_0 = gravity_parameter(s, 0)
    gp_0 = HamiltonianState(gm_0, primary(ss))
    h_0 = 0.0 # norm(gp_0)^2 / 2gm_total

    for (i, state) in enumerate(secondaries(ss))
        gm_i = gravity_parameter(s, i)
        p_i = HamiltonianState(gm_i, state)
        p_mag_i = norm(p_i)
        h_0 += p_mag_i^2 / 2gm_0

        h_0 += p_mag_i^2 / 2gm_i -
               gm_0 * gm_i / norm(coordinates(state)) / GRAVITATIONAL_CONSTANT

        if i != length(secondaries(ss))
            for j in 2:length(secondaries(ss))
                x_j = secondaries(ss)[j]
                gm_j = gravity_parameter(s, j)
                h_0 -= gm_i * gm_j / norm(coordinates(state) - coordinates(x_j))
            end
        end
    end

    return h_0
end

end