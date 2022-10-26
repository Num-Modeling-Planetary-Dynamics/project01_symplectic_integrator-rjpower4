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

export SUN, SOLAR_SYSTEM_BARYCENTER
export EARTH_BARYCENTER, MARS_BARYCENTER, JUPITER_BARYCENTER, SATURN_BARYCENTER
export URANUS_BARYCENTER, NEPTUNE_BARYCENTER, PLUTO_BARYCENTER

# ----------------------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------------------
include("constants.jl")

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
    # pos = SVector{3}(posv)
    # vel = SVector{3}(velv)
    # unit = SVector{3}(one(eltype(pos)), zero(eltype(pos)), zero(eltype(pos)))

    # normalize(v) = v ./ norm(v)
    # angle_between(v1, v2) = acos(clamp(dot(v1, v2), -1, 1))
    # wrap_if(val, condition) = condition ? 2pi - val : val

    # r = norm(pos)
    # v = norm(vel)

    # rhat = pos ./ r

    # sma = 1 / (2 / r - v^2 / gm)
    # h = cross(pos, vel)
    # hmag = norm(h)
    # ecc = sqrt(1 - hmag^2 / (gm * sma))
    # inc = inclination(h)

    # ecc_vector = cross(vel, h) / gm - rhat
    # e_hat = ecc != 0.0 ? ecc_vector ./ ecc : unit

    # nodal_axis = SVector{3}(-h[2], h[1], 0.0)
    # nodal_axis_hat = nodal_axis ./ norm(nodal_axis)

    # raan = wrap_if(acos(nodal_axis_hat[1]), nodal_axis_hat[2] < 0.0)
    # aop = wrap_if(angle_between(nodal_axis_hat, e_hat), ecc_vector[3] < 0.0)
    # ta = acos(clamp(dot(e_hat, rhat), -1.0, 1.0))

    # return KeplerianElements(sma, ecc, inc, aop, raan, ta)

    khat = SVector{3}(0.0, 0.0, 1.0)

    rvec = SVector{3}(posv)
    vvec = SVector{3}(velv)
    hvec = cross(rvec, vvec)

    r = norm(rvec)
    v = norm(vvec)
    h = norm(hvec)

    nvec = cross(khat, hvec)
    n = norm(nvec)

    energy = v^2 / 2 - gm / r

    evec = cross(vvec, hvec) / gm - rvec / r
    ecc = norm(evec)

    sma = - gm / (2energy)
    inc = acos(hvec[3] / h)

    raan = let raan_0 = acos(nvec[1] / n)
        if nvec[2] >= 0
            raan_0
        else
            2π - raan_0
        end
    end

    aop = let aop_0 = acos(dot(nvec, evec) / (n * ecc))
        if evec[3] < 0
            2π - aop_0
        else
            aop_0
        end
    end

    ta = let ta_0 = acos(dot(evec, rvec) / (ecc * r))
        if dot(rvec, vvec) >= 0
            ta_0
        else
            2π - ta_0
        end
    end

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

include("kepler.jl")

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
            dx_i += LagrangianState(SVector{3}(0.0, 0.0, 0.0), (gm_i + gm_j) / r^3 * q_ij)
        end
        ss.secondary_states[i] += dt * dx_i
    end
end

function kepler_drift!(ss::SimulationState, s::Simulation, dt)
    states = secondaries(ss)
    gm_0 = gravity_parameter(s, 0)

    for (i, x) in enumerate(states)
        x_new = kepler_propagate(gm_0 + gravity_parameter(s, i), dt, x)
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
    gm_0 = gravity_parameter(s, 0)
    gp_0 = HamiltonianState(gm_0, primary(ss))

    h_0 = 0.0

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
                h_0 -= gm_i * gm_j / norm(coordinates(state) - coordinates(x_j)) / GRAVITATIONAL_CONSTANT
            end
        end
    end

    return h_0
end

end
