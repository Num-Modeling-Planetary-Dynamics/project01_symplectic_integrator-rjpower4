# ========================================================================================
# File: symple.jl
# Author: Rolfe Power
# ========================================================================================
module Symple

# ----------------------------------------------------------------------------------------
# External packages
# ----------------------------------------------------------------------------------------
using LinearAlgebra: norm, cross, dot
using StaticArrays: FieldVector, SVector, MVector
using SPICE: furnsh, spkezr, str2et

# ----------------------------------------------------------------------------------------
# Internal Includes
# ----------------------------------------------------------------------------------------
include("constants.jl")
include("dimensions.jl")

# ----------------------------------------------------------------------------------------
# Utility Functions
# ----------------------------------------------------------------------------------------
"""
    angular_difference(a, b)

Return the correctly signed counter-clockwise angular difference from a to b.

## Examples

* `angular_difference(0.0, π/2) = π/2`
* `angular_difference(π/2, 0.0) = 3π/2`
"""
function angular_difference(a, b)
    ay, ax = sincos(a)
    by, bx = sincos(b)
    theta = atan(ax * by - ay * bx, ax * bx + ay * by)
    return theta > 0 ? theta : 2π + theta
end

# ----------------------------------------------------------------------------------------
# Celestial Body Type
# ----------------------------------------------------------------------------------------
"""
    Body

Celestial body abstraction.
"""
struct Body
    name::String
    naif_id::Int32
    gm::Float64

    # Check for following invariants:
    #   1. Name is not empty
    #   2. GM is positive number
    function Body(name, naif_id, gm)
        isempty(name) && throw(ArgumentError("empty body name"))
        (gm <= 0.0) && throw(ArgumentError("non-positive gravitational parameter"))
        return new(name, naif_id, gm)
    end
end

# Accessors
name(b::Body) = b.name
naif_id(b::Body) = b.naif_id
gravitational_parameter(b::Body) = b.gm

# Derived Values
mass(b::Body) = gravitational_parameter(b) / GRAVITATIONAL_CONSTANT

# ----------------------------------------------------------------------------------------
# Pre defined Bodies
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

# ----------------------------------------------------------------------------------------
# State representations
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
LagrangianState(r, v) = LagrangianState(r[1], r[2], r[3], v[1], v[2], v[3])
velocities(ls::LagrangianState) = SVector{3}(ls[4], ls[5], ls[6])

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

HamiltonianState(r, p) = HamiltonianState(r[1], r[2], r[3], p[1], p[2], p[3])
momenta(hs::HamiltonianState) = SVector{3}(hs[4], hs[5], hs[6])

# ----------------------------------------------------------------------------------------
# Hamiltonian <--> Lagrangian Conversions
# ----------------------------------------------------------------------------------------
function LagrangianState(m, s::HamiltonianState)
    return LagrangianState(coordinates(s), momenta(s) / m)
end

function LagrangianState(b::Body, s::HamiltonianState)
    return LagrangianState(mass(b), s)
end

function HamiltonianState(m, s::LagrangianState)
    return HamiltonianState(coordinates(s), velocities(s) * m)
end

function HamiltonianState(b::Body, s::LagrangianState)
    return HamiltonianState(mass(b), s)
end

# ----------------------------------------------------------------------------------------
# Spice Wrapping Functions
# ----------------------------------------------------------------------------------------
"""
    load_kernels(dirpath::String)

Load all of the SPICE kernels in the directory pointed to by `dirpath`.
"""
function load_kernels(dirpath::String)
    valid_endings = ("bsp", "tls", "tpc")
    if !isdir(dirpath)
        return
    end

    for pth in readdir(dirpath; join = true)
        if isfile(pth)
            furnsh(pth)
        end
    end
end

function ephemeris_time(time_string::String)
    return str2et(time_string)
end

function LagrangianState(b::Body, et; frame = "J2000", observer = "SOLAR_SYSTEM_BARYCENTER")
    out_raw, _light_time = spkezr(name(b), et, frame, "NONE", observer)
    return LagrangianState(out_raw)
end

# ----------------------------------------------------------------------------------------
# Two body utilities
# ----------------------------------------------------------------------------------------
function specific_energy(gm, x::LagrangianState)
    pos = coordinates(x)
    vel = velocities(x)
    r = norm(pos)
    v2 = sum(x -> x^2, vel)
    return (1 // 2) * v2 - gm / r
end

# ----------------------------------------------------------------------------------------
# Keplerian Element utilities
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

# Accessor methods
semi_major_axis(ke::KeplerianElements) = ke.sma
eccentricity(ke::KeplerianElements) = ke.ecc
inclination(ke::KeplerianElements) = ke.inc
argument_of_periapsis(ke::KeplerianElements) = ke.aop
right_ascension(ke::KeplerianElements) = ke.right_ascension
true_anomaly(ke::KeplerianElements) = ke.ta

# Derived
function mean_motion(gm, ke::KeplerianElements)
    return sqrt(gm / semi_major_axis(ke)^3)
end

function period(gm, ke::KeplerianElements)
    return 2π / mean_motion(gm, ke)
end

function eccentric_anomaly(gm, ke::KeplerianElements)
    ecc = eccentricity(ke)

    ta = true_anomaly(ke)
    ta = mod(ta, 2π)
    alpha = sqrt((1 + ecc) / (1 - ecc))
    tane2 = tan(ta / 2) / alpha

    return mod(2.0 * atan(tane2), 2π)
end

"""
    shape_elements(gm, x::LagrangianState)

Compute the semi-major axis and eccentricity from the gravitational paramter and state.
"""
function shape_elements(gm, s::LagrangianState)
    r_vec = coordinates(s)
    r_mag = norm(r_vec)
    r_hat = r_vec / r_mag

    v_vec = velocities(s)
    v_mag_squared = sum(x -> x^2, v_vec)

    energy = v_mag_squared / 2 - gm / r_mag
    sma = -gm / (2energy)

    e_vec = (1 / gm) * (v_mag_squared * r_vec - dot(r_vec, v_vec) * v_vec) - r_hat
    ecc = norm(e_vec)

    return sma, ecc
end

inclination(h::AbstractVector) = acos(h[3] / norm(h))

"""
    KeplerianElements(gm, state::LagrangianState)

Construct a new set of Keplerian elements from the state and central body GM value.
"""
function KeplerianElements(gm, state::LagrangianState)
    wrap_if(val, condition) = condition ? 2pi - val : val
    angle_between(v1, v2) = acos(dot(v1, v2))
    normalize(v) = v ./ norm(v)

    r = coordinates(state)
    v = velocities(state)

    r_mag = norm(r)
    r_hat = r / r_mag
    v_mag_squared = sum(x -> x^2, v)

    # Semi-Major Axis
    spec_energy = v_mag_squared / 2 - gm / r_mag
    sma = -gm / (2 * spec_energy)

    # Inclination
    h = cross(r, v)
    inc = inclination(h)

    # Eccentricity
    e_vector = (1 / gm) * (v_mag_squared * r - dot(r, v) * v) - r_hat
    ecc = norm(e_vector)
    e_hat = e_vector ./ ecc

    nodal_axis = (-h[2], h[1], 0.0)
    nodal_axis_hat = h[2] == 0 ? (1.0, 0.0, 0.0) : normalize(nodal_axis)

    # Right Ascension, Argument of Periapsis, and True anomaly
    raan = wrap_if(acos(nodal_axis_hat[1]), nodal_axis_hat[2] < 0.0)
    aop = wrap_if(angle_between(nodal_axis_hat, e_hat), e_vector[3] < 0.0)
    ta = clamp(dot(e_hat, r_hat), -1.0, 1.0) |> acos

    return KeplerianElements(sma, ecc, inc, aop, raan, ta)
end

# ----------------------------------------------------------------------------------------
# Solution to Kepler's Problem
# ----------------------------------------------------------------------------------------

"""
    kepler_initial_guess(eccentricity, mean_anomaly)

Return an initial guess to the eccentric anomaly given the mean anomaly and eccentricity.

## References
* Danby, J. M. A., and Burkardt, T. M. “The Solution of Kepler’s Equation, I.” Celestial
  Mechanics, Vol. 31, No. 2, 1983, pp. 95-107. https://doi.org/10.1007/BF01686811.
"""
function kepler_initial_guess(eccentricity, mean_anomaly)
    return mean_anomaly
end

"""
    _kepler_error(ecc, ma, ea)

Compute error in Kepler problem along with first 3 derivatives w.r.t. eccentric anomaly.
"""
function _kepler_error(ecc, ma, ea)
    ese, ece = ecc .* sincos(ea)
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
    ea = kepler_initial_guess(eccentricity, mean_anomaly)

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

    ea_initial = eccentric_anomaly(gm, els)
    ma_initial = ea_initial - eccentricity(els) * sin(ea_initial)

    ma_final = ma_initial + n * dt
    ea_final, did_converge = kepler_danby(eccentricity(els), ma_final)
    if !did_converge
        error("unable to converge eccentric anomaly")
    end

    de = angular_difference(ea_initial, ea_final)
    sde, cde = sincos(de)

    q0 = coordinates(x0)
    dq0 = velocities(x0)

    r0 = norm(coordinates(x0))

    f = a / r0 * (cde - 1) + 1
    g = dt + (1 / n) * (sde - de)

    new_coordinates = f * q0 + g * dq0
    r = norm(new_coordinates)

    df = -a^2 / (r * r0) * n * sde
    dg = a / r * (cde - 1) + 1

    new_velocities = df * q0 + dg * dq0

    return LagrangianState(new_coordinates, new_velocities)
end

# ----------------------------------------------------------------------------------------
# Democratic Heliocentric
# ----------------------------------------------------------------------------------------
function democratic_heliocentric(bodies, states::Vector{HamiltonianState})
    out = similar(states)

    q_dem_i = Vector{SVector{3, Float64}}(undef, length(states))
    p_dem_i = Vector{SVector{3, Float64}}(undef, length(states))

    masses = mass.(bodies)
    total_mass = sum(masses)

    q_dem_i[1] = SVector{3}(0.0, 0.0, 0.0)
    for ind in eachindex(states)
        q_dem_i[1] += masses[ind] * coordinates(states[ind])
    end

    q_dem_i[1] /= total_mass
    p_dem_i[1] = sum(momenta, states)

    out[1] = HamiltonianState(q_dem_i[1], p_dem_i[1])

    for ind in 2:length(states)
        q_dem_i[ind] = coordinates(states[ind]) - coordinates(states[1])
        p_dem_i[ind] = momenta(states[ind]) - mass(bodies[ind]) / total_mass * p_dem_i[1]

        out[ind] = HamiltonianState(q_dem_i[ind], p_dem_i[ind])
    end

    return out
end

function democratic_heliocentric(bodies, states::Vector{LagrangianState})
    return democratic_heliocentric(bodies, HamiltonianState.(bodies, states))
end

function from_democratic_heliocentric(bodies::Vector{Body},
                                      b_states::Vector{HamiltonianState})
    m_total = sum(gravitational_parameter, bodies) / GRAVITATIONAL_CONSTANT

    q_out = fill(SVector{3}(0.0, 0.0, 0.0), length(b_states))
    p_out = fill(SVector{3}(0.0, 0.0, 0.0), length(b_states))

    masses = mass.(bodies)

    coords = coordinates.(b_states)
    moms = momenta.(b_states)

    mq = SVector{3}(0.0, 0.0, 0.0)
    for k in 2:length(bodies)
        mq += gravitational_parameter(bodies[k]) * coords[k]
    end
    mq /= GRAVITATIONAL_CONSTANT

    out = Vector{HamiltonianState}(undef, length(b_states))

    q_out[1] = coords[1] - mq / m_total

    for k in 2:length(bodies)
        q_out[k] = coords[k] + q_out[1] - mq / m_total
        p_out[k] = moms[k] + masses[k] / m_total * moms[1]

        out[k] = HamiltonianState(q_out[k], p_out[k])
    end

    p_out[1] = moms[1] - sum(p_out)

    out[1] = HamiltonianState(q_out[1], p_out[1])

    return out
end

# ----------------------------------------------------------------------------------------
# Simulation
# ----------------------------------------------------------------------------------------
mutable struct Simulation
    time_span::Tuple{Float64, Float64}
    time_step::Float64
    bodies::Vector{Body}
    states::Matrix{HamiltonianState}
    step_index::Int32
end

function Simulation(tspan, bodies; max_step)
    tdiff = tspan[2] - tspan[1]
    step = tdiff / ceil(tdiff / max_step)
    n_steps = Int(ceil(tdiff / step)) + 1

    return Simulation(tspan,
                      step,
                      bodies,
                      Matrix{HamiltonianState}(undef, n_steps, length(bodies)),
                      -1)
end

Base.length(s::Simulation) = length(s.states)
bodies(s::Simulation) = s.bodies
time_step(s::Simulation) = s.time_step
body_count(s::Simulation) = length(s.bodies)
step_count(s::Simulation) = s.step_index
time_span(s::Simulation) = b.time_span
total_gm(s::Simulation) = sum(gravitational_parameter, bodies(s))
total_mass(s::Simulation) = total_gm(s) / GRAVITATIONAL_CONSTANT

function init!(s::Simulation, x0::Vector{HamiltonianState})
    s.step_index = 1
    s.states[1, :] .= x0

    return nothing
end

function linear_drift!(s::Simulation, dt)
    sind = s.step_index
    m0 = mass(bodies(s)[1])
    z = SVector{3}(0.0, 0.0, 0.0)
    for k in 2:body_count(s)
        cs = s.states[sind, k]
        mom = momenta(cs)
        s.states[sind, k] += dt * HamiltonianState(mom / m0, z)
    end
end

function interaction_kick!(s::Simulation, dt)
    z = SVector{3}(0.0, 0.0, 0.0)
    sind = s.step_index

    for i in 2:(body_count(s) - 1)
        xi = s.states[sind, i]
        ci = coordinates(xi)
        gmi = gravitational_parameter(bodies(s)[i])
        dp = z

        for j in (i + 1):body_count(s)
            xj = s.states[sind, j]
            cj = coordinates(xj)
            mj = mass(bodies(s)[j])

            cdiff = (cj - ci)
            cdiffmag = norm(cdiff)

            dp += -gmi * mj / cdiffmag^3 * cdiff
        end

        s.states[sind, i] -= dt * HamiltonianState(z, dp)
    end
end

function kepler_drift!(s::Simulation, dt)
    sind = s.step_index
    gm0 = gravitational_parameter(bodies(s)[1])
    for k in 2:body_count(s)
        b = bodies(s)[k]
        le = LagrangianState(b, s.states[sind, k])
        le2 = kepler_propagate(gm0, dt, le)
        s.states[sind, k] = HamiltonianState(mass(b), le2)
    end
end

function do_step!(s::Simulation)
    if s.step_index == size(s.states)[1]
        return -1
    end

    s.states[s.step_index + 1, :] .= s.states[s.step_index, :]
    s.step_index += 1

    dt = time_step(s)

    interaction_kick!(s, dt / 2)
    linear_drift!(s, dt / 2)
    kepler_drift!(s, dt)
    linear_drift!(s, dt / 2)
    interaction_kick!(s, dt / 2)

    s.step_index
end

# ----------------------------------------------------------------------------------------
# Energy
# ----------------------------------------------------------------------------------------
function system_energy(s::Simulation)
    gm0 = gravitational_parameter(bodies(s)[1])
    gms = gravitational_parameter.(bodies(s))
    masses = mass.(bodies(s))
    n_time_steps = size(s.states, 1)
    m_tot = total_mass(s)
    coords = coordinates.(s.states)
    moms = momenta.(s.states)

    moms_2 = zeros(size(moms))
    for k in eachindex(moms)
        moms_2[k] = sum(x -> x^2, moms[k])
    end

    out = zeros(n_time_steps)

    for k in 1:n_time_steps
        # Hsun
        out[k] += sum(moms_2[k, 2:end]) / (2masses[1])

        for j in 2:length(masses)
            # Hkep
            out[k] += moms_2[k, j] / (2masses[j]) - masses[j] * gm0 / norm(coords[k, j])

            if j < length(masses)
                for i in (j + 1):length(masses)

                    # Hint
                    out[k] -= gms[j] * masses[i] / norm(coords[k,i] - coords[k,j])
                end
            end
        end
    end

    out
end

end
