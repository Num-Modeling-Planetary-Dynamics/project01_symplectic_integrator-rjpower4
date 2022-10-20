# ========================================================================================
# File: elements.jl
# Brief: Keplerian elements
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

# ----------------------------------------------------------------------------------------
# Type definition and external constructors
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
    ecc > 0.0 || throw(ArgumentError("invalid eccentricity: $ecc"))

    if sma .* (ecc - 1.0) > 0.0
        throw(ArgumentError("invalid semi-major axis and eccentricity"))
    elseif ecc == 1.0 && !isinf(sma)
        throw(ArgumentError("invalid semi-major axis and eccentricity"))
    end

    return KeplerianElements(sma, ecc, inc, aop, raan, ta)
end

# ----------------------------------------------------------------------------------------
# Accessors
# ----------------------------------------------------------------------------------------
semi_major_axis(ke::KeplerianElements) = ke.sma
eccentricity(ke::KeplerianElements) = ke.ecc
inclination(ke::KeplerianElements) = ke.inc
argument_of_periapsis(ke::KeplerianElements) = ke.aop
right_ascension(ke::KeplerianElements) = ke.raan
true_anomaly(ke::KeplerianElements) = ke.ta

# ----------------------------------------------------------------------------------------
# Extension of the two-body methods
# ----------------------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------------------
# Conversion from state to Keplerian elements
# ----------------------------------------------------------------------------------------
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

    energy = specific_energy(gm, r, v)
    sma = semi_major_axis(gm, energy)
    h = cross(pos, vel)
    inc = inclination(h)

    ecc_vector = (1/gm) * (v^2 * pos - dot(r, v) * vel) - rhat
    ecc = norm(ecc_vector)
    e_hat = ecc != 0.0  ? ecc_vector ./ ecc : unit

    nodal_axis = SVector{3}(-h[2], h[1], 0.0)
    nodal_axis_hat = h[2] == 0 ? unit : normalize(nodal_axis)

    raan = wrap_if(acos(nodal_axis_hat[1]), nodal_axis_hat[2] < 0.0)
    aop = wrap_if(angle_between(nodal_axis_hat, e_hat), ecc_vector[3] < 0.0)
    ta = clamp(dot(e_hat, rhat), -1.0, 1.0) |> acos

    return KeplerianElements(sma, ecc, inc, aop, raan, ta)
end

function KeplerianElements(gm, x::LagrangianState)
    return KeplerianElements(gm, coordinates(x), velocities(x))
end

# ----------------------------------------------------------------------------------------
# Conversion from Keplerian elements to State
# ----------------------------------------------------------------------------------------
function LagrangianState(gm, ke::KeplerianElements)
    z = zero(eltype(ke))
    rot_mat = perifocal_to_inertial_rotation(ke)
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

