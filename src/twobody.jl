# ========================================================================================
# File: twobody.jl
# Brief: Basic two-body routines
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

function circular_speed(gm, r)
    return sqrt(gm / r)
end

# ----------------------------------------------------------------------------------------
# Specific Energy
# ----------------------------------------------------------------------------------------
function specific_energy(gm, a::AbstractFloat)
    return -gm / 2a
end

function specific_energy(gm, r::AbstractFloat, v::AbstractFloat)
    return (1 // 2) * v^2 - gm / r
end

function specific_energy(gm, pos::AbstractVector, vel::AbstractVector)
    return specific_energy(gm, norm(pos), norm(vel))
end

function specific_energy(gm, x::LagrangianState)
    return specific_energy(gm, coordinates(x), velocities(x))
end

function specific_energy(b::Body, args...)
    return specific_energy(gravity_parameter(b), args...)
end

# ----------------------------------------------------------------------------------------
# Semi-Major Axis
# ----------------------------------------------------------------------------------------
function semi_major_axis(gm, spec_energy)
    return -gm / 2spec_energy
end

# ----------------------------------------------------------------------------------------
# Semi-Minor Axis
# ----------------------------------------------------------------------------------------
function semi_minor_axis(sma, ecc)
    return sma * sqrt(1 - ecc^2)
end

# ----------------------------------------------------------------------------------------
# Semi-Latus Rectum
# ----------------------------------------------------------------------------------------
function semi_latus_rectum(sma, ecc)
    return sma * (1 - ecc^2)
end

# ----------------------------------------------------------------------------------------
# Radius
# ----------------------------------------------------------------------------------------
function radius(sma, ecc, ta)
    p = semi_latus_rectum(sma, ecc)
    return p / (1 + ecc * cos(ta))
end

# ----------------------------------------------------------------------------------------
# Mean motion
# ----------------------------------------------------------------------------------------
function mean_motion(gm, sma)
    return sqrt(gm / sma^3)
end

# ----------------------------------------------------------------------------------------
# Eccentric Anomaly
# ----------------------------------------------------------------------------------------
function eccentric_anomaly(ecc, ta)
    ta_bound = mod(ta, 2π)
    alpha = sqrt((1 + ecc) / (1 - ecc))
    tane2 = tan(ta_bound / 2) / alpha
    return mod(2atan(tane2), 2π)
end

# ----------------------------------------------------------------------------------------
# Mean anomaly
# ----------------------------------------------------------------------------------------
function mean_anomaly(ecc, eccentric_anomaly)
    return eccentric_anomaly - ecc * sin(eccentric_anomaly)
end

# ----------------------------------------------------------------------------------------
# Period
# ----------------------------------------------------------------------------------------
function period(gm, sma)
    return 2π / mean_motion(gm, sma)
end

# ----------------------------------------------------------------------------------------
# Inclination
# ----------------------------------------------------------------------------------------
function inclination(angular_momentum::AbstractVector)
    return acos(angular_momentum[3] / norm(angular_momentum))
end
