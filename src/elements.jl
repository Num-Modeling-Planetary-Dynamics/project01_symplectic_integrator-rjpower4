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

function specific_energy(gm, x::State)
    return specific_energy(gm, coordinates(x), velocities(x))
end

function specific_energy(b::Body, x::State)
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

    return SMatrix{3,3}(r11, r21, r31, r12, r22, r32, r13, r23, r33)
end

"""
    KeplerianElements

Collection of Keplerian elements
"""
struct KeplerianElements <: FieldVector{6,Float64}
    sma::Float64
    ecc::Float64
    inc::Float64
    aop::Float64
    raan::Float64
    ta::Float64
end

function KeplerianElements(sma; ecc=0.0, inc=0.0, aop=0.0, raan=0.0, ta=0.0)
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

function eccentric_anomaly(ke::KeplerianElements; ta=true_anomaly(ke))
    return eccentric_anomaly(eccentricity(ke), ta)
end

function mean_anomaly(ke::KeplerianElements; ta=true_anomaly(ke)) 
    ea = eccentric_anomaly(ke; ta=ta)
    ma = mean_anomaly(eccentricity(ke), ea)
end

function longitude_of_periapsis(ke::KeplerianElements)
    return argument_of_periapsis(ke) + right_ascension(ke)
end

function radius(ke::KeplerianElements; ta=true_anomaly(ke))
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

    sma = -gm / (2energy)
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

function KeplerianElements(gm, x::State)
    return KeplerianElements(gm, coordinates(x), velocities(x))
end

function State(gm, ke::KeplerianElements)
    z = zero(eltype(ke))
    rot_mat = perifocal_to_inertial_rotation(
        inclination(ke), argument_of_periapsis(ke) + true_anomaly(ke), right_ascension(ke)
    )
    ta = true_anomaly(ke)
    p = semi_latus_rectum(ke)
    h = sqrt(gm * p)

    r_peri = SVector{3}(radius(ke), z, z)
    r_inert = rot_mat * r_peri

    v_transverse = h / norm(r_peri)
    v_radial = h / p * eccentricity(ke) * sin(ta)
    v_peri = SVector{3}(v_radial, v_transverse, z)
    v_inert = rot_mat * v_peri

    return State(; pos=r_inert, vel=v_inert)
end
