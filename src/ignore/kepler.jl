""" kepler.jl

Kepler problem and propagation.
"""

# ----------------------------------------------------------------------------------------
# Kepler Error evaluation
# ----------------------------------------------------------------------------------------
kepler_error(e, M, E) = E - e * sin(E) - M
kepler_error_d1(e, M, E) = 1 - e * cos(E)
kepler_error_d2(e, M, E) = e * sin(E)
kepler_error_d3(e, M, E) = e * cos(E)

# ----------------------------------------------------------------------------------------
# Update Methods
#
#   Currently Implemented:
#       * Simple Iteration
#       * Newton-Raphson
#       * Danby (Murray and Dermott Section 2.4)
# ----------------------------------------------------------------------------------------
abstract type KeplerUpdater end

initial_guess(::KeplerUpdater, ecc, ma) = ma

struct DanbyKeplerUpdater <: KeplerUpdater end

function update(::DanbyKeplerUpdater, e, M, E)
    f = kepler_error(e, M, E)
    fp = kepler_error_d1(e, M, E)
    fpp = kepler_error_d2(e, M, E)
    fppp = kepler_error_d3(e, M, E)

    delta_1 = -f / fp
    delta_2 = -f / (fp + (1 // 2) * delta_1 * fpp)
    delta_3 = -f / (fp + (1 // 2) * delta_2 * fpp + (1 // 6) * delta_2^2 * fppp)

    return E + delta_3
end

struct SimpleKeplerUpdate <: KeplerUpdater end

function update(::SimpleKeplerUpdate, e, M, E)
    return M + e * sin(E)
end

struct NewtonRaphsonKeplerUpdate <: KeplerUpdater end

function update(::NewtonRaphsonKeplerUpdate, e, M, E)
    return E - kepler_error(e, M, E) / kepler_error_d1(e, M, E)
end

# ----------------------------------------------------------------------------------------
# Solve Kepler's Problem: Find the eccentric anomaly associated with a given mean anomaly
# ----------------------------------------------------------------------------------------
function kepler_solve(updater, ecc, ma; tolerance = deg2rad(1e-6), max_iter = Inf)
    guess = initial_guess(updater, ecc, ma)
    err_val = kepler_error(ecc, ma, guess)

    iter_count = 0
    while abs(err_val) > tolerance && iter_count < max_iter
        guess = update(updater, ecc, ma, guess)
        err_val = kepler_error(ecc, ma, guess)
    end

    return guess
end

# ----------------------------------------------------------------------------------------
# Universal Variable Formulation
# ----------------------------------------------------------------------------------------
function stumpff(z)
    c, s = if z > 0
        sz = sqrt(z)
        ((1 - cos(sz)) / z,
         (sz - sin(sz)) / sqrt(z^3))
    elseif z < 0
        sz = sqrt(-z)
        ((cosh(sqrt(-z)) - one(z)) / (-z),
         (sinh(sz) - sz) / sqrt((-z)^3))
    else
        (one(z) / 2, one(z) / 6)
    end

    (s, c)
end
function stumpff_c(z)
    zz = zero(z)
    if z > zz
        return (one(z) - cos(sqrt(z))) / z
    elseif z < zz
        return (cosh(sqrt(-z)) - one(z)) / (-z)
    else
        return one(z) / 2
    end
end

function stumpff_s(z)
    zz = zero(z)
    if z > zz
        sz = sqrt(z)
        return (sz - sin(sz)) / sqrt(z^3)
    elseif z < zz
        sz = sqrt(-z)
        return (sinh(sz) - sz) / sqrt((-z)^3)
    else
        return one(z) / 6
    end
end

_universal_kepler_z(alpha, chi) = alpha * chi^2

function _universal_kepler_error(gm, x0, dt, alpha, chi)
    sgm = sqrt(gm)
    pos, vel = posvel(x0)
    r0 = norm(pos)
    rdv = dot(pos, vel)
    z = _universal_kepler_z(alpha, chi)
    s, c = stumpff(z)

    err = rdv / sgm * chi^2 * c + (1 - alpha * r0) * chi^3 * s + r0 * chi - sgm * dt
    derr = rdv / sgm * chi * (1 - z * s) + (1 - alpha * r0) * chi^2 * c + r0
    return (err, derr)
end

"""
    kepler_propagate(gm, dt, x0)

Propagate the state for the given time duration in the two-body keplerian model.
"""
function kepler_propagate(gm, dt, x0; tol = 1e-6)
    pos0, vel0 = posvel(x0)
    r = norm(pos0)
    v = norm(vel0)
    alpha = (2 / r - v^2 / gm)
    chi = sqrt(gm) * abs(alpha) * dt

    err, derr = _universal_kepler_error(gm, x0, dt, alpha, chi)
    while abs(err) > tol
        chi -= err / derr
        err, derr = _universal_kepler_error(gm, x0, dt, alpha, chi)
    end

    z = _universal_kepler_z(alpha, chi)
    s, c = stumpff(z)

    f = 1 - chi^2 / r * c
    g = dt  - chi^3 / sqrt(gm) * s

    pos = f * pos0 + g * vel0

    df = sqrt(gm) / (r * norm(pos)) * chi * (z * s - 1)
    dg = 1 - chi^2 / norm(pos) * c

    vel = df * pos0 + dg * vel0


    return State(pos..., vel...)
end
