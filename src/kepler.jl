# ========================================================================================
# File   : src/kepler.jl
# Brief  : Solution of Kepler problem and kepler propagation
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================

# ----------------------------------------------------------------------------------------
# Solving Kepler's problem
# ----------------------------------------------------------------------------------------
"""
    kepler_initial_guess(eccentricity, mean_anomaly)

Return an initial guess to the eccentric anomaly given the mean anomaly and eccentricity.

## References

  - Danby, J. M. A., and Burkardt, T. M. “The Solution of Kepler’s Equation, I.” Celestial
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

  - Danby, J. M. A. “The Solution of Kepler's Equation, III.” Celestial Mechanics, Vol. 40,
    Nos. 3-4, 1987, pp. 303-312. https://doi.org/10.1007/BF01235847.
  - Danby, J. M. A., and Burkardt, T. M. “The Solution of Kepler's Equation, I.” Celestial
    Mechanics, Vol. 31, No. 2, 1983, pp. 95-107. https://doi.org/10.1007/BF01686811.
  - Burkardt, T. M., and Danby, J. M. A. “The Solution of Kepler's Equation, II.” Celestial
    Mechanics, Vol. 31, No. 3, 1983, pp. 317-328. https://doi.org/10.1007/BF01844230.
"""
function kepler_danby(eccentricity, mean_anomaly; tolerance=1e-12, max_iter=10)
    ea = _kepler_initial_guess(eccentricity, mean_anomaly)

    (f, f1, f2, f3) = _kepler_error(eccentricity, mean_anomaly, ea)

    iter_count = 0
    while abs(f) > tolerance && iter_count < max_iter
        d1 = -f / f1
        d2 = -f / (f1 + (1//2) * d1 * f2)
        d3 = -f / (f1 + (1//2) * d2 * f2 + (1//6) * d2^2 * f3)

        ea += d3
        iter_count += 1

        (f, f1, f2, f3) = _kepler_error(eccentricity, mean_anomaly, ea)
    end

    return (ea, abs(f) <= tolerance)
end

# ----------------------------------------------------------------------------------------
# Kepler Propagation
# ----------------------------------------------------------------------------------------
function kepler_propagate(gm, dt, x0::State)
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

    return State(pos=rv, vel=vv)
end
