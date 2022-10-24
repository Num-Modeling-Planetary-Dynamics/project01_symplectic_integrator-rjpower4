# ========================================================================================
# File: kepler.jl
# Brief: solution to Kepler's problem as well as Kepler Propagation
# Author: Rolfe Power
# ========================================================================================

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

    de = angular_difference(ea_initial, ea_final)
    sde, cde = sincos(de)

    q0 = coordinates(x0)
    dq0 = velocities(x0)

    r0 = norm(q0)

    f = a / r0 * (cde - 1) + 1
    g = dt + (1 / n) * (sde - (ea_final - ea_initial))

    new_coordinates = f * q0 + g * dq0
    r = norm(new_coordinates)

    df = -a^2 / (r * r0) * n * sde
    dg = a / r * (cde - 1) + 1

    new_velocities = df * q0 + dg * dq0

    return LagrangianState(new_coordinates, new_velocities)
end

