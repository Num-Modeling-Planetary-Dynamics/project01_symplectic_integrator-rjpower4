# ========================================================================================
# File   : src/simulation.jl
# Brief  : Definition of main simulation type
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================

# ----------------------------------------------------------------------------------------
# Simulation Type Definition and Constructors
# ----------------------------------------------------------------------------------------
struct Simulation
    primary::Body
    secondaries::Vector{Body}
    gm_primary::Float64
    gm_secondaries::Vector{Float64}
    unit_set::UnitSet
end

function Simulation(primary, secondaries; units=UnitSet())
    gm_scale = gm_unit(units)
    gm0 = gravity_parameter(primary) / gm_scale
    gmi = @. gravity_parameter(secondaries) / gm_scale
    return Simulation(primary, secondaries, gm0, gmi, units)
end

# ----------------------------------------------------------------------------------------
# Accessors and Basic Derived Quantities
# ----------------------------------------------------------------------------------------
secondary_body_count(s::Simulation) = length(s.secondaries)

gravity_parameter(s::Simulation) = s.gm_primary
gravity_parameter(s::Simulation, index) = s.gm_secondaries[index]
gravity_parameters(s::Simulation) = (s.gm_primary, s.gm_secondaries)

function total_gravity_parameter(s::Simulation)
    let (gm0, gms) = gravity_parameters(s)
        return gm0 + sum(gms)
    end
end

unit_set(s::Simulation) = s.unit_set

# ----------------------------------------------------------------------------------------
# Convert Heliocentric Coordinates to Barycentric
# ----------------------------------------------------------------------------------------
function helio_to_bary(s::Simulation, x0_helio, xi_helio)
    gm_tot = total_gravity_parameter(s)
    _, gmi = gravity_parameters(s)

    x0_bary = -sum(zip(gmi, xi_helio)) do (gm, xi)
        gm * xi
    end ./ gm_tot

    xi_bary = [xi + x0_bary for xi in xi_helio]

    return x0_bary, xi_bary
end

function bary_to_helio(s::Simulation, x0_bary, xi_bary)
    x0_helio = x0_bary - x0_bary
    xi_helio = [xi - x0_bary for xi in xi_bary]
    return (x0_helio, xi_helio)
end

# ----------------------------------------------------------------------------------------
# Simulate
# ----------------------------------------------------------------------------------------
function linear_drift!(xi_barys, s::Simulation, dt)
    gm_tot = total_gravity_parameter(s)
    _, gmi = gravity_parameters(s)
    pt = sum(zip(gmi, xi_barys)) do (gm, xi)
        v = velocities(xi)
        gm * v
    end ./ gm_tot

    for k in eachindex(xi_barys)
        xi_barys[k] += State(; pos=dt * pt)
    end

    return nothing
end

function kepler_drift!(xi_barys, s::Simulation, dt)
    gm0 = gravity_parameter(s)
    for k in eachindex(xi_barys)
        xi_barys[k] = kepler_propagate(gm0, dt, xi_barys[k])
    end
    return nothing
end

function interaction_kick!(xi_barys, s::Simulation, dt)
    coords = coordinates.(xi_barys)
    for k in eachindex(coords)
        for j in eachindex(coords)
            if k == j
                continue
            end

            dr = coords[j] - coords[k]
            dr_mag = norm(dr)
            irij3 = dr_mag^3
            xi_barys[k] += State(; vel=dr * dt * gravity_parameter(s, j) / irij3)
        end
    end
end

function simulate(s::Simulation, dt, x0_helio::State, xi_helio::Vector{State}; n_steps)
    dt_scaled = dt / time_unit(s)

    x0_helios = Vector{State}(undef, n_steps)
    xi_helios = Matrix{State}(undef, n_steps, secondary_body_count(s))

    x0_barys = Vector{State}(undef, n_steps)
    xi_barys = Matrix{State}(undef, n_steps, secondary_body_count(s))

    x_scale = state_unit(s)

    x0_helios[1] = x0_helio ./ x_scale
    xi_helios[1, :] .= [xi ./ x_scale for xi in xi_helio]

    x0_barys[1], xi_barys[1, :] = helio_to_bary(s, x0_helios[1], xi_helios[1, :])

    for k in 2:n_steps
        x0_barys[k] = x0_barys[k - 1]
        xi_barys[k, :] = xi_barys[k - 1, :]

        linear_drift!(selectdim(xi_barys, 1, k), s, dt_scaled / 2)
        interaction_kick!(selectdim(xi_barys, 1, k), s, dt_scaled / 2)
        kepler_drift!(selectdim(xi_barys, 1, k), s, dt_scaled)
        interaction_kick!(selectdim(xi_barys, 1, k), s, dt_scaled / 2)
        linear_drift!(selectdim(xi_barys, 1, k), s, dt_scaled / 2)

        x0_helios[k], xi_helios[k, :] = bary_to_helio(s, x0_barys[k], xi_barys[k, :])
    end

    return x0_helios, xi_helios, x0_barys, xi_barys
end