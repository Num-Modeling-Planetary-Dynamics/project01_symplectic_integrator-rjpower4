# ========================================================================================
# File: simulation.jl
# Brief: specification of the simulation
# Author: Rolfe Power
# ========================================================================================

mutable struct Simulation
    time::Float64
    step_size::Float64
    primary::Body
    secondaries::Vector{Body}
    dimension_set::DimensionSet
    primary_state::HamiltonianState
    secondary_states::Vector{HamiltonianState}
end

function Simulation(dt, primary, secondaries; dimensions = DimensionSet())
    z_state = HamiltonianState(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    Simulation(NaN,
               dt,
               primary,
               secondaries,
               dimensions,
               z_state,
               Vector{HamiltonianState}(undef, length(secondaries)))
end

primary(s::Simulation) = s.primary
secondaries(s::Simulation) = s.secondaries
secondary_states(s::Simulation) = s.secondary_states
secondary_body_count(s::Simulation) = length(secondaries(s))
body_count(s::Simulation) = 1 + secondary_body_count(s)
current_index(s::Simulation) = s.current_index
current_time(s::Simulation) = s.time
step_size(s::Simulation) = s.step_size

function total_gravity_parameter(s::Simulation)
    return gravity_parameter(primary(s)) + sum(gravity_parameter, secondaries(s))
end

total_mass(s::Simulation) = total_gravity_parameter(s) / GRAVITATIONAL_CONSTANT

function init!(s::Simulation, t0, x0::HamiltonianState, xn::Vector{HamiltonianState})
    s.primary_state = x0
    s.secondary_states .= xn
    s.time = t0
    return current_time(s)
end

# ----------------------------------------------------------------------------------------
# Dynamics
# ----------------------------------------------------------------------------------------
function linear_drift!(s::Simulation, dt)
    gm0 = gravity_parameter(primary(s))
    update = Symple.GRAVITATIONAL_CONSTANT * sum(momenta, secondary_states(s))
    s.primary_state += HamiltonianState(dt * update / gm0, SVector{3}(0.0, 0.0, 0.0))
end

function interaction_kick!(s::Simulation, dt)
    z = SVector{3}(0.0, 0.0, 0.0)

    function gm_coord(s, k)
        gmk = secondaries(s)[k] |> gravity_parameter
        qk = secondary_states(s)[k] |> coordinates
        return (gmk, qk)
    end

    for i in 1:(secondary_body_count(s) - 1)
        gmi, qi = gm_coord(s, i)
        dp = z

        for j in (i+1):secondary_body_count(s)
            gmj, qj = gm_coord(s, j)
            dq = qj - qi
            dp += gmi * gmj / norm(dq)^3 * dq / GRAVITATIONAL_CONSTANT
        end

        s.secondary_states[i] += dt * HamiltonianState(z, dp)
    end
end

function kepler_drift!(s::Simulation, dt)
    m0 = gravity_parameter(primary(s))
    for (i, (body, state)) in enumerate(zip(secondaries(s), secondary_states(s)))
        # @show state
        le_init = LagrangianState(gravity_parameter(body) / GRAVITATIONAL_CONSTANT, state)
        le_final =  kepler_propagate(m0, dt, le_init)
        s.secondary_states[i] = HamiltonianState(gravity_parameter(body) / GRAVITATIONAL_CONSTANT, le_final)
    end
end

# ----------------------------------------------------------------------------------------
# Stepping
# ----------------------------------------------------------------------------------------
function do_step!(s::Simulation)
    dt = step_size(s)

    linear_drift!(s, dt/2)
    # interaction_kick!(s, dt/2)
    kepler_drift!(s, dt)
    # interaction_kick!(s, dt/2)
    linear_drift!(s, dt/2)

    return nothing
end

# ----------------------------------------------------------------------------------------
# Energy
# ----------------------------------------------------------------------------------------
function energy(s::Simulation)
    gm0 = gravity_parameter(primary(s))
    total_momenta = sum(momenta, secondary_states(s))
    h_sun = norm(total_momenta)^2 / 2gm0 / GRAVITATIONAL_CONSTANT

    h_kep = 0.0
    h_int = 0.0

    for k in 1:secondary_body_count(s)
        gm_i = secondaries(s)[k] |> gravity_parameter
        q_i = coordinates(secondary_states(s)[k])
        p_i = momenta(secondary_states(s)[k])

        h_kep += norm(p_i)^2 / (2gm_i) / GRAVITATIONAL_CONSTANT - gm0 * gm_i / norm(q_i) / GRAVITATIONAL_CONSTANT

        if k != secondary_body_count(s)
            for j in (k + 1):secondary_body_count(s)
                gm_j = secondaries(s)[j] |> gravity_parameter
                q_j = coordinates(secondary_states(s)[j])
                h_int -= gm_i * gm_j / norm(q_i - q_j) / GRAVITATIONAL_CONSTANT 
            end
        end
    end

    return h_kep + h_int + h_sun + h_sun
end