# ========================================================================================
# File: simulation.jl
# Brief: specification of the simulation
# Author: Rolfe Power
# ========================================================================================

struct Sim{T}
    primary::Body
    secondaries::Vector{Body}
    dimension_set::DimensionSet
    times::T
    primary_states::Vector{HamiltonianState}
    secondary_states::Matrix{HamiltonianState}
    current_index::Int32
end

# ----------------------------------------------------------------------------------------
# N-Body System Structure Definitions and Constructors
# ----------------------------------------------------------------------------------------
struct NBodySystem
    primary::Body
    secondaries::Vector{Body}
end

# ----------------------------------------------------------------------------------------
# N-Body System Accessors
# ----------------------------------------------------------------------------------------
primary_body(nb::NBodySystem) = nb.primary
secondary_bodies(nb::NBodySystem) = nb.secondaries
body_count(nb::NBodySystem) = 1 + length(secondary_bodies(nb))

# ----------------------------------------------------------------------------------------
# N-Body System Derived Values
# ----------------------------------------------------------------------------------------
function total_gravity_parameter(nb::NBodySystem)
    gm0 = gravity_parameter(primary_body(nb))
    gms = sum(gravity_parameter, secondary_bodies(nb))
    return gm0 + gms
end

function total_mass(nb::NBodySystem)
    return total_gravity_parameter(nb) / GRAVITATIONAL_CONSTANT
end

# ----------------------------------------------------------------------------------------
# Simulation Structure Definition and Constructors
# ----------------------------------------------------------------------------------------
struct Simulation{T}
    system::NBodySystem
    time_steps::T
    dimension_set::DimensionSet
end

function Simulation(sys::NBodySystem, tspan; dimensions=DimensionSet())
    return Simulation(sys, tspan, dimensions)
end

# ----------------------------------------------------------------------------------------
# Simulation Structure Accessors
# ----------------------------------------------------------------------------------------
n_body_system(s::Simulation) = s.system
time_steps(s::Simulation) = s.time_steps

# ----------------------------------------------------------------------------------------
# Simulation Initialization
# ----------------------------------------------------------------------------------------
mutable struct SimulationSolution{T}
    sim::Simulation{T}
    times::Vector{Float64}
    primary_states::Vector{HamiltonianState}
    secondary_States::Matrix{HamiltonianState}
    current_index::Int32
end

simulation(ss::SimulationSolution) = ss.sim
current_index(ss::SimulationSolution) = ss.current_index

function init(s::Simulation, x0::Vector{HamiltonianState})
    t = time_steps(s)
    n = length(t)
    n_bodies = body_count(n_body_system(s))

    pri_states = Vector{HamiltonianState}(undef, n)
    sec_states = Matrix{HamiltonianState}(undef, n, n_bodies - 1)

    if length(x0) != n_bodies
        throw(ArgumentError("invalid number of bodies"))
    end

    pri_states[1] = x0[1]
    sec_states[1, :] .= x0[2:end]

    return SimulationSolution(s, t, pri_states, sec_states, 1)
end

# ----------------------------------------------------------------------------------------
# Dynamics
# ----------------------------------------------------------------------------------------
function linear_drift!(ss::SimulationSolution, dt)
    ci = current_index(ss)
end

function interaction_kick!(ss::SimulationSolution, dt)
end

function kepler_drift!(ss::SimulationSolution, dt)
end
