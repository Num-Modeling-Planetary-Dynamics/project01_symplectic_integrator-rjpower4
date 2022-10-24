# ========================================================================================
# File: state.jl
# Brief: state representations
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

# ----------------------------------------------------------------------------------------
# Base state
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

# ----------------------------------------------------------------------------------------
# Lagrangian State
# ----------------------------------------------------------------------------------------
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

"""
    LagrangianState(r::AbstractVector, v::AbstractVector)

Construct a Lagrangian state from the position and velocity vectors
"""
function LagrangianState(r::AbstractVector, v::AbstractVector)
    return LagrangianState(r[1], r[2], r[3], v[1], v[2], v[3])
end

"""
    velocities(ls::LagrangianState)

Return the velocity measure numbers of the state.
"""
velocities(ls::LagrangianState) = SVector{3}(ls[4], ls[5], ls[6])

"""
    radial_velocity(ls::LagrangianState)

Return the signed scalar radial velocity of the state.
"""
function radial_velocity(ls::LagrangianState)
    pos = coordinates(ls)
    vel = velocities(ls)
    return dot(pos, vel) / norm(pos)
end

"""
    specific_angular_momentum(ls::LagrangianState)

Return the specific angular momentum.
"""
function specific_angular_momentum(ls::LagrangianState)
    pos = coordinates(ls)
    vel = coordinates(ls)
    return cross(pos, vel)
end

# ----------------------------------------------------------------------------------------
# Hamiltonian State
# ----------------------------------------------------------------------------------------
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

"""
    HamiltonianState(r::AbstractVector, p::AbstractVector)

Construct a Hamiltonian state from the position and momenta vectors
"""
function HamiltonianState(r::AbstractVector, p::AbstractVector)
    return HamiltonianState(r[1], r[2], r[3], p[1], p[2], p[3])
end

"""
    momenta(::HamiltonianState)

Return the momenta components of the Hamiltonian state.
"""
function momenta(hs::HamiltonianState)
    return SVector{3}(hs[4], hs[5], hs[6])
end

# ----------------------------------------------------------------------------------------
# Conversion between Lagrangian and Hamiltonian states
# ----------------------------------------------------------------------------------------
"""
    LagrangianState(m, s::HamiltonianState)

Convert the Hamiltonian state to a Lagrangian state.
"""
function LagrangianState(m, s::HamiltonianState)
    return LagrangianState(coordinates(s), momenta(s) / m)
end

"""
    HamiltonianState(m, s::HamiltonianState)

Convert the Lagrangian state to a Hamiltonian state.
"""
function HamiltonianState(m, s::LagrangianState)
    return HamiltonianState(coordinates(s), velocities(s) * m)
end
