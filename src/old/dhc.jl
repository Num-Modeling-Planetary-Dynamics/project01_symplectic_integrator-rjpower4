# ========================================================================================
# File: dhc.jl
# Brief: (De)construction of the democratic Heliocentric coordinates
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

"""
    democratic_heliocentric(s::NBodySystem, x0::HamiltonianState, xn::Vector{HamiltonianState})

Convert Hamiltonian states in the inertial frame to states in the democratic Heliocentric frame.
"""
function democratic_heliocentric(p::Body, s::Vector{Body}, x0::HamiltonianState,
                                 xn::Vector{HamiltonianState})
    n_bodies = 1 + length(s)
    gm_total = gravity_parameter(p) + sum(gravity_parameter, s)

    q_bary = Vector{SVector{3, Float64}}(undef, n_bodies)
    p_bary = Vector{SVector{3, Float64}}(undef, n_bodies)

    q_bary[1] = gravity_parameter(p) * coordinates(x0)
    p_bary[1] = momenta(x0) + sum(momenta, xn)

    for (body, state) in zip(s, xn)
        q_bary[1] += gravity_parameter(body) * coordinates(state)
    end
    q_bary[1] /= gm_total

    for (i, (body, state)) in enumerate(zip(s, xn))
        q_bary[i + 1] = coordinates(state) - q_bary[1]
        p_bary[i + 1] = momenta(state) - gravity_parameter(body) / gm_total * p_bary[1]
    end

    return (HamiltonianState(q_bary[1], p_bary[1]),
            HamiltonianState.(q_bary[2:end], p_bary[2:end]))
end

function from_democratic_heliocentric(p::Body, s::Vector{Body}, x0::HamiltonianState,
                                      xn::Vector{HamiltonianState})
    n_bodies = 1 + length(s)
    gm_total = gravity_parameter(p) + sum(gravity_parameter, s)

    q_secondary = fill(SVector{3}(0.0, 0.0, 0.0), length(s))
    p_out = fill(SVector{3}(0.0, 0.0, 0.0), n_bodies)

    masses = mass.(bodies)

    q_primary = coordinates(x0)
    for (body, state) in zip(s, xn)
        q_primary -= gravity_parameter(body) * coordinates(state)
    end
    for k in eachindex(q_secondary)
        q_secondary[k] = coordinates(xn[k]) - q_primary
    end

    coords = coordinates.(b_states)
    moms = momenta.(b_states)

    mq = SVector{3}(0.0, 0.0, 0.0)
    for k in 2:length(bodies)
        mq += gravitational_parameter(bodies[k]) * coords[k]
    end
    mq /= GRAVITATIONAL_CONSTANT

    out = Vector{HamiltonianState}(undef, length(b_states))

    q_out[1] = coords[1] - mq / m_total

    for k in 2:length(bodies)
        q_out[k] = coords[k] + q_out[1] - mq / m_total
        p_out[k] = moms[k] + masses[k] / m_total * moms[1]

        out[k] = HamiltonianState(q_out[k], p_out[k])
    end

    p_out[1] = moms[1] - sum(p_out)

    out[1] = HamiltonianState(q_out[1], p_out[1])

    return out
end
