function democratic_heliocentric(s::NBodySystem, pstate::HamiltonianState,
                                 sstates::Vector{HamiltonianState})
    gm_total = total_gravity_parameter(s)
end

function democratic_heliocentric(bodies, states::Vector{HamiltonianState})
    out = similar(states)

    q_dem_i = Vector{SVector{3, Float64}}(undef, length(states))
    p_dem_i = Vector{SVector{3, Float64}}(undef, length(states))

    masses = mass.(bodies)
    total_mass = sum(masses)

    q_dem_i[1] = SVector{3}(0.0, 0.0, 0.0)
    for ind in eachindex(states)
        q_dem_i[1] += masses[ind] * coordinates(states[ind])
    end

    q_dem_i[1] /= total_mass
    p_dem_i[1] = sum(momenta, states)

    out[1] = HamiltonianState(q_dem_i[1], p_dem_i[1])

    for ind in 2:length(states)
        q_dem_i[ind] = coordinates(states[ind]) - coordinates(states[1])
        p_dem_i[ind] = momenta(states[ind]) - mass(bodies[ind]) / total_mass * p_dem_i[1]

        out[ind] = HamiltonianState(q_dem_i[ind], p_dem_i[ind])
    end

    return out
end

function democratic_heliocentric(bodies, states::Vector{LagrangianState})
    return democratic_heliocentric(bodies, HamiltonianState.(bodies, states))
end

function from_democratic_heliocentric(bodies::Vector{Body},
                                      b_states::Vector{HamiltonianState})
    m_total = sum(gravitational_parameter, bodies) / GRAVITATIONAL_CONSTANT

    q_out = fill(SVector{3}(0.0, 0.0, 0.0), length(b_states))
    p_out = fill(SVector{3}(0.0, 0.0, 0.0), length(b_states))

    masses = mass.(bodies)

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
