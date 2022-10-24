mutable struct Simulation
    time_span::Tuple{Float64, Float64}
    time_step::Float64
    bodies::Vector{Body}
    states::Matrix{HamiltonianState}
    step_index::Int32
end

function Simulation(tspan, bodies; max_step)
    tdiff = tspan[2] - tspan[1]
    step = tdiff / ceil(tdiff / max_step)
    n_steps = Int(ceil(tdiff / step)) + 1

    return Simulation(tspan,
                      step,
                      bodies,
                      Matrix{HamiltonianState}(undef, n_steps, length(bodies)),
                      -1)
end

Base.length(s::Simulation) = length(s.states)
bodies(s::Simulation) = s.bodies
time_step(s::Simulation) = s.time_step
body_count(s::Simulation) = length(s.bodies)
step_count(s::Simulation) = s.step_index
time_span(s::Simulation) = b.time_span
total_gm(s::Simulation) = sum(gravitational_parameter, bodies(s))
total_mass(s::Simulation) = total_gm(s) / GRAVITATIONAL_CONSTANT



function init!(s::Simulation, x0::Vector{HamiltonianState})
    s.step_index = 1
    s.states[1, :] .= x0

    return nothing
end

"""
    linear_drift!(::Simulation, dt)

Perform the linear drift part of the drift-kick-drift method for the given time step.

This is the component of the Hamiltonian associated with the linear motion of the system
barycenter. Here, we drift the bodies by the velocity of the barycenter times the time
step.
"""
function linear_drift!(s::Simulation, dt)
    sind = s.step_index
    m0 = mass(bodies(s)[1])
    z = SVector{3}(0.0, 0.0, 0.0)
    for k in 2:body_count(s)
        cs = s.states[sind, k]
        mom = momenta(cs)
        s.states[sind, k] += dt * HamiltonianState(mom / m0, z)
    end
end

function interaction_kick!(s::Simulation, dt)
    z = SVector{3}(0.0, 0.0, 0.0)
    sind = s.step_index

    for i in 2:(body_count(s) - 1)
        xi = s.states[sind, i]
        ci = coordinates(xi)
        gmi = gravitational_parameter(bodies(s)[i])
        dp = z

        for j in (i + 1):body_count(s)
            xj = s.states[sind, j]
            cj = coordinates(xj)
            mj = mass(bodies(s)[j])

            cdiff = (cj - ci)
            cdiffmag = norm(cdiff)

            dp += gmi * mj / cdiffmag^3 * cdiff
            # dp += gravitational_parameter(bodies(s)[j]) / cdiffmag^3 * cdiff
        end

        println("A: $(s.states[sind, i])")
        s.states[sind, i] -= dt * HamiltonianState(z, dp)
        println("Update: $(dt * HamiltonianState(z, dp))")
        println("B: $(s.states[sind, i])")
    end
end

function kepler_drift!(s::Simulation, dt)
    sind = s.step_index
    gm0 = gravitational_parameter(bodies(s)[1])
    for k in 2:body_count(s)
        b = bodies(s)[k]
        le = LagrangianState(b, s.states[sind, k])
        le2 = kepler_propagate(gm0, dt, le)
        s.states[sind, k] = HamiltonianState(mass(b), le2)
    end
end

function do_step!(s::Simulation)
    if s.step_index == size(s.states)[1]
        return -1
    end

    s.states[s.step_index + 1, :] .= s.states[s.step_index, :]
    s.step_index += 1

    dt = time_step(s)

    # linear_drift!(s, dt / 2)
    # interaction_kick!(s, dt / 2)
    kepler_drift!(s, dt)
    # interaction_kick!(s, dt / 2)
    # linear_drift!(s, dt / 2)

    s.step_index
end

# ----------------------------------------------------------------------------------------
# Energy
# ----------------------------------------------------------------------------------------
function system_energy(s::Simulation)
    gm0 = gravitational_parameter(bodies(s)[1])
    gms = gravitational_parameter.(bodies(s))
    masses = mass.(bodies(s))
    n_time_steps = size(s.states, 1)
    m_tot = total_mass(s)
    coords = coordinates.(s.states)
    moms = momenta.(s.states)

    moms_2 = zeros(size(moms))
    for k in eachindex(moms)
        moms_2[k] = sum(x -> x^2, moms[k])
    end

    out = zeros(n_time_steps)

    for k in 1:n_time_steps
        # Hsun
        out[k] += sum(moms_2[k, 2:end]) / (2masses[1])

        for j in 2:length(masses)
            # Hkep
            out[k] += moms_2[k, j] / (2masses[j]) - masses[j] * gm0 / norm(coords[k, j])

            if j < length(masses)
                for i in (j + 1):length(masses)

                    # Hint
                    out[k] -= gms[j] * masses[i] / norm(coords[k, i] - coords[k, j])
                end
            end
        end
    end

    out
end

end
