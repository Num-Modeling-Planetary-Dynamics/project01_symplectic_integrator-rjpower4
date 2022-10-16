module EAPS

using LinearAlgebra: norm, dot
using StaticArrays: FieldVector, SVector, MVector
using SPICE: spkezr, str2et

include("constants.jl")
include("pathutils.jl")
include("spice.jl")
include("state.jl")
include("bodies.jl")
include("twobody.jl")
include("kepler.jl")

# ----------------------------------------------------------------------------------------
# Canonical State
# ----------------------------------------------------------------------------------------
struct CanonicalState <: FieldVector{6, Float64}
    x::Float64
    y::Float64
    z::Float64
    px::Float64
    py::Float64
    pz::Float64
end

coordinates(cs::CanonicalState) = SVector{3}(cs.x, cs.y, cs.z)
momenta(cs::CanonicalState) = SVector{3}(cs.px, cs.py, cs.pz)

function CanonicalState(b::Body, s::State)
    pos, vel = posvel(s)
    m = mass(b)
    return CanonicalState(pos..., (m * vel)...)
end

function State(b::Body, cs::CanonicalState)
    pos = coordinates(cs)
    mom = momenta(cs)
    return State(pos...,
                 (mom / mass(b))...)
end

# ----------------------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------------------
function to_democratic_heliocentric(bodies::Vector{Body}, c_states::Vector{CanonicalState})
    out = similar(c_states)

    q_dem_i = Vector{SVector{3, Float64}}(undef, length(c_states))
    p_dem_i = Vector{SVector{3, Float64}}(undef, length(c_states))

    masses = mass.(bodies)
    total_mass = sum(masses)
    q_dem_i[1] = SVector{3}(0.0, 0.0, 0.0)
    for ind in eachindex(c_states)
        q_dem_i[1] += masses[ind] * coordinates(c_states[ind])
    end

    q_dem_i[1] /= total_mass
    p_dem_i[1] = sum(momenta, c_states)

    out[1] = CanonicalState(q_dem_i[1][1],
                            q_dem_i[1][2],
                            q_dem_i[1][3],
                            p_dem_i[1][1],
                            p_dem_i[1][2],
                            p_dem_i[1][3])

    for ind in 2:length(c_states)
        q_dem_i[ind] = coordinates(c_states[ind]) - coordinates(c_states[1])
        p_dem_i[ind] = momenta(c_states[ind]) - mass(bodies[ind]) / total_mass * p_dem_i[1]

        out[ind] = CanonicalState(q_dem_i[ind][1],
                                  q_dem_i[ind][2],
                                  q_dem_i[ind][3],
                                  p_dem_i[ind][1],
                                  p_dem_i[ind][2],
                                  p_dem_i[ind][3])
    end

    return out
end

function from_democratic_heliocentric(bodies::Vector{Body}, b_states::Vector{CanonicalState})
    m_total = total_mass(bodies)
    
    q_out = fill(SVector{3}(0.0, 0.0, 0.0), length(b_states))
    p_out = fill(SVector{3}(0.0, 0.0, 0.0), length(b_states))

    masses = mass.(bodies)

    coords = coordinates.(b_states)
    moms = momenta.(b_states)

    mq = sum(gravitational_parameter.(bodies) .* coords) / GRAVITATIONAL_CONSTANT

    out = Vector{CanonicalState}(undef, length(b_states))

    q_out[1] = coords[1] - mq / m_total

    for k in 2:length(bodies)
        q_out[k] = coords[k] + q_out[1] - mq / m_total
        p_out[k] = moms[k] + masses[k] / m_total * moms[1]

        out[k] = CanonicalState(q_out[k]..., p_out[k]...)
    end

    p_out[1] = moms[1] - sum(p_out)

    out[1] = CanonicalState(q_out[1]..., p_out[1]...)

    return out
end

struct SimulationResult
    start_et::Float64
    stop_et::Float64
    time_step::Float64
    bodies::Vector{Body}
    time_history::Vector{Float64}
    state_history::Matrix{State}
    canonical_history::Matrix{CanonicalState}
    barycentric_history::Matrix{CanonicalState}
end

total_mass(bs::Vector{Body}) = sum(gravitational_parameter, bs) / GRAVITATIONAL_CONSTANT
total_mass(s::SimulationResult) = total_mass(s.bodies)
time_step(s::SimulationResult) = s.time_step
step_count(s::SimulationResult) = length(s.time_history)
body_count(s::SimulationResult) = length(s.bodies)

function simulation_result(start::String, stop::String, dt, bodies)
    n_bodies = length(bodies)

    et0 = ephemeris_time(start)
    etf = ephemeris_time(stop)

    time_history = collect(et0:dt:etf)
    n_times = length(time_history)

    time_history[1] = et0

    s_hist = fill(State(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), n_times, n_bodies)
    c_hist = fill(CanonicalState(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), n_times, n_bodies)
    b_hist = fill(CanonicalState(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), n_times, n_bodies)

    init_states = State.(bodies, Ref(et0))

    s_hist[1, :] = init_states
    c_hist[1, :] = CanonicalState.(bodies, init_states)
    b_hist[1, :] = to_democratic_heliocentric(bodies, c_hist[1, :])

    return SimulationResult(et0,
                            etf,
                            dt,
                            bodies,
                            time_history,
                            s_hist,
                            c_hist,
                            b_hist)
end

# Update coordinates with result of linear drift
function _linear_drift!(sim::SimulationResult, dt, a, b)
    m0 = mass(sim.bodies[1])

    for k in 2:length(sim.bodies)
        cs = sim.barycentric_history[a, k]
        mom = momenta(cs)
        sim.barycentric_history[b, k] = sim.barycentric_history[a, k] +
                                        CanonicalState((mom / m0 * dt)..., 0.0, 0.0,
                                                       0.0)
    end
end

function _interaction_kick!(sim::SimulationResult, dt, a, b)
    for i in 2:(body_count(sim) - 1)
        xi = sim.barycentric_history[a, i]
        ci = coordinates(xi)
        gmi = gravitational_parameter(sim.bodies[i])
        dp = SVector{3}(0.0, 0.0, 0.0)
        for j in (i + 1):body_count(sim)
            xj = sim.barycentric_history[a, j]
            cj = coordinates(xj)
            mj = mass(sim.bodies[j])
            cdiff = ci - cj

            cdiffmag = norm(cdiff)

            dp += -gmi * mj / cdiffmag^3 * cdiff
        end
        sim.barycentric_history[a, i] += dt * CanonicalState(0.0, 0.0, 0.0,
                                                        dp[1], dp[2], dp[3])
    end
end

function _kepler_drift!(sim::SimulationResult, dt, a, b)
    gm = gravitational_parameter(sim.bodies[1])
    for k in 2:length(sim.bodies)
        s = State(sim.bodies[k], sim.barycentric_history[a, k])
        new_state = kepler_propagate(gm, dt, s; tol=1e-5)
        c = CanonicalState(sim.bodies[k], new_state)
        sim.barycentric_history[b, k] = c
    end
end

function _set_state_canon!(sim::SimulationResult, step_num)
    sim.canonical_history[step_num, :] .= from_democratic_heliocentric(sim.bodies, sim.barycentric_history[step_num, :])
    sim.state_history[step_num, :] .= State.(sim.bodies, sim.canonical_history[step_num, :])
end

function simulate(start::String, stop::String, dt, bodies)
    sim = simulation_result(start, stop, dt, bodies)

    for step_num in 2:step_count(sim)
        @info "Running step $step_num"
        tstep = time_step(sim)
        _linear_drift!(sim, tstep / 2, step_num-1, step_num)
        _interaction_kick!(sim, tstep / 2, step_num, step_num)
        _kepler_drift!(sim, tstep, step_num, step_num)
        _linear_drift!(sim, tstep / 2, step_num, step_num)
        _interaction_kick!(sim, tstep / 2, step_num, step_num)
        _set_state_canon!(sim, step_num)
    end

    return sim
end

end