using Symple
using GLMakie

function step_size(primary, secondary, x_secondary, n_per_period)
    gm_0 = Symple.gravity_parameter(primary)
    kes = Symple.KeplerianElements.(Ref(gm_0), x_secondary)
    periods = Symple.period.(Ref(gm_0), kes)
    p_min = minimum(periods)
    dt = p_min / n_per_period
    return dt
end

function main(n_years)
    epoch_init = "JAN 01 2023"
    et_init = Symple.ephemeris_time(epoch_init)

    primary = SUN
    secondaries = [NEPTUNE_BARYCENTER, PLUTO_BARYCENTER]

    x_h_0 = LagrangianState(spice_state(naif_id(primary), et_init; frame="ECLIPJ2000"))
    x_h_i = LagrangianState.(spice_state.(naif_id.(secondaries), et_init; frame="ECLIPJ2000"))

    s = Simulation(primary, secondaries)
    ss_h = SimulationState(x_h_0, x_h_i)
    ss = Symple.heliocentric_to_barycentric(s, ss_h)

    dt = step_size(primary, secondaries, x_h_i, 40)

    times = et_init:dt:(et_init + n_years * Symple.SECONDS_PER_YEAR)

    x1 = Float64[]
    y1 = Float64[]
    z1 = Float64[]

    x2 = Float64[]
    y2 = Float64[]
    z2 = Float64[]

    energies = Float64[]

    for t in times
        push!(x1, ss.secondary_states[1][1])# - ss.primary_state[1])
        push!(y1, ss.secondary_states[1][2])# - ss.primary_state[2])
        push!(z1, ss.secondary_states[1][3])# - ss.primary_state[3])

        push!(x2, ss.secondary_states[2][1])# - ss.primary_state[1])
        push!(y2, ss.secondary_states[2][2])# - ss.primary_state[2])
        push!(z2, ss.secondary_states[2][3])# - ss.primary_state[3])

        push!(energies, Symple.energy(ss, s))

        Symple.do_step!(ss, s, dt)
    end

    f = Figure();
    ax = Axis3(f[1,1])
    lines!(ax, x1, y1, z1)
    lines!(ax, x2, y2, z2)
    scale!(ax.scene, 1.0, 1.0, 1.0)


    ax2 = Axis(f[1, 2])
    lines!(ax2, (energies .- energies[1]) / energies[1])

    f
end