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

function main(n_years, n=40)
    epoch_init = "JAN 01 2023"
    et_init = Symple.ephemeris_time(epoch_init)

    primary = SUN
    secondaries = [NEPTUNE_BARYCENTER, PLUTO_BARYCENTER]

    x_h_0 = LagrangianState(spice_state(naif_id(primary), et_init; frame="ECLIPJ2000"))
    x_h_i = LagrangianState.(spice_state.(naif_id.(secondaries), et_init; frame="ECLIPJ2000"))

    s = Simulation(primary, secondaries)
    ss_h = SimulationState(x_h_0, x_h_i)
    ss = Symple.heliocentric_to_barycentric(s, ss_h)

    dt = step_size(primary, secondaries, x_h_i, n)

    times = et_init:dt:(et_init + n_years * Symple.SECONDS_PER_YEAR)

    x1 = Float64[]
    y1 = Float64[]
    z1 = Float64[]

    x2 = Float64[]
    y2 = Float64[]
    z2 = Float64[]

    energies = Float64[]
    kes1 = Symple.KeplerianElements[]
    kes2 = Symple.KeplerianElements[]

    for t in times
        push!(x1, ss.secondary_states[1][1] - ss.primary_state[1])
        push!(y1, ss.secondary_states[1][2] - ss.primary_state[2])
        push!(z1, ss.secondary_states[1][3] - ss.primary_state[3])

        push!(x2, ss.secondary_states[2][1] - ss.primary_state[1])
        push!(y2, ss.secondary_states[2][2] - ss.primary_state[2])
        push!(z2, ss.secondary_states[2][3] - ss.primary_state[3])

        push!(energies, Symple.energy(ss, s))

        push!(kes1, Symple.KeplerianElements(primary.gm, ss.secondary_states[1] - ss.primary_state))
        push!(kes2, Symple.KeplerianElements(primary.gm, ss.secondary_states[2] - ss.primary_state))

        Symple.do_step!(ss, s, dt)
    end

    ma1 = Symple.mean_anomaly.(Symple.eccentricity.(kes1), Symple.eccentric_anomaly.(kes1))
    ma2 = Symple.mean_anomaly.(Symple.eccentricity.(kes2), Symple.eccentric_anomaly.(kes2))

    varpi1 = Symple.right_ascension.(kes1) + Symple.argument_of_periapsis.(kes1)
    varpi2 = Symple.right_ascension.(kes2) + Symple.argument_of_periapsis.(kes2)

    lambda1 = ma1 + varpi1
    lambda2 = ma2 + varpi2

    f = Figure();
    ax = Axis3(f[1:2,1])
    lines!(ax, x1, y1, z1)
    lines!(ax, x2, y2, z2)
    scale!(ax.scene, 1.0, 1.0, 1.0)


    ax2 = Axis(f[1, 2])
    lines!(ax2, times[2:end] / Symple.SECONDS_PER_YEAR, (energies[2:end] .- energies[1:end-1]) ./ energies[1:end-1])

    ax3 = Axis(f[2, 2])
    lines!(ax3, times / Symple.SECONDS_PER_YEAR , 3lambda2 - 2lambda1 - varpi2)
    f
end