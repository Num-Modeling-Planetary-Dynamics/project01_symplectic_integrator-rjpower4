using Symple
using GLMakie

function main(start_et, n_years, steps_per_year)
    TU = Symple.SECONDS_PER_YEAR
    DU = Symple.ASTRONOMICAL_UNIT
    GU = TU^2 / DU^3

    et_init = Symple.ephemeris_time(start_et)
    x0 = Symple.LagrangianState(Symple.spice_state(naif_id(Symple.PLUTO_BARYCENTER), et_init; frame="ECLIPJ2000", observer=10))

    @show dt = Symple.SECONDS_PER_YEAR / steps_per_year
    @show times = et_init:dt:(et_init + n_years*Symple.SECONDS_PER_YEAR)

    states = Vector{Symple.LagrangianState}(undef, length(times))
    elements = Vector{Symple.KeplerianElements}(undef, length(states))
    states[1] = x0;
    elements[1] = Symple.KeplerianElements(Symple.SUN.gm, states[1])

    for i in 2:length(states)
        states[i] = Symple.kepler_propagate2(Symple.SUN.gm, dt, states[i - 1])
        elements[i] = Symple.KeplerianElements(Symple.SUN.gm, states[i])
    end

    tscale = 1 #./ Symple.SECONDS_PER_YEAR

    f2 = Figure(fontname="Roboto");
    ax = Axis3(f2[1:4,1])
    lines!(ax, first.(states), [s[2] for s in states], [s[3] for s in states])
    scale!(ax.scene, 1.0, 1.0, 1.0)


    ax = Axis(f2[5, 1], xlabel="Time [years]", ylabel="Specific Energy")
    lines!(ax, times .* tscale, Symple.specific_energy.(Ref(Symple.SUN), states))

    # f2 = Figure()
    ax = Axis(f2[1, 2], ylabel="Semi-Major Axis")
    lines!(ax, times .* tscale, Symple.semi_major_axis.(elements))
    ax = Axis(f2[2, 2], ylabel="Eccentricity")
    lines!(ax, times .* tscale, Symple.eccentricity.(elements))
    ax = Axis(f2[3, 2], ylabel="Inclination")
    lines!(ax, times .* tscale, Symple.inclination.(elements))
    ax = Axis(f2[4, 2], ylabel="Argument of Periapsis")
    lines!(ax, times .* tscale, Symple.argument_of_periapsis.(elements))
    ax = Axis(f2[5, 2], xlabel="Time [years]", ylabel="Rigth Ascension")
    lines!(ax, times .* tscale, Symple.right_ascension.(elements))

    # di = Symple.inclination.(elements)
    # @show di[2:end] - di[1:end-1]

    f2
end