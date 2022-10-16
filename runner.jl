using Symple
using DataFrames

function run()
    n_years = 100_000

    start_epoch = "OCT 17 2022"

    start_et = Symple.ephemeris_time(start_epoch)
    final_et = Symple.SECONDS_PER_YEAR * n_years + start_et

    bodies = [Symple.SUN, Symple.NEPTUNE, Symple.PLUTO]
    init_lag_states = Symple.LagrangianState.(bodies, Ref(start_et))

    ke_neptune = Symple.KeplerianElements(bodies[1].gm, init_lag_states[2])
    period_neptune = Symple.period(bodies[1].gm, ke_neptune)
    dt = period_neptune / 30

    simulation = Symple.Simulation((start_et, final_et), bodies; max_step=dt)

    dem_states = Symple.democratic_heliocentric(bodies, init_lag_states)

    Symple.init!(simulation, dem_states)

    while Symple.do_step!(simulation) > 0
    end

    ham_states = similar(simulation.states)
    lag_states = Matrix{Symple.LagrangianState}(undef, size(ham_states)[1], length(bodies))
    for k in 1:size(ham_states)[1]
        ham_states[k, :] = Symple.from_democratic_heliocentric(bodies, simulation.states[k, :])

        lag_states[k, :] .= Symple.LagrangianState.(bodies, ham_states[k, :])
    end

    t = range(start_et, final_et; length=size(ham_states)[1])
    energy = Symple.system_energy(simulation)


    out = DataFrame(
        :t => t,
        :x_sun => lag_states[:,1],
        :x_neptune => lag_states[:, 2],
        :x_pluto => lag_states[:, 3],
        :energy => energy,
    )

    # t, Symple.specific_energy.(Ref(Symple.gravitational_parameter(bodies[1])), lag_states[:, 3]), simulation
    return  out
end
