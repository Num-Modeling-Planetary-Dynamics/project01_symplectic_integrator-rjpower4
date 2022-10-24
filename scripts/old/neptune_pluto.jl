# ========================================================================================
# File: neptune_pluto.jl
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

using Symple
using Symple: LagrangianState, HamiltonianState, Simulation, democratic_heliocentric
using Symple: SUN, NEPTUNE, PLUTO, gravity_parameter
using Symple: spice_state, naif_id, period

function spice_init()
    kernels = Symple.list_kernels("data/kernels")
    Symple.load_kernels(kernels)
end

function determine_step(start, stop, max_step)
    return (stop - start) / ceil((stop - start) / max_step)
end

function main(n_years; save_every=1000)
    # Timing
    start_epoch = "JAN 01 2023"
    start_et = Symple.ephemeris_time(start_epoch)
    stop_et = n_years * Symple.SECONDS_PER_YEAR + start_et

    # Initial Lagrangian States
    x0_sun = LagrangianState(spice_state(naif_id(SUN), start_et; observer=10))
    x0_nep = LagrangianState(spice_state(naif_id(NEPTUNE), start_et; observer=10))
    x0_plu = LagrangianState(spice_state(naif_id(PLUTO), start_et; observer=10))

    # Hamiltonian
    p0_sun = HamiltonianState(gravity_parameter(SUN) / Symple.GRAVITATIONAL_CONSTANT, x0_sun)
    p0_nep = HamiltonianState(gravity_parameter(NEPTUNE) / Symple.GRAVITATIONAL_CONSTANT, x0_nep)
    p0_plu = HamiltonianState(gravity_parameter(PLUTO) / Symple.GRAVITATIONAL_CONSTANT, x0_plu)

    # DHC
    (s0, sn) = democratic_heliocentric(SUN, [NEPTUNE, PLUTO], p0_sun, [p0_nep, p0_plu])

    # Determine the step size
    max_step = period(SUN, Symple.KeplerianElements(SUN, x0_nep)) / 30
    dt = determine_step(start_et, stop_et, max_step)

    # Build time span
    tspan = start_et:dt:stop_et

    # Build a dimension set corresponding to AU/YEAR/KG
    dimensions = Symple.DimensionSet(mass=1.0, length=Symple.ASTRONOMICAL_UNIT, time=Symple.SECONDS_PER_YEAR)

    # Build the simulation
    sim = Simulation(dt, SUN, [NEPTUNE, PLUTO]; dimensions=dimensions)
    Symple.init!(sim, start_et, s0, sn)

    e0 = Symple.energy(sim)
    eccs = Float64[]

    last_e = e0
    times = start_et:dt:stop_et
    des = Float64[]
    count = 1;

    for (i, t) in enumerate(times)
        Symple.do_step!(sim)
        enew = Symple.energy(sim)
        
        if mod(i, save_every) == 0
            push!(des, (enew - last_e) / last_e)
            count += 1
        end
        last_e = enew
    end
    return des
end
