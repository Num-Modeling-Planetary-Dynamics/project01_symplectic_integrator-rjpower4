# ========================================================================================
# File: neptune_pluto.jl
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

using Symple
using Symple: LagrangianState, HamiltonianState, Simulation, NBodySystem
using Symple: SUN, NEPTUNE, PLUTO
using Symple: spice_state, naif_id, period

function spice_init()
    kernels = Symple.list_kernels("data/kernels")
    Symple.load_kernels(kernels)
end

function determine_step(start, stop, max_step)
    return ceil((stop - start) / max_step)
end

function main()
    # Timing
    start_epoch = "JAN 01 2023"
    start_et = Symple.ephemeris_time(start_epoch)
    stop_et = 10_000 * Symple.SECONDS_PER_YEAR + start_et

    # Initial Lagrangian States
    x0_sun = LagrangianState(spice_state(naif_id(SUN), start_et))
    x0_nep = LagrangianState(spice_state(naif_id(NEPTUNE), start_et))
    x0_plu = LagrangianState(spice_state(naif_id(PLUTO), start_et))

    # Determine the step size
    max_step = period(SUN, Symple.KeplerianElements(SUN, x0_nep)) / 30
    dt = determine_step(start_et, stop_et, max_step)

    # Build time span
    tspan = start_et:dt:stop_et

    # Build a dimension set corresponding to AU/YEAR/KG
    dimensions = Symple.DimensionSet(mass=1.0, length=Symple.ASTRONOMICAL_UNIT, time=Symple.SECONDS_PER_YEAR)

    # Build Simulation
    sim = Simulation(NBodySystem(SUN, [NEPTUNE, PLUTO]), tspan; dimensions=dimensions)
    Symple.init(sim, [x0_sun, x0_nep, x0_plu])
end
