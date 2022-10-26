# ========================================================================================
# File   : pluto_orbit_demo.jl
# Brief  : Demo visualization comparing the orbit of Pluto with a Kepler-propagated one
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================
using Symple: SUN, PLUTO_BARYCENTER, EARTH_BARYCENTER
using Symple: ephemeris_time, SECONDS_PER_DAY, LagrangianState, spice_state
using Symple: kepler_propagate, kepler_propagate2
using GLMakie

function main(start_epoch, stop_epoch, time_step_days)
    primary = SUN
    secondary = EARTH_BARYCENTER

    et0 = ephemeris_time(start_epoch);
    etf = ephemeris_time(stop_epoch);

    dt_sec = time_step_days * SECONDS_PER_DAY

    times = et0:dt_sec:etf

    x0_sun = LagrangianState(spice_state(primary.naif_id, et0))
    x0_pluto = LagrangianState(spice_state(secondary.naif_id, et0))

    x_pluto_ephem = spice_state.(Ref(secondary.naif_id), times)

    x_pluto_kep = LagrangianState[]

    for time in times
        out = kepler_propagate(primary.gm, time - times[1], x0_pluto)
        push!(x_pluto_kep, out)
    end

    f = Figure()
    ax = Axis(f[1, 1])
    lines!(ax, times, first.(x_pluto_ephem))
    lines!(ax, times, first.(x_pluto_kep))

    f
end
