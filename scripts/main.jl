# ========================================================================================
# File   : scripts/main.jl
# Brief  : Main solution script
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================
using Symple
using SPICE: str2et, spkez, furnsh, kclear
using GLMakie
using LinearAlgebra

# ----------------------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------------------
const CENTRAL_BODY = Symple.SUN
const SECONDARY_BODIES = [Symple.NEPTUNE_BARYCENTER, Symple.PLUTO_BARYCENTER]

# ----------------------------------------------------------------------------------------
# Spice Utilities
# ----------------------------------------------------------------------------------------
function load_kernels()
    kernels = readdir(Symple.data_dir("kernels"); join=true)
    for kernel in kernels
        furnsh(kernel)
    end
end

function clear_kernels()
    return kclear()
end

# ----------------------------------------------------------------------------------------
# Get the Heliocentric states at an epoch
# ----------------------------------------------------------------------------------------
function helio_states_at_epoch(epoch::String; frame="ECLIPJ2000")
    et = str2et(epoch)
    id0 = Symple.naif_id(CENTRAL_BODY)
    idi = Symple.naif_id.(SECONDARY_BODIES)

    x0 = Symple.State(spkez(id0, et, frame, "NONE", id0)[1])

    xi = Vector{Symple.State}(undef, length(idi))

    for (i, id) in enumerate(idi)
        xi[i] = Symple.State(spkez(id, et, frame, "NONE", id0)[1])
    end

    return et, x0, xi
end

# ----------------------------------------------------------------------------------------
# Energy Computation
# ----------------------------------------------------------------------------------------
function energy(s::Symple.Simulation, x0_helios, xi_helios, x0_barys, xi_barys)
    g_const = Symple.gravitational_constant(s)
    gm0, gmi = Symple.gravity_parameters(s)

    x0_bary_vel_mag = norm.(Symple.velocities.(x0_barys))
    bary_vel_mags = norm.(Symple.velocities.(xi_barys))

    # Kinetic
    energy = (1//2) .* (gm0 .* x0_bary_vel_mag)
    for (i, k) in enumerate(eachrow(bary_vel_mags))
        energy[i] += (1//2) * sum(gmi .* k .^ 2)
    end
    energy ./= g_const

    # Potential
    for (row_num, r) in enumerate(eachrow(xi_barys))
        for i in eachindex(r)
            gmi = Symple.gravity_parameter(s, i)
            for j in eachindex(r)
                gmj = Symple.gravity_parameter(s, j)
                if i == j
                    continue
                end

                dr = Symple.coordinates(r[j]) - Symple.coordinates(r[i])
                energy[row_num] -= gmi * gmj / norm(dr) / g_const
            end

            energy[row_num] -= gm0 * gmi / norm(Symple.coordinates(r[i])) / g_const
        end
    end

    return energy
end

# ----------------------------------------------------------------------------------------
# Energy Computation
# ----------------------------------------------------------------------------------------
function resonance_angle(s::Symple.Simulation, x0_helios, xi_helios, x0_barys, xi_barys)
    g_const = Symple.gravitational_constant(s)
    gm0, gmi = Symple.gravity_parameters(s)

    elements = Symple.KeplerianElements.(Ref(gm0), xi_helios)
    lops = Symple.longitude_of_periapsis.(elements)
    mean_anoms = Symple.mean_anomaly.(elements)
    mlons = mean_anoms .+ lops

    return 3 .* mlons[:, 2] - 2 .* mlons[:, 1] - lops[:, 2]
end

# ----------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------
function main()
    load_kernels()
    # Construct Scaling
    units = Symple.UnitSet(;
        mass=Symple.gravity_parameter(CENTRAL_BODY) / Symple.GRAVITATIONAL_CONSTANT,
        length=Symple.ASTRONOMICAL_UNIT,
        time=Symple.SECONDS_PER_YEAR,
    )

    # Get initial states
    et0, x0_helio, xi_helio = helio_states_at_epoch("JAN 01 2022")
    etf = et0 + 100_000 * Symple.SECONDS_PER_YEAR
    dt = 1 * Symple.SECONDS_PER_YEAR
    nsteps = Int32(round((etf - et0) / dt))

    # Build Simulation
    s = Symple.Simulation(CENTRAL_BODY, SECONDARY_BODIES; units=units)

    clear_kernels()

    x0_helios, xi_helios, x0_barys, xi_barys = Symple.simulate(
        s, dt, x0_helio, xi_helio; n_steps=nsteps
    )

    en = energy(s, x0_helios, xi_helios, x0_barys, xi_barys)
    ra = resonance_angle(s, x0_helios, xi_helios, x0_barys, xi_barys)
    ts = (dt .* (0:(length(en) - 1))) ./ Symple.SECONDS_PER_YEAR

    f = Figure(; font="Roboto", fontsize=25, resolution=(1500, 1000))
    # ax = Axis3(f[1:2, 1]; aspect=:data, title="Heliocentric Frame")

    # x1 = [a[1] for a in xi_helios[:, 1]]
    # y1 = [a[2] for a in xi_helios[:, 1]]
    # z1 = [a[3] for a in xi_helios[:, 1]]

    # x2 = [a[1] for a in xi_helios[:, 2]]
    # y2 = [a[2] for a in xi_helios[:, 2]]
    # z2 = [a[3] for a in xi_helios[:, 2]]

    # lines!(ax, x1, y1, z1)
    # lines!(ax, x2, y2, z2)

    ax = Axis(f[1, 1], xlabel="Time [year]", ylabel=L"$\frac{\mathcal{E} - \mathcal{E}_0}{\mathcal{E}_0}$")
    lines!(ax, ts, (en .- en[1]) ./ en[1])

    ax = Axis(f[2, 1], xlabel="Time [year]", ylabel="Resonance Angle")
    lines!(ax, ts, rad2deg.(mod.(ra, 2π)))
    return f
end