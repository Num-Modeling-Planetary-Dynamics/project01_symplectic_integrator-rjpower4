using Symple
using GLMakie


function main()
    central_body = Symple.EARTH
    elements_0 = Symple.KeplerianElements(1e6; ecc=0.8)
    @show x0 = Symple.LagrangianState(central_body.gm, elements_0)
    period = Symple.period(central_body.gm, elements_0)
    n_steps = 1000

    times = range(0.0, period; length=n_steps)

    states = Symple.kepler_propagate.(Ref(central_body.gm), times, Ref(x0))


    x_values = [s[1] for s in states]
    y_values = [s[2] for s in states]
    
    vx_values = [s[4] for s in states]
    vy_values = [s[5] for s in states]

    f = Figure(fontname="IBM Plex Sans", fontsize=20)
    ax = Axis(f[1,1], aspect=DataAspect(), xlabel="X-Position", ylabel="Y-Position")
    lines!(ax, x_values, y_values, linewidth=3)
    ax = Axis(f[1,2], aspect=DataAspect(), xlabel="X-Velocity", ylabel="Y-Velocity")
    lines!(ax, vx_values, vy_values, linewidth=3)
    f
end
