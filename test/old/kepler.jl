using Symple: kepler_danby, LagrangianState, mean_motion, kepler_propagate, KeplerianElements, period
using Test

@test all(x -> first(kepler_danby(x, π)) ≈ π, 0.1:0.05:0.9)
@test all(x -> first(kepler_danby(x, 0)) ≈ 0, 0.1:0.05:0.9)

let (result, did_converge) = kepler_danby(0.1, 0.3)
    @test result ≈ 3.32655400424575864449e-01
    @test did_converge
end

let (result, did_converge) = kepler_danby(0.9, 0.3)
    @test result ≈ 1.10351772030308703521
    @test did_converge
end

let gm = 398600.4354360959
    r = 10_000.0
    v = 4.9332898193068
    n = mean_motion(gm, r)
    p = 2π / n
    x0 = LagrangianState(r, 0.0, 0.0, 0.0, v, 0.0)
    els = KeplerianElements(gm, x0)
    final = kepler_propagate(gm, period(gm, els), x0)
    @test x0 ≈ final
end
