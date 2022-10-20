using Symple: kepler_danby
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


