# ========================================================================================
# File   : test/kepler.jl
# Brief  : Tests for kepler solution
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================
using Symple: kepler_danby

@testset "Solving Kepler's Problem" begin
    @test kepler_danby(0.0, 0.3)[1] ≈ 0.3
    @test kepler_danby(0.0, 2π - 0.01)[1] ≈ (2π - 0.01)

    let res = kepler_danby(0.3, 0.3)[1]
        @test (res - 0.3 * sin(res)) ≈ 0.3
    end

    let ecc = 0.99, res = kepler_danby(ecc, 0.3)[1]
        @test (res - ecc * sin(res)) ≈ 0.3
    end
end