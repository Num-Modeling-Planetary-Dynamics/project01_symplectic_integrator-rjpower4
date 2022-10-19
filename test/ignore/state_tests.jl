using Test
using EAPS
using StaticArrays: @SVector

let s = EAPS.State(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    @test s[1] == 1.0
    @test s[2] == 2.0
    @test s[3] == 3.0
    @test s[4] == 4.0
    @test s[5] == 5.0
    @test s[6] == 6.0

    r = @SVector [1.0, 2.0, 3.0]
    v = @SVector [4.0, 5.0, 6.0]

    @test EAPS.position(s) == r
    @test EAPS.velocity(s) == v

    @test (r, v) == EAPS.posvel(s)
end
