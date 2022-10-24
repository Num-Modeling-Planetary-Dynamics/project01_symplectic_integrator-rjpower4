using Symple: is_kernel_path
using Test

@testset "Is Kernel Path" begin
    @test is_kernel_path("de430.bsp"; require_exists = false)
    @test is_kernel_path("naif0012.tls"; require_exists = false)
    @test is_kernel_path("gm_de431.tpc"; require_exists = false)
    @test !is_kernel_path("de430.txt"; require_exists = false)
    @test !is_kernel_path("de430.bsp"; require_exists = true)
end
