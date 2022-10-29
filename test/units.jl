# ========================================================================================
# File   : test/units.jl
# Brief  : Tests for unit set
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================
using Test
using Symple: UnitSet, mass_unit, length_unit, time_unit, velocity_unit
using Symple: acceleration_unit, force_unit, gm_unit, energy_unit
using Symple

@test UnitSet(1.0, 2.0, 3.0) == UnitSet(1.0, 2.0, 3.0)

@testset "Bad Constructions" begin
    @test_throws ArgumentError UnitSet(0.0, 2.0, 3.0)
    @test_throws ArgumentError UnitSet(1.0, 0.0, 3.0)
    @test_throws ArgumentError UnitSet(1.0, 2.0, 0.0)

    @test_throws ArgumentError UnitSet(-1.0, 2.0, 3.0)
    @test_throws ArgumentError UnitSet(1.0, -2.0, 3.0)
    @test_throws ArgumentError UnitSet(1.0, 2.0, -3.0)
end

@testset "Constant Scaling" begin
    @test UnitSet(1) == UnitSet(1.0, 1.0, 1.0)
    @test UnitSet() == UnitSet(1.0, 1.0, 1.0)
    @test UnitSet(12.3) == UnitSet(12.3, 12.3, 12.3)
end

@testset "Custom Scaling" begin
    @test UnitSet(; mass=3.0, length=5.0, time=2.0) == UnitSet(3.0, 5.0, 2.0)
end

@testset "Basic Accessors" begin
    let unit_set = UnitSet(; mass=3.1, length=4.2, time=5.3)
        @test mass_unit(unit_set) == 3.1
        @test length_unit(unit_set) == 4.2
        @test time_unit(unit_set) == 5.3
    end
end

@testset "Velocity Computation" begin
    @test velocity_unit(UnitSet()) ≈ 1.0
    @test velocity_unit(UnitSet(; length=2.0, time=2.0)) ≈ 1.0
    @test velocity_unit(UnitSet(; length=5.0, time=2.0)) ≈ 2.5

    # The velocity unit should be independent of the mass
    let f = x -> velocity_unit(UnitSet(; mass=x, length=1.2, time=4.3))
        @test f(12.1) == f(5.3)
    end
end

@testset "Acceleration Computation" begin
    @test acceleration_unit(UnitSet()) ≈ 1.0
    @test acceleration_unit(UnitSet(; length=2.0, time=2.0)) ≈ 0.5
    @test acceleration_unit(UnitSet(; length=5.0, time=2.0)) ≈ 1.25

    # The acceleration unit should be independent of the mass
    let f = x -> acceleration_unit(UnitSet(; mass=x, length=1.2, time=4.3))
        @test f(12.1) == f(5.3)
    end
end

@testset "Force Computation" begin
    @test force_unit(UnitSet()) ≈ 1.0
    @test force_unit(UnitSet(; mass=8.0, length=2.0, time=2.0)) ≈ 4.0
    @test force_unit(UnitSet(; mass=0.5, length=5.0, time=2.0)) ≈ 0.625
end

@testset "GM Computation" begin
    @test gm_unit(UnitSet()) ≈ 1.0
    @test gm_unit(UnitSet(; mass=5.0, length=2.0, time=9.3)) ≈ 0.09249624234015491

    # The GM unit should be independent of the mass
    let f = x -> gm_unit(UnitSet(; mass=x, length=1.2, time=4.3))
        @test f(12.1) == f(5.3)
    end
end

@testset "Energy Computation" begin
    @test energy_unit(UnitSet()) ≈ 1.0
    @test energy_unit(UnitSet(; mass=5.0, length=2.0, time=9.3)) ≈ 0.2312406058503873
end

@testset "Unit Set Interface" begin
    struct UnitSetSampleType end
    Symple.unit_set(::UnitSetSampleType) = UnitSet(; mass=5.0, length=2.0, time=9.3)

    @test mass_unit(UnitSetSampleType()) == 5.0
    @test length_unit(UnitSetSampleType()) == 2.0
    @test time_unit(UnitSetSampleType()) == 9.3

    @test velocity_unit(UnitSetSampleType()) ≈ 0.2150537634408602
    @test acceleration_unit(UnitSetSampleType()) ≈ 0.0231240605850387
    @test force_unit(UnitSetSampleType()) ≈ 0.1156203029251937
    @test gm_unit(UnitSetSampleType()) ≈ 0.09249624234015491
    @test energy_unit(UnitSetSampleType()) ≈ 0.2312406058503873
end
