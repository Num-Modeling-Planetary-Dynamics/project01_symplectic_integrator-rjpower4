# ========================================================================================
# File   : test/simulation.jl
# Brief  : Tests for the simulation type
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================
using Symple: Simulation, UnitSet, unit_set, velocity_unit, acceleration_unit, force_unit
using Symple: gm_unit, energy_unit, gravity_parameter
using Symple: SUN, EARTH_BARYCENTER, MARS_BARYCENTER

@testset "Unit Set Interface" begin
    us = UnitSet(; mass=5.0, length=2.0, time=9.3)
    s = Simulation(SUN, [EARTH_BARYCENTER, MARS_BARYCENTER]; units=us)

    @test unit_set(s) == us

    for m in (velocity_unit, acceleration_unit, force_unit, gm_unit, energy_unit)
        @test m(s) == m(us)
    end
end

@testset "Gravity Parameter" begin
    us = UnitSet(; mass=5.0, length=2.0, time=9.3)
    s = Simulation(SUN, [EARTH_BARYCENTER, MARS_BARYCENTER])
    @test gravity_parameter(s) == gravity_parameter(SUN)
    @test gravity_parameter(s, 1) == gravity_parameter(EARTH_BARYCENTER)
    @test gravity_parameter(s, 2) == gravity_parameter(MARS_BARYCENTER)
    @test_throws BoundsError gravity_parameter(s, 3)
end
