using Symple: KeplerianElements
using Symple: semi_major_axis, eccentricity, inclination, argument_of_periapsis
using Symple: right_ascension, true_anomaly
using Test

@test_throws ArgumentError KeplerianElements(0.0)
@test_throws ArgumentError KeplerianElements(-1.0; ecc = 0.5)
@test_throws ArgumentError KeplerianElements(1.0; ecc = 1.5)
@test_throws ArgumentError KeplerianElements(1.0; ecc = 1.0)

let ke = KeplerianElements(1.5; ecc = 0.3, inc = 0.4, aop = 0.5, raan = 0.6, ta = 0.7)
    @test semi_major_axis(ke) == 1.5
    @test eccentricity(ke) == 0.3
    @test inclination(ke) == 0.4
    @test argument_of_periapsis(ke) == 0.5
    @test right_ascension(ke) == 0.6
    @test true_anomaly(ke) == 0.7
end
