using Test
using Symple
using Symple: KeplerianElements, eccentric_anomaly, mean_motion

function random_angle(args...)
    return rand(args...) * 2π
end

@testset "Mean Motion" begin
    @test mean_motion(1, 1) ≈ 1
    @test mean_motion(3, 1) ≈ sqrt(3)
    @test mean_motion(5.0, 2.0) ≈ sqrt(5//8)

    @test all(rand(100) * 10_000) do a
        mean_motion(0.3, a) ≈ mean_motion(0.3, KeplerianElements(a; ecc=0.3))
    end
end

@testset "Eccentric anomaly calculations" begin
    @test eccentric_anomaly(0.4, 0.0) ≈ 0.0
    @test eccentric_anomaly(0.9, 0.0) ≈ 0.0
    @test eccentric_anomaly(0.4, π) ≈ π
    @test eccentric_anomaly(0.9, π) ≈ π

    N = 100

    @test all(random_angle(N)) do value
        eccentric_anomaly(0.0, value) ≈ value
    end

    let ecc = 0.3
        ratio = sqrt((1 + ecc) / (1 - ecc))

        @test all(random_angle(N)) do ta
            ea = eccentric_anomaly(ecc, ta)
            ke = KeplerianElements(rand()*100.0 + 1.0; ecc=ecc, ta=ta)
            ea2 = eccentric_anomaly(ke)
            ((tan(ta / 2) / tan(ea / 2)) ≈ ratio) && (ea2 ≈ ea)
        end
    end

    @test eccentric_anomaly(KeplerianElements(1.1; ecc=0.4, ta=0.0)) ≈ 0.0
    @test eccentric_anomaly(KeplerianElements(1.1; ecc=0.9, ta=0.0)) ≈ 0.0
    @test eccentric_anomaly(KeplerianElements(1.1; ecc=0.4, ta=π)) ≈ π
    @test eccentric_anomaly(KeplerianElements(1.1; ecc=0.9, ta=π)) ≈ π
end

@testset "Keplerian Elements" begin
    let gm = 403503.0
        ke = KeplerianElements(12e3; ecc=0.3, inc=0.4, aop=0.5, raan=0.6, ta=0.7)
        x = LagrangianState(gm, ke)
        @test all(ke .≈ KeplerianElements(gm, x))

        x = LagrangianState(5e3, 6e3, 2e3, 3.5, 2.1, 0.3)
        ke = KeplerianElements(gm, x)
        @test all(x .≈ LagrangianState(gm, ke))
    end
end