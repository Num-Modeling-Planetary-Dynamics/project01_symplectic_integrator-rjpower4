using Symple: KeplerianElements, State
using Symple: semi_major_axis, eccentricity, inclination, right_ascension
using Symple: argument_of_periapsis, true_anomaly

@testset "Elements 1" begin
    gm = 1.3271927514442485E+11
    EC = 1.414884439617544E-02
    QR = 4.469325595222297E+09
    IN = 1.763253042672047E+00
    OM = 1.315860830709717E+02
    W = 2.473162281457197E+02
    Tp = 2463985.471495890059
    N = 6.838229398826735E-08
    MA = 3.339744143752608E+02
    TA = 3.332515540507960E+02
    A = 4.533468941855504E+09
    AD = 4.597612288488711E+09
    PR = 5.264520667612683E+09

    X = 4.433073704290709E+09
    Y = -6.119493983258691E+08
    Z = -8.956937853178886E+07
    VX = 7.124342362655471E-01
    VY = 5.431693496590404E+00
    VZ = -1.273893581279648E-01
    LT = 1.493035543390672E+04
    RG = 4.476007954344552E+09
    RR = -3.445880888460209E-02

    x0 = State(X, Y, Z, VX, VY, VZ)
    ke = KeplerianElements(gm, x0)

    @test semi_major_axis(ke) ≈ A
    @test eccentricity(ke) ≈ EC
    @test inclination(ke) ≈ deg2rad(IN)
    @test argument_of_periapsis(ke) ≈ deg2rad(W)
    @test right_ascension(ke) ≈ deg2rad(OM)
    @test true_anomaly(ke) ≈ deg2rad(TA)

    x0_out = State(gm, ke)

    @test x0_out[1] ≈ X
    @test x0_out[2] ≈ Y
    @test x0_out[3] ≈ Z
    @test x0_out[4] ≈ VX
    @test x0_out[5] ≈ VY
    @test x0_out[6] ≈ VZ
end

@testset "Elements 2" begin
    gm = 1.3271927514442485E+11
    EC = 1.446228889250229E-02
    QR = 4.468165734677928E+09
    IN = 1.768588425072905E+00
    OM = 1.317337657552168E+02
    W = 2.491768455221591E+02
    Tp = 2464314.475393503904
    N = 6.837629943742825E-08
    MA = 3.330373465614567E+02
    TA = 3.322736558090996E+02
    A = 4.533733904161646E+09
    AD = 4.599302073645364E+09
    PR = 5.264982208190999E+09

    X = 4.442809943627258E+09
    Y = -5.320757054778891E+08
    Z = -9.143575357269061E+07
    VX = 6.127870378712618E-01
    VY = 5.444652951534985E+00
    VZ = -1.260306289806319E-01
    LT = 1.492863313086395E+04
    RG = 4.475491620881938E+09
    RR = -3.640885144227313E-02

    x0 = State(X, Y, Z, VX, VY, VZ)
    ke = KeplerianElements(gm, x0)

    @test semi_major_axis(ke) ≈ A
    @test eccentricity(ke) ≈ EC
    @test inclination(ke) ≈ deg2rad(IN)
    @test argument_of_periapsis(ke) ≈ deg2rad(W)
    @test right_ascension(ke) ≈ deg2rad(OM)
    @test true_anomaly(ke) ≈ deg2rad(TA)

    x0_out = State(gm, ke)

    @test x0_out[1] ≈ X
    @test x0_out[2] ≈ Y
    @test x0_out[3] ≈ Z
    @test x0_out[4] ≈ VX
    @test x0_out[5] ≈ VY
    @test x0_out[6] ≈ VZ
end
