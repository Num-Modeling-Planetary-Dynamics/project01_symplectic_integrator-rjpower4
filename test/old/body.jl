using Symple: Body, name, gravity_parameter, naif_id
using Test

@test_throws ArgumentError Body("", 199, 1.1)
@test_throws ArgumentError Body("", 199, 0.0)
@test_throws ArgumentError Body("body", 199, 0.0)
@test_throws ArgumentError Body("body", 199, -10.0)

let b = Body(gm = 1.1, name = "name", naif_id = 12)
    @test name(b) == "name"
    @test gravity_parameter(b) == 1.1
    @test naif_id(b) == 12
end
