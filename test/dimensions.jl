using Symple: DimensionSet, InvalidDimensionError, mass_unit, length_unit, time_unit
using Symple: velocity_unit, acceleration_unit, force_unit
using Test

# Valid Constructions
@test DimensionSet(2.0) == DimensionSet(2.0, 2.0, 2.0)
@test DimensionSet(mass = 2.0, length = 2.0, time = 2.0) == DimensionSet(2.0, 2.0, 2.0)

# Default Construction
@test DimensionSet() == DimensionSet(mass = 1.0, length = 1.0, time = 1.0)
@test DimensionSet() == DimensionSet(1.0)

# Invalid Constructions
@test_throws InvalidDimensionError DimensionSet(mass = 0.0)
@test_throws InvalidDimensionError DimensionSet(length = 0.0)
@test_throws InvalidDimensionError DimensionSet(time = 0.0)
@test_throws InvalidDimensionError DimensionSet(mass = -10.0)
@test_throws InvalidDimensionError DimensionSet(length = -10.0)
@test_throws InvalidDimensionError DimensionSet(time = -10.0)
@test_throws InvalidDimensionError DimensionSet(0.0)
@test_throws InvalidDimensionError DimensionSet(-10.0)

# Assigning values
let ds = DimensionSet(mass = 4.1, length = 5.2, time = 6.3)
    @test mass_unit(ds) == 4.1
    @test length_unit(ds) == 5.2
    @test time_unit(ds) == 6.3
end

# Velocity
@test velocity_unit(DimensionSet()) == 1.0
@test velocity_unit(DimensionSet(3.0)) == 1.0
@test velocity_unit(DimensionSet(length = 4.0, time = 2.0)) ≈ 2.0
@test velocity_unit(DimensionSet(length = 3.0, time = 6.0)) ≈ 0.5

# Acceleration
@test acceleration_unit(DimensionSet()) == 1.0
@test acceleration_unit(DimensionSet(4.0)) ≈ 1 // 4
@test acceleration_unit(DimensionSet(length = 9.0, time = 3.0)) ≈ 1.0

# Force
@test force_unit(DimensionSet()) == 1.0
@test force_unit(DimensionSet(5.0)) == 1.0
@test force_unit(DimensionSet(mass = 2.0, length = 3.0, time = 4.0)) ≈ 0.375
