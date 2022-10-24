using Symple: angular_difference
using Test

@test angular_difference(0, 0.1) ≈ 0.1
@test angular_difference(0.1, 0) ≈ 2π - 0.1
@test angular_difference(0, 3π / 4) ≈ 3π / 4
@test angular_difference(3π / 2, π / 2) ≈ π
@test angular_difference(7π / 2, π / 2) ≈ π
@test angular_difference(π / 2, 7π / 2) ≈ π
