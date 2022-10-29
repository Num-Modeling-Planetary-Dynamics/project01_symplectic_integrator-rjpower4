# ========================================================================================
# File   : test/state.jl
# Brief  : Tests for states
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================
using Test
using StaticArrays: SVector
using Symple: State, coordinates, velocities

@test coordinates(State(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)) == SVector{3}(1.0, 2.0, 3.0)
@test velocities(State(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)) == SVector{3}(4.0, 5.0, 6.0)