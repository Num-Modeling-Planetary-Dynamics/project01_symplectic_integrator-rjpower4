# ========================================================================================
# File   : test/runtests.jl
# Brief  : Top-level test running entry point for Symple
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================
using SafeTestsets: @safetestset

@safetestset "Path utility tests" begin
    include("path.jl")
end

@safetestset "Constants tests" begin
    include("constants.jl")
end

@safetestset "Units tests" begin
    include("units.jl")
end

@safetestset "Kepler tests" begin
    include("kepler.jl")
end

@safetestset "State tests" begin
    include("state.jl")
end

@safetestset "Elements tests" begin
    include("elements.jl")
end

@safetestset "Simulation tests" begin
    include("simulation.jl")
end
