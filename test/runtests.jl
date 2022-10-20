using SafeTestsets
using Test

@safetestset "Utility Tests" begin include("util.jl") end
@safetestset "Dimension Tests" begin include("dimensions.jl") end
@safetestset "Body Tests" begin include("body.jl") end
@safetestset "Spice Tests" begin include("spice.jl") end
@safetestset "Elements Tests" begin include("elements.jl") end
@safetestset "Kepler Tests" begin include("kepler.jl") end
