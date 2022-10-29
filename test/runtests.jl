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
