# ========================================================================================
# File   : src/Symple.jl
# Brief  : Top-level Symple module definition
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================
module Symple

using LinearAlgebra
using StaticArrays: FieldVector, @SVector, SVector, SMatrix

include("path.jl")
include("constants.jl")
include("body.jl")
include("units.jl")
include("state.jl")
include("kepler.jl")
include("elements.jl")
include("simulation.jl")

end
