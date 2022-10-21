# ========================================================================================
# File: Symple.jl
# Author: Rolfe Power
# ========================================================================================
module Symple

# ----------------------------------------------------------------------------------------
# External packages
# ----------------------------------------------------------------------------------------
using LinearAlgebra: norm, cross, dot
using StaticArrays: FieldVector, SVector, SMatrix, MVector
using SPICE: furnsh, spkez, spkezr, str2et

# ----------------------------------------------------------------------------------------
# Internal Includes
# ----------------------------------------------------------------------------------------
include("util.jl")
include("constants.jl")
include("dimensions.jl")
include("body.jl")
include("spice.jl")
include("state.jl")
include("twobody.jl")
include("elements.jl")
include("kepler.jl")
include("simulation.jl")
include("dhc.jl")

end
