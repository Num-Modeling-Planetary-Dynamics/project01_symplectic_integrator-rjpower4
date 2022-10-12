module EAPS

using LinearAlgebra: norm, dot
using StaticArrays: FieldVector, SVector, MVector

include("constants.jl")
include("pathutils.jl")
include("state.jl")
include("twobody.jl")
include("kepler.jl")

end
