struct State <: FieldVector{6, Float64}
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
end

Base.position(s::State) = SVector{3}(s[1], s[2], s[3])
velocity(s::State) = SVector{3}(s[4], s[5], s[6])
