struct State <: FieldVector{6,Float64}
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
end

const _SZ3 = SVector{3}(0.0, 0.0, 0.0)

function State(; pos=_SZ3, vel=_SZ3)
    State(pos[1], pos[2], pos[3], vel[1], vel[2], vel[3])
end
coordinates(s::State) = @SVector [s.x, s.y, s.z]
velocities(s::State) = @SVector [s.vx, s.vy, s.vz]

function state_unit(us)
    return State(
        length_unit(us),
        length_unit(us),
        length_unit(us),
        velocity_unit(us),
        velocity_unit(us),
        velocity_unit(us),
    )
end