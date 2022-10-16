struct Body
    name::String
    gm::Float64
end

struct BodyHistory
    body::Body
    times::Vector{Float64}
    states::Vector{State}
end

BodyHistory(b::Body) = BodyHistory(b, Vector{Float64}(), Vector{State}())

struct Simulation
    histories::Vector{BodyHistory}
end
