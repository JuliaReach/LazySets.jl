"""
# Extended help

    rand(::Type{HPolyhedron}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

### Algorithm

We first create a random polytope and then for each constraint randomly (50%)
decide whether to include it.
"""
function rand(::Type{HPolyhedron};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    P = rand(HPolytope; N=N, dim=dim, rng=rng)
    constraints_P = constraints_list(P)
    constraints_Q = Vector{eltype(constraints_P)}()
    for i in eachindex(constraints_P)
        if rand(rng, Bool)
            push!(constraints_Q, constraints_P[i])
        end
    end
    return HPolyhedron(constraints_Q)
end
