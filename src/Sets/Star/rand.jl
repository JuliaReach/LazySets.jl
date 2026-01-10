"""
# Extended help

    rand(::Type{Star}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

### Algorithm

A predicate `P` can be passed directly. If `P` is `nothing` (default), we
generate a random `HPolyhedron` of dimension `dim`.

All numbers are normally distributed with mean 0 and standard deviation 1.
"""
function rand(::Type{Star};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing,
              P::Union{AbstractPolyhedron,Nothing}=nothing)
    rng = reseed!(rng, seed)
    c = randn(rng, N, dim)
    if isnothing(P)
        P = rand(HPolyhedron; N=N, dim=dim, rng=rng, seed=seed)
        # may have no constraints, which is better represented as a Universe
        if isempty(P.constraints)
            P = rand(Universe; N=N, dim=dim, rng=rng, seed=seed)
        end
    end
    V = randn(rng, N, dim, LazySets.API.dim(P))
    return Star(c, V, P)
end
