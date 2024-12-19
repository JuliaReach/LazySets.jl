"""
# Extended help

    rand(::Type{LineSegment}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
"""
function rand(::Type{LineSegment};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    @assert dim == 2 "cannot create a random LineSegment of dimension $dim"
    rng = reseed!(rng, seed)
    p = randn(rng, N, dim)
    q = randn(rng, N, dim)
    return LineSegment(p, q)
end
