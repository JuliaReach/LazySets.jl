"""
# Extended help

    rand(::Type{Line}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
"""
function rand(::Type{Line};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    d = randn(rng, N, dim)
    while iszero(d)
        d = randn(rng, N, dim)  # COV_EXCL_LINE
    end  # COV_EXCL_LINE
    p = randn(rng, N, dim)
    return Line(p, d)
end
