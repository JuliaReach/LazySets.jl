"""
# Extended help

    rand(::Type{ZeroSet}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

### Output

The (only) zero set of the given numeric type and dimension.
"""
function rand(::Type{ZeroSet};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    return ZeroSet{N}(dim)
end
