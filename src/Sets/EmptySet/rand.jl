"""
# Extended help

    rand(::Type{EmptySet}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

### Output

The (unique) empty set of the given numeric type and dimension.
"""
function rand(::Type{EmptySet};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    rng = reseed!(rng, seed)
    return EmptySet{N}(dim)
end
