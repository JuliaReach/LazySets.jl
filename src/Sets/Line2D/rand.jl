"""
# Extended help

    rand(::Type{Line2D}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
Additionally, the constraint `a` is nonzero.
"""
function rand(::Type{Line2D};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    @assert dim == 2 "cannot create a random Line2D of dimension $dim"
    rng = reseed!(rng, seed)
    a = randn(rng, N, dim)
    while iszero(a)
        a = randn(rng, N, dim)  # COV_EXCL_LINE
    end  # COV_EXCL_LINE
    b = randn(rng, N)
    return Line2D(a, b)
end
