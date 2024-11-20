"""
# Extended help

    rand(::Type{Interval}; [N]::Type{<:Real}=Float64, [dim]::Int=1,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing)

### Notes

All numbers are normally distributed with mean 0 and standard deviation 1.
"""
function rand(::Type{Interval};
              N::Type{<:Real}=Float64,
              dim::Int=1,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int,Nothing}=nothing)
    @assert dim == 1 "cannot create a random Interval of dimension $dim"
    rng = reseed!(rng, seed)
    x = randn(rng, N)
    y = randn(rng, N)
    return x < y ? Interval(x, y) : Interval(y, x)
end
