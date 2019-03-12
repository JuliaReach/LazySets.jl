import Base.rand

export Hyperrectangle

"""
    Hyperrectangle{N<:Real} <: AbstractHyperrectangle{N}

Type that represents a hyperrectangle.

A [hyperrectangle](https://en.wikipedia.org/wiki/Hyperrectangle) is the
Cartesian product of one-dimensional intervals.

### Fields

- `center` -- center of the hyperrectangle as a real vector
- `radius` -- radius of the ball as a real vector, i.e., half of its width along
              each coordinate direction

### Examples

There is also a constructor from lower and upper bounds with keyword arguments
`high` and `low`.
The following two constructions are equivalent:

```jldoctest
julia> c = ones(2);

julia> r = [0.1, 0.2];

julia> l = [0.9, 0.8];

julia> h = [1.1, 1.2];

julia> Hyperrectangle(c, r)
Hyperrectangle{Float64}([1.0, 1.0], [0.1, 0.2])
julia> Hyperrectangle(low=l, high=h)
Hyperrectangle{Float64}([1.0, 1.0], [0.1, 0.2])
```
"""
struct Hyperrectangle{N<:Real} <: AbstractHyperrectangle{N}
    center::Vector{N}
    radius::Vector{N}

    # default constructor with length comparison & domain constraint for radius
    function Hyperrectangle{N}(center::Vector{N},
                               radius::Vector{N}) where {N<:Real}
        @assert length(center) == length(radius) "length of center and " *
            "radius must be equal"
        @assert all(v -> v >= zero(N), radius) "radius must not be negative"
        return new{N}(center, radius)
    end
end

# convenience constructor without type parameter
Hyperrectangle(center::Vector{N}, radius::Vector{N}) where {N<:Real} =
    Hyperrectangle{N}(center, radius)

# constructor from keyword arguments (lower and upper bounds)
@static if VERSION < v"0.7-"
@eval begin

function Hyperrectangle(;kwargs...)
    dict = Dict{Symbol, Any}(kwargs)
    if haskey(dict, :high) && haskey(dict, :low)
        # compute center and radius from high and low vectors
        high = dict[:high]
        center = (high .+ dict[:low]) ./ 2
        radius = abs.(high .- center)
        return Hyperrectangle(center, radius)
    end
    throw(ArgumentError("invalid arguments for Hyperrectangle: " *
        "use 'high' and 'low'"))
end

end # @eval
else
@eval begin

function Hyperrectangle(;
                        high::AbstractVector{N},
                        low::AbstractVector{N}) where {N<:Real}
    # compute center and radius from high and low vectors
    center = (high .+ low) ./ 2
    radius = abs.(high .- center)
    return Hyperrectangle(center, radius)
end

end # @eval
end # if


# --- AbstractHyperrectangle interface functions ---


"""
    radius_hyperrectangle(H::Hyperrectangle{N}, i::Int)::N where {N<:Real}

Return the box radius of a hyperrectangle in a given dimension.

### Input

- `H` -- hyperrectangle
- `i` -- dimension of interest

### Output

The radius in the given dimension.
"""
function radius_hyperrectangle(H::Hyperrectangle{N}, i::Int)::N where {N<:Real}
    return H.radius[i]
end

"""
    radius_hyperrectangle(H::Hyperrectangle{N})::Vector{N} where {N<:Real}

Return the box radius of a hyperrectangle in every dimension.

### Input

- `H` -- hyperrectangle

### Output

The box radius of the hyperrectangle.
"""
function radius_hyperrectangle(H::Hyperrectangle{N})::Vector{N} where {N<:Real}
    return H.radius
end


# --- AbstractCentrallySymmetric interface functions ---


"""
    center(H::Hyperrectangle{N})::Vector{N} where {N<:Real}

Return the center of a hyperrectangle.

### Input

- `H` -- hyperrectangle

### Output

The center of the hyperrectangle.
"""
function center(H::Hyperrectangle{N})::Vector{N} where {N<:Real}
    return H.center
end


# --- LazySet interface functions ---


"""
    rand(::Type{Hyperrectangle}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::Hyperrectangle{N}

Create a random hyperrectangle.

### Input

- `Hyperrectangle` -- type for dispatch
- `N`              -- (optional, default: `Float64`) numeric type
- `dim`            -- (optional, default: 2) dimension
- `rng`            -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`           -- (optional, default: `nothing`) seed for reseeding

### Output

A random hyperrectangle.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
Additionally, the radius is nonnegative.
"""
function rand(::Type{Hyperrectangle};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::Hyperrectangle{N}
    rng = reseed(rng, seed)
    center = randn(rng, N, dim)
    radius = abs.(randn(rng, N, dim))
    return Hyperrectangle(center, radius)
end

"""
    translate(H::Hyperrectangle{N}, v::AbstractVector{N}; share::Bool=false
             ) where {N<:Real}

Translate (i.e., shift) a hyperrectangle by a given vector.

### Input

- `H`     -- hyperrectangle
- `v`     -- translation vector
- `share` -- (optional, default: `false`) flag for sharing unmodified parts of
             the original set representation

### Output

A translated hyperrectangle.

### Notes

The radius vector is shared with the original hyperrectangle if `share == true`.

### Algorithm

We add the vector to the center of the hyperrectangle.
"""
function translate(H::Hyperrectangle{N}, v::AbstractVector{N}; share::Bool=false
                  ) where {N<:Real}
    @assert length(v) == dim(H) "cannot translate a $(dim(H))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    c = center(H) + v
    radius = share ? H.radius : copy(H.radius)
    return Hyperrectangle(c, radius)
end
