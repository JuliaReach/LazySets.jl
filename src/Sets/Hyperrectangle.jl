import Base.rand

export Hyperrectangle

"""
    Hyperrectangle{N<:Real, VNC<:AbstractVector{N}, VNR<:AbstractVector{N}
                  } <: AbstractHyperrectangle{N}

Type that represents a hyperrectangle.

A [hyperrectangle](https://en.wikipedia.org/wiki/Hyperrectangle) is the
Cartesian product of one-dimensional intervals.

### Fields

- `center` -- center of the hyperrectangle as a real vector
- `radius` -- radius of the ball as a real vector, i.e., half of its width along
              each coordinate direction

### Examples

The `Hyperrectangle` type stores a vector representing the center and another
vector representing the radius. The default constructor `Hyperrectangle(c, r)`
receives the center and radius, in that order. For instance,

```jldoctest hyperrectangle_constructor
julia> c = [-1.0, 1.0];

julia> r = [2.0, 1.0];

julia> H = Hyperrectangle(c, r)
Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}([-1.0, 1.0], [2.0, 1.0])
```

Which creates the hyperrectangle with vertices:

```jldoctest hyperrectangle_constructor
julia> vertices_list(H)
4-element Array{Array{Float64,1},1}:
 [1.0, 2.0]
 [-3.0, 2.0]
 [1.0, 0.0]
 [-3.0, 0.0]
```

The getter functions for the center and the radius are `center` and
`radius_hyperrectangle` (since `radius` corresponds to the radius of the
enclosing ball of minimal volume):

```jldoctest hyperrectangle_constructor
julia> center(H)
2-element Array{Float64,1}:
 -1.0
  1.0

julia> radius_hyperrectangle(H)
2-element Array{Float64,1}:
 2.0
 1.0
```

There is also a constructor from lower and upper bounds with keyword arguments
`high` and `low`. The following construction results in the same hyperrectangle
as in the previous paragraph:

```jldoctest hyperrectangle_constructor
julia> l = [-3.0, 0.0];

julia> h = [1.0, 2.0];

julia> Hyperrectangle(low=l, high=h)
Hyperrectangle{Float64,Array{Float64,1},Array{Float64,1}}([-1.0, 1.0], [2.0, 1.0])
```
"""
struct Hyperrectangle{N<:Real, VNC<:AbstractVector{N}, VNR<:AbstractVector{N}
                     } <: AbstractHyperrectangle{N}
    center::VNC
    radius::VNR

    # default constructor with length comparison & domain constraint for radius
    function Hyperrectangle(center::VNC, radius::VNR) where
            {N<:Real, VNC<:AbstractVector{N}, VNR<:AbstractVector{N}}
        @assert length(center) == length(radius) "length of center and " *
            "radius must be equal"
        @assert all(v -> v >= zero(N), radius) "radius must not be negative"
        return new{N, VNC, VNR}(center, radius)
    end
end

isoperationtype(::Type{<:Hyperrectangle}) = false
isconvextype(::Type{<:Hyperrectangle}) = true

# constructor from keyword arguments (lower and upper bounds)
function Hyperrectangle(;
                        high::AbstractVector{N},
                        low::AbstractVector{N}) where {N<:Real}
    @assert all(i -> low[i] <= high[i], eachindex(low)) "lower bound must be lower than upper bound"
    # compute center and radius from high and low vectors
    center = (high .+ low) ./ 2
    radius = high .- center
    return Hyperrectangle(center, radius)
end


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

# Particular dispatch for SingleEntryVector
function ρ(d::SingleEntryVector{N}, H::Hyperrectangle{N}) where {N<:Real}
    return _ρ_sev_hyperrectangle(d, H)
end
