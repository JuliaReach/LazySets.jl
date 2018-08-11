export Hyperrectangle,
       low,
       high

"""
    Hyperrectangle{N<:Real} <: AbstractHyperrectangle{N}

Type that represents a hyperrectangle.

A [hyperrectangle](https://en.wikipedia.org/wiki/Hyperrectangle) is the
Cartesian product of one-dimensional intervals.

### Fields

- `center` -- center of the hyperrectangle as a real vector
- `radius` -- radius of the ball as a real vector, i.e., half of its width along
              each coordinate direction
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

"""
    Hyperrectangle(;kwargs...)

Construct a hyperrectangle from keyword arguments.

### Input

- `kwargs` -- keyword arguments; two combinations are allowed:
  1. `center`, `radius` -- vectors
  2. `high`, `low`      -- vectors (if both `center` and `radius` are also
                           defined, those are chosen instead)

### Output

A hyperrectangle.

### Examples

The following three constructions are equivalent:

```jldoctest
julia> c = ones(2);

julia> r = [0.1, 0.2];

julia> l = [0.9, 0.8];

julia> h = [1.1, 1.2];

julia> H1 = Hyperrectangle(c, r)
LazySets.Hyperrectangle{Float64}([1.0, 1.0], [0.1, 0.2])
julia> H2 = Hyperrectangle(center=c, radius=r)
LazySets.Hyperrectangle{Float64}([1.0, 1.0], [0.1, 0.2])
julia> H3 = Hyperrectangle(low=l, high=h)
LazySets.Hyperrectangle{Float64}([1.0, 1.0], [0.1, 0.2])
```
"""
function Hyperrectangle(;kwargs...)
    dict = Dict{Symbol, Any}(kwargs)
    if haskey(dict, :center) && haskey(dict, :radius)
        return Hyperrectangle(dict[:center], dict[:radius])
    elseif haskey(dict, :high) && haskey(dict, :low)
        # compute center and radius from high and low vectors
        center = (dict[:high] .+ dict[:low]) ./ 2
        radius = abs.(dict[:high] .- center)
        return Hyperrectangle(center, radius)
    end
    throw(ArgumentError("invalid arguments for Hyperrectangle: Use either " *
        "'center' and 'radius' or 'high' and 'low'."))
end


# --- AbstractHyperrectangle interface functions ---


"""
    radius_hyperrectangle(H::Hyperrectangle{N}, i::Int)::N where {N<:Real}

Return the box radius of a hyperrectangle in a given dimension.

### Input

- `H` -- hyperrectangle

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


# --- AbstractPointSymmetric interface functions ---


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


# --- Hyperrectangle functions ---


"""
    high(H::Hyperrectangle{N})::Vector{N} where {N<:Real}

Return the higher coordinates of a hyperrectangle.

### Input

- `H` -- hyperrectangle

### Output

A vector with the higher coordinates of the hyperrectangle, one entry per
dimension.
"""
function high(H::Hyperrectangle{N})::Vector{N} where {N<:Real}
    return H.center .+ H.radius
end

"""
    low(H::Hyperrectangle{N})::Vector{N} where {N<:Real}

Return the lower coordinates of a hyperrectangle.

### Input

- `H` -- hyperrectangle

### Output

A vector with the lower coordinates of the hyperrectangle, one entry per
dimension.
"""
function low(H::Hyperrectangle{N})::Vector{N} where {N<:Real}
    return H.center .- H.radius
end
