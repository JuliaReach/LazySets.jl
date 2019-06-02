import Base: rand,
             ∈,
             isempty

export Universe

"""
    Universe{N<:Real} <: AbstractPolyhedron{N}

Type that represents the universal set, i.e., the set of all elements.
"""
struct Universe{N<:Real} <: AbstractPolyhedron{N}
    dim::Int
end

# default constructor of type Float64
Universe(dim::Int) = Universe{Float64}(dim)


# --- AbstractPolyhedron interface functions ---


"""
    constraints_list(U::Universe{N}) where {N<:Real}

Return the list of constraints defining a universe.

### Input

- `U` -- universe

### Output

The empty list of constraints, as the universe is unconstrained.
"""
function constraints_list(U::Universe{N}) where {N<:Real}
    return LinearConstraint{N, Vector{N}}[]
end

"""
    constrained_dimensions(U::Universe)::Vector{Int}

Return the indices in which a universe is constrained.

### Input

- `U` -- universe

### Output

The empty vector, as the universe is unconstrained in every dimension.
"""
function constrained_dimensions(U::Universe)::Vector{Int}
    return Int[]
end


# --- LazySet interface functions ---


"""
    dim(U::Universe)

Return the dimension of a universe.

### Input

- `U` -- universe

### Output

The dimension of a universe.
"""
function dim(U::Universe)::Int
    return U.dim
end

"""
    ρ(d::AbstractVector{N}, U::Universe{N}) where {N<:Real}

Return the support function of a universe.

### Input

- `d` -- direction
- `U` -- universe

### Output

The support function in the given direction.

### Algorithm

If the direction is all zero, the result is zero.
Otherwise, the result is `Inf`.
"""
function ρ(d::AbstractVector{N}, U::Universe{N}) where {N<:Real}
    return iszero(d) ? zero(N) : N(Inf)
end

"""
    σ(d::AbstractVector{N}, U::Universe{N}) where {N<:Real}

Return the support vector of a universe.

### Input

- `d` -- direction
- `U` -- universe

### Output

A vector with infinity values, except in dimensions where the direction is zero.
"""
function σ(d::AbstractVector{N}, U::Universe{N}) where {N<:Real}
    return [v == zero(N) ? v : v > zero(N) ? N(Inf) : N(-Inf) for v in d]
end

"""
    ∈(x::AbstractVector{N}, U::Universe{N})::Bool where {N<:Real}

Check whether a given point is contained in a universe.

### Input

- `x` -- point/vector
- `U` -- universe

### Output

The output is always `true`.

### Examples

```jldoctest
julia> ∈([1.0, 0.0], Universe(2))
true
```
"""
function ∈(x::AbstractVector{N}, U::Universe{N})::Bool where {N<:Real}
    @assert length(x) == dim(U)
    return true
end

"""
    an_element(U::Universe{N}) where {N<:Real}

Return some element of a universe.

### Input

- `U` -- universe

### Output

The origin.
"""
function an_element(U::Universe{N}) where {N<:Real}
    return zeros(N, dim(U))
end

"""
    rand(::Type{Universe}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::Universe{N}

Create a universe (note that there is nothing to randomize).

### Input

- `Universe` -- type for dispatch
- `N`        -- (optional, default: `Float64`) numeric type
- `dim`      -- (optional, default: 2) dimension
- `rng`      -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`     -- (optional, default: `nothing`) seed for reseeding

### Output

The (only) universe of the given numeric type and dimension.
"""
function rand(::Type{Universe};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::Universe{N}
    rng = reseed(rng, seed)
    return Universe{N}(dim)
end

"""
    isempty(U::Universe)::Bool

Return if a universe is empty or not.

### Input

- `U` -- universe

### Output

`false`.
"""
function isempty(U::Universe)::Bool
    return false
end

"""
    isbounded(U::Universe)::Bool

Determine whether a universe is bounded.

### Input

- `S` -- universe

### Output

`false` as the universe is unbounded.
"""
function isbounded(U::Universe)::Bool
    return false
end

"""
    norm(U::Universe, [p]::Real=Inf)

Return the norm of a universe.
It is the norm of the enclosing ball (of the given ``p``-norm) of minimal volume
that is centered in the origin.

### Input

- `U` -- universe
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function norm(U::Universe, p::Real=Inf)
    error("a universe does not have a norm")
end

"""
    radius(U::Universe, [p]::Real=Inf)

Return the radius of a universe.
It is the radius of the enclosing ball (of the given ``p``-norm) of minimal
volume with the same center.

### Input

- `U` -- universe
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function radius(U::Universe, p::Real=Inf)
    error("a universe does not have a radius")
end

"""
    diameter(U::Universe, [p]::Real=Inf)

Return the diameter of a universe.
It is the maximum distance between any two elements of the set, or,
equivalently, the diameter of the enclosing ball (of the given ``p``-norm) of
minimal volume with the same center.

### Input

- `U` -- universe
- `p` -- (optional, default: `Inf`) norm

### Output

An error.
"""
function diameter(U::Universe, p::Real=Inf)
    error("a universe does not have a diameter")
end

"""
    translate(U::Universe{N}, v::AbstractVector{N}) where {N<:Real}

Translate (i.e., shift) a universe by a given vector.

### Input

- `U` -- universe
- `v` -- translation vector

### Output

The universe.
"""
function translate(U::Universe{N}, v::AbstractVector{N}) where {N<:Real}
    @assert length(v) == dim(U) "cannot translate a $(dim(U))-dimensional " *
                                "set by a $(length(v))-dimensional vector"
    return U
end
