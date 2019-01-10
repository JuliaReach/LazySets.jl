import Base: rand,
             ∈,
             isempty,
             convert

export HalfSpace, LinearConstraint,
       an_element,
       constrained_dimensions,
       halfspace_left, halfspace_right

"""
    HalfSpace{N<:Real} <: LazySet{N}

Type that represents a (closed) half-space of the form ``a⋅x ≤ b``.

### Fields

- `a` -- normal direction
- `b` -- constraint

### Examples

The set ``y ≥ 0`` in the plane:

```jldoctest
julia> HalfSpace([0, -1.], 0.)
HalfSpace{Float64}([0.0, -1.0], 0.0)
```
"""
struct HalfSpace{N<:Real} <: LazySet{N}
    a::AbstractVector{N}
    b::N
end

function convert(::Type{HalfSpace{N}}, hs::HalfSpace) where {N<:Real}
    return HalfSpace{N}(hs.a, hs.b)
end

"""
    LinearConstraint

Alias for `HalfSpace`
"""
const LinearConstraint = HalfSpace


# --- LazySet interface functions ---


"""
    dim(hs::HalfSpace)::Int

Return the dimension of a half-space.

### Input

- `hs` -- half-space

### Output

The ambient dimension of the half-space.
"""
function dim(hs::HalfSpace)::Int
    return length(hs.a)
end

"""
    ρ(d::AbstractVector{N}, hs::HalfSpace{N})::N where {N<:Real}

Evaluate the support function of a half-space in a given direction.

### Input

- `d`  -- direction
- `hs` -- half-space

### Output

The support function of the half-space.
If the set is unbounded in the given direction, the result is `Inf`.
"""
function ρ(d::AbstractVector{N}, hs::HalfSpace{N})::N where {N<:Real}
    v, unbounded = σ_helper(d, Hyperplane(hs.a, hs.b); error_unbounded=false,
                            halfspace=true)
    if unbounded
        return N(Inf)
    end
    return dot(d, v)
end

"""
    σ(d::AbstractVector{N}, hs::HalfSpace{N}) where {N<:Real}

Return the support vector of a half-space.

### Input

- `d`  -- direction
- `hs` -- half-space

### Output

The support vector in the given direction, which is only defined in the
following two cases:
1. The direction has norm zero.
2. The direction is the half-space's normal direction.
In both cases the result is any point on the boundary (the defining hyperplane).
Otherwise this function throws an error.
"""
function σ(d::AbstractVector{N}, hs::HalfSpace{N}) where {N<:Real}
    v, unbounded = σ_helper(d, Hyperplane(hs.a, hs.b); error_unbounded=true,
                            halfspace=true)
    return v
end

"""
    isbounded(hs::HalfSpace)::Bool

Determine whether a half-space is bounded.

### Input

- `hs` -- half-space

### Output

`false`.
"""
function isbounded(::HalfSpace)::Bool
    return false
end

"""
    an_element(hs::HalfSpace{N})::Vector{N} where {N<:Real}

Return some element of a half-space.

### Input

- `hs` -- half-space

### Output

An element on the defining hyperplane.
"""
function an_element(hs::HalfSpace{N})::Vector{N} where {N<:Real}
    return an_element_helper(Hyperplane(hs.a, hs.b))
end

"""
    ∈(x::AbstractVector{N}, hs::HalfSpace{N})::Bool where {N<:Real}

Check whether a given point is contained in a half-space.

### Input

- `x` -- point/vector
- `hs` -- half-space

### Output

`true` iff ``x ∈ hs``.

### Algorithm

We just check if ``x`` satisfies ``a⋅x ≤ b``.
"""
function ∈(x::AbstractVector{N}, hs::HalfSpace{N})::Bool where {N<:Real}
    return dot(x, hs.a) <= hs.b
end

"""
    rand(::Type{HalfSpace}; [N]::Type{<:Real}=Float64, [dim]::Int=2,
         [rng]::AbstractRNG=GLOBAL_RNG, [seed]::Union{Int, Nothing}=nothing
        )::HalfSpace{N}

Create a random half-space.

### Input

- `HalfSpace` -- type for dispatch
- `N`         -- (optional, default: `Float64`) numeric type
- `dim`       -- (optional, default: 2) dimension
- `rng`       -- (optional, default: `GLOBAL_RNG`) random number generator
- `seed`      -- (optional, default: `nothing`) seed for reseeding

### Output

A random half-space.

### Algorithm

All numbers are normally distributed with mean 0 and standard deviation 1.
Additionally, the constraint `a` is nonzero.
"""
function rand(::Type{HalfSpace};
              N::Type{<:Real}=Float64,
              dim::Int=2,
              rng::AbstractRNG=GLOBAL_RNG,
              seed::Union{Int, Nothing}=nothing
             )::HalfSpace{N}
    rng = reseed(rng, seed)
    a = randn(rng, N, dim)
    while iszero(a)
        a = randn(rng, N, dim)
    end
    b = randn(rng, N)
    return HalfSpace(a, b)
end

"""
    isempty(hs::HalfSpace)::Bool

Return if a half-space is empty or not.

### Input

- `hs` -- half-space

### Output

`false`.
"""
function isempty(hs::HalfSpace)::Bool
    return false
end

"""
    constraints_list(hs::HalfSpace{N})::Vector{LinearConstraint{N}}
        where {N<:Real}

Return the list of constraints of a half-space.

### Input

- `hs` -- half-space

### Output

A singleton list containing the half-space.
"""
function constraints_list(hs::HalfSpace{N}
                         )::Vector{LinearConstraint{N}} where {N<:Real}
    return [hs]
end

"""
    constrained_dimensions(hs::HalfSpace{N})::Vector{Int} where {N<:Real}

Return the indices in which a half-space is constrained.

### Input

- `hs` -- half-space

### Output

A vector of ascending indices `i` such that the half-space is constrained in
dimension `i`.

### Examples

A 2D half-space with constraint ``x1 ≥ 0`` is constrained in dimension 1 only.
"""
function constrained_dimensions(hs::HalfSpace{N})::Vector{Int} where {N<:Real}
    return nonzero_indices(hs.a)
end

"""
    halfspace_left(p::AbstractVector{N},
                   q::AbstractVector{N})::HalfSpace{N} where {N<:Real}

Return a half-space describing the 'left' of a two-dimensional line segment
through two points.

### Input

- `p` -- first point
- `q` -- second point

### Output

The half-space whose boundary goes through the two points `p` and `q` and which
describes the left-hand side of the directed line segment `pq`.

### Algorithm

The implementation is simple: the half-space ``a⋅x ≤ b`` is calculated as
`a = [dy, -dx]`, where ``d = (dx, dy)`` denotes the line segment
`pq`, that is, ``\\vec{d} = \\vec{p} - \\vec{q}``, and `b = dot(p, a)`.

### Examples

The left half-space of the "east" and "west" directions in two-dimensions are
the upper and lower half-spaces:

```jldoctest halfspace_left
julia> import LazySets.halfspace_left

julia> halfspace_left([0.0, 0.0], [1.0, 0.0])
HalfSpace{Float64}([0.0, -1.0], 0.0)

julia> halfspace_left([0.0, 0.0], [-1.0, 0.0])
HalfSpace{Float64}([0.0, 1.0], 0.0)
```

We create a box from the sequence of line segments that describe its edges:

```jldoctest halfspace_left
julia> H1 = halfspace_left([-1.0, -1.0], [1.0, -1.0]);

julia> H2 = halfspace_left([1.0, -1.0], [1.0, 1.0]);

julia> H3 = halfspace_left([1.0, 1.0], [-1.0, 1.0]);

julia> H4 = halfspace_left([-1.0, 1.0], [-1.0, -1.0]);

julia> H = HPolygon([H1, H2, H3, H4]);

julia> B = BallInf([0.0, 0.0], 1.0);

julia> B ⊆ H && H ⊆ B
true
```
"""
function halfspace_left(p::AbstractVector{N},
                        q::AbstractVector{N})::HalfSpace{N} where {N<:Real}
    @assert length(p) == length(q) == 2 "the points must be two-dimensional"
    @assert p != q "the points must not be equal"
    a = [q[2] - p[2], p[1] - q[1]]
    return HalfSpace(a, dot(p, a))
end

"""
    halfspace_right(p::AbstractVector{N},
                    q::AbstractVector{N})::HalfSpace{N} where {N<:Real}

Return a half-space describing the 'right' of a two-dimensional line segment
through two points.

### Input

- `p` -- first point
- `q` -- second point

### Output

The half-space whose boundary goes through the two points `p` and `q` and which
describes the right-hand side of the directed line segment `pq`.

### Algorithm

See the documentation of `halfspace_left`. 
"""
function halfspace_right(p::AbstractVector{N},
                         q::AbstractVector{N})::HalfSpace{N} where {N<:Real}
    return halfspace_left(q, p)
end
