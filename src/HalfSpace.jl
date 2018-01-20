import Base.∈

export HalfSpace,
       an_element

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
LazySets.HalfSpace{Float64}([0.0, -1.0], 0.0)
```
"""
struct HalfSpace{N<:Real} <: LazySet{N}
    a::Vector{N}
    b::N
end


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
    σ(d::AbstractVector{<:Real}, hs::HalfSpace)::AbstractVector{<:Real}

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
function σ(d::AbstractVector{N},
           hs::HalfSpace)::AbstractVector{<:Real} where {N<:Real}
    return σ_helper(d, Hyperplane(hs.a, hs.b), "half-space")
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
