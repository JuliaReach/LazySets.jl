import Base.∈

export Singleton

"""
    Singleton{N<:Real} <: LazySet

Type that represents a singleton, that is, a set with a unique element.

### Fields

- `element` -- the only element of the set
"""
struct Singleton{N<:Real} <: LazySet
    element::Vector{N}
end

"""
    dim(S::Singleton)::Int

Return the dimension of a singleton.

### Input

- `S` -- singleton

### Output

The ambient dimension of the singleton.
"""
function dim(S::Singleton)::Int
    return length(S.element)
end

"""
    σ(d::AbstractVector{<:Real}, S::LazySets.Singleton{N})::Vector{N} where {N<:Real}

Return the support vector of a singleton.

### Input

- `d` -- direction
- `B` -- singleton

### Output

The support vector, which is the singleton's vector itself, irrespective of the
given direction.
"""
function σ(d::AbstractVector{<:Real},
           S::LazySets.Singleton{N})::Vector{N} where {N<:Real}
    return S.element
end

"""
    ∈(x::AbstractVector{N}, S::Singleton{N})::Bool where {N<:Real}

Check whether a given point is contained in a singleton.

### Input

- `x` -- point/vector
- `S` -- singleton

### Output

`true` iff ``x ∈ S``.

### Notes

This implementation performs an exact comparison, which may be insufficient with
floating point computations.

### Examples

```jldoctest
julia> S = Singleton([1., 1.]);

julia> ∈([0.9, 1.1], S)
false
julia> ∈([1.0, 1.0], S)
true
```
"""
function ∈(x::AbstractVector{N}, S::Singleton{N})::Bool where {N<:Real}
    return x == S.element
end

"""
    ∈(x::Singleton, set::LazySet)::Bool

Check whether a given singleton is contained in a convex set.

### Input

- `x`   -- singleton
- `set` -- convex set

### Output

`true` iff ``x ∈ \text{set}``.
"""
function ∈(x::Singleton, set::LazySet)::Bool
    return ∈(x.element, set)
end
