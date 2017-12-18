import Base: ∈, ⊆

export AbstractSingleton,
       element,
       an_element

"""
    AbstractSingleton{N<:Real} <: AbstractHyperrectangle{N}

Abstract type for sets with a single value.

### Notes

Every concrete `AbstractSingleton` must define the following functions:
- `element(::AbstractSingleton{N})::Vector{N}` -- return the single element
- `element(::AbstractSingleton{N}, i::Int)::N` -- return the single element's
    entry in the `i`-th dimension

```jldoctest
julia> subtypes(AbstractSingleton)
2-element Array{Union{DataType, UnionAll},1}:
 LazySets.Singleton
 LazySets.ZeroSet
```
"""
abstract type AbstractSingleton{N<:Real} <: AbstractHyperrectangle{N} end


# --- AbstractHyperrectangle interface functions ---


"""
    radius_hyperrectangle(S::AbstractSingleton{N}, i::Int)::N where {N<:Real}

Return the box radius of a set with a single value in a given dimension.

### Input

- `S` -- set with a single value

### Output

Zero.
"""
function radius_hyperrectangle(S::AbstractSingleton{N}, i::Int
                              )::N where {N<:Real}
    return zero(N)
end


"""
    radius_hyperrectangle(S::AbstractSingleton{N})::Vector{N} where {N<:Real}

Return the box radius of a set with a single value in every dimension.

### Input

- `S` -- set with a single value

### Output

The zero vector.
"""
function radius_hyperrectangle(S::AbstractSingleton{N}
                              )::Vector{N} where {N<:Real}
    return zeros(N, dim(S))
end


# --- AbstractPointSymmetric interface functions ---


"""
    center(S::AbstractSingleton{N})::Vector{N} where {N<:Real}

Return the center of a set with a single value.

### Input

- `S` -- set with a single value

### Output

The only element of the set.
"""
function center(S::AbstractSingleton{N})::Vector{N} where {N<:Real}
    return element(S)
end


# --- AbstractPolytope interface functions ---


"""
    vertices_list(S::AbstractSingleton{N})::Vector{Vector{N}} where {N<:Real}

Return the list of vertices of a set with a single value.

### Input

- `S` -- set with a single value

### Output

A list containing only a single vertex.
"""
function vertices_list(S::AbstractSingleton{N}
                      )::Vector{Vector{N}} where {N<:Real}
    return [element(S)]
end


# --- LazySet interface functions ---


"""
    σ(d::AbstractVector{N}, S::AbstractSingleton{N})::Vector{N} where {N<:Real}

Return the support vector of a set with a single value.

### Input

- `d` -- direction
- `S` -- set with a single value

### Output

The support vector, which is the set's vector itself, irrespective of the given
direction.
"""
function σ(d::AbstractVector{N},
           S::AbstractSingleton{N})::Vector{N} where {N<:Real}
    return element(S)
end


"""
    an_element(S::AbstractSingleton{N})::Vector{N} where {N<:Real}

Return some element of a set with a single value.

### Input

- `S` -- set with a single value

### Output

The only element in the set.
"""
function an_element(S::AbstractSingleton{N})::Vector{N} where {N<:Real}
    return element(S)
end


"""
    ∈(x::AbstractVector{N}, S::AbstractSingleton{N})::Bool where {N<:Real}

Check whether a given point is contained in a set with a single value.

### Input

- `x` -- point/vector
- `S` -- set with a single value

### Output

`true` iff ``x ∈ S``.

### Notes

This implementation performs an exact comparison, which may be insufficient with
floating point computations.
"""
function ∈(x::AbstractVector{N}, S::AbstractSingleton{N})::Bool where {N<:Real}
    return x == element(S)
end


"""
    ⊆(S::AbstractSingleton, set::LazySet)::Bool

Check whether a given set with a single value is contained in a convex set.

### Input

- `S`   -- set with a single value
- `set` -- convex set

### Output

`true` iff ``S ⊆ \\text{set}``.
"""
function ⊆(S::AbstractSingleton, set::LazySet)::Bool
    return ∈(element(S), set)
end


"""
    ⊆(S1::AbstractSingleton, S2::AbstractSingleton)::Bool

Check whether a given set with a single value is contained in another set with a
single value.

### Input

- `S1` -- first set with a single value (containee?)
- `S2` -- second set with a single value (containee?)

### Output

`true` iff ``S1 ⊆ S2`` iff ``S1 == S2``.
"""
function ⊆(S1::AbstractSingleton, S2::AbstractSingleton)::Bool
    return element(S1) == element(S2)
end
