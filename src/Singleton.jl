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
    dim(S::Singleton)

Return the dimension of a singleton.

### Input

- `S` -- singleton

### Output

The ambient dimension of the singleton.
"""
function dim(S::Singleton)
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
