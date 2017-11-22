export Singleton

"""
    Singleton <: LazySet

Type that represents a singleton, that is, a set with a unique element.

### Fields

- `element` -- the only element of the set
"""
struct Singleton{N<:Real} <: LazySet
    element::Vector{N}
end

# dimension of a singleton
function dim(s::Singleton)::Int64
    return length(s.element)
end

# support vector of a singleton
function σ(d::AbstractVector{<:Real}, s::LazySets.Singleton)::Vector{<:Real}
    return s.element
end
