export Singleton

"""
    Singleton <: LazySet

Type that represents a singleton, that is, a set with a unique element.

### Fields

- `element` -- the only element of the set
"""
struct Singleton <: LazySet
    element::Vector{Float64}
end

# dimension of a singleton
function dim(s::Singleton)::Int64
    return length(s.element)
end

# support vector of a singleton
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, s::LazySets.Singleton)::Vector{Float64}
    return s.element
end
