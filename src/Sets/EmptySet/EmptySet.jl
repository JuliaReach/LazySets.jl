"""
    EmptySet{N} <: ConvexSet{N}

Type that represents the empty set, i.e., the set with no elements.
"""
struct EmptySet{N} <: ConvexSet{N}
    dim::Int
end

# default constructor of type Float64
EmptySet(n::Int) = EmptySet{Float64}(n)

"""
    ∅

Alias for `EmptySet{Float64}`.
"""
const ∅ = EmptySet{Float64}
