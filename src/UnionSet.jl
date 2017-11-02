export UnionSet, is_contained

"""
    UnionSet <: LazySet

Type that represents a union of convex sets.

### Fields

- `sets` --  a list of convex sets
"""
struct UnionSet{T<:LazySet} <: LazySet
    sets::Vector{T}
end
UnionSet(sets::Vector{T}) where {T<:LazySet} = UnionSet{T}(sets)

"""
    is_contained(x, U)

States if a given vector belongs to a given union of sets.

### Input

- `x` -- a vector
- `U` -- a union

### Output

Return true iff ``x ∈ U``.
"""
function is_contained(x::Vector{Float64}, U::UnionSet)::Bool
    res = false
    for s in U.sets
        res = res || in(x, s)
    end
    return res
end

