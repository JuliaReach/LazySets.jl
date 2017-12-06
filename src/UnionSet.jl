export UnionSet,
       is_contained

"""
    UnionSet{S<:LazySet} <: LazySet

Type that represents a union of convex sets.

### Fields

- `sets` --  a list of convex sets
"""
struct UnionSet{S<:LazySet} <: LazySet
    sets::Vector{S}
end

"""
    is_contained(x::AbstractVector{<:Real}, U::UnionSet)::Bool

Return whether a given vector is contained in a union of convex sets.

### Input

- `x` -- a vector
- `U` -- a union

### Output

Return `true` iff ``x ∈ U``.
"""
function is_contained(x::AbstractVector{<:Real}, U::UnionSet)::Bool
    for s in U.sets
        if in(x, s)
            return true
        end
    end
    return false
end
