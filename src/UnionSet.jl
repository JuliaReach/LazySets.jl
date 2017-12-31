import Base.∈

export UnionSet,
       ∈

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
    ∈(x::AbstractVector{<:Real}, U::UnionSet)::Bool

Return whether a given vector is contained in a union of convex sets.

### Input

- `x` -- a vector
- `U` -- a union

### Output

`true` iff ``x ∈ U``.
"""
function ∈(x::AbstractVector{<:Real}, U::UnionSet)::Bool
    for S in U.sets
        if ∈(x, S)
            return true
        end
    end
    return false
end
