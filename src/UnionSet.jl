import Base.∈

export UnionSet

"""
    UnionSet{N<:Real, S<:LazySet{N}}

Type that represents a union of sets.

### Fields

- `sets` --  a list of sets
"""
struct UnionSet{N<:Real, S<:LazySet{N}}
    sets::Vector{S}
end

"""
    ∈(x::AbstractVector{<:Real}, U::UnionSet)::Bool

Return whether a given vector is contained in a union of sets.

### Input

- `x` -- vector
- `U` -- union

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
