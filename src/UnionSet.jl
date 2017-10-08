"""
    UnionSet

Type that represents a union of convex sets.

FIELDS:

- ``sets`` --  a list of convex sets
"""
struct UnionSet
    sets::Array{LazySet, 1}
end

"""
    is_contained(x, U)

States if a given vector belongs to a given union of sets.

INPUT:

- ``x`` -- a vector
- ``U`` -- a union

OUTPUT:

Return true iff x ∈ U.
"""
function is_contained(x::Vector{Float64}, U::UnionSet)::Bool
    res = false
    for s in U.sets
        res = res || in(x, s)
    end
    return res
end

export UnionSet, is_contained
