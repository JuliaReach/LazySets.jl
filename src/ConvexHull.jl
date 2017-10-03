"""
    ConvexHull <: LazySet

Type that represents the convex hull of the union of two convex sets.

FIELDS:

- ``s1`` -- a convex set
- ``s2`` -- another convex set
"""
struct ConvexHull <: LazySet
    s1::LazySet
    s2::LazySet
    ConvexHull(s1, s2) = dim(s1) != dim(s2) ? throw(DimensionMismatch) : new(s1, s2)
end
# function alias
CH = ConvexHull

"""
    dim(P)

Return the ambient dimension of the convex hull of two sets.

INPUT:

- ``ch`` -- the convex hull of two sets
"""
function dim(ch::ConvexHull)::Int64
    dim(ch.s1)
end

"""
    σ(d, P)

Return the support vector of a convex hull in a given direction.

INPUT:

- ``d``  -- direction
- ``ch`` -- the convex hull of two sets
"""
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, ch::ConvexHull)::Vector{Float64}
    σ1 = σ(d, ch.s1)
    σ2 = σ(d, ch.s2)
    ρ1 = dot(d, σ1)::Float64
    ρ2 = dot(d, σ2)::Float64
    if ρ1 >= ρ2
        res = σ1
    else
        res = σ2
    end
    return res
end

export ConvexHull, CH
