"""
    ConvexHull <: LazySet

Type that represents the convex hull of the union of two convex sets.

FIELDS:

- ``X`` -- a convex set
- ``Y`` -- another convex set
"""
struct ConvexHull <: LazySet
    X::LazySet
    Y::LazySet
    ConvexHull(X, Y) = dim(X) != dim(Y) ? throw(DimensionMismatch) : new(X, Y)
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
    dim(ch.X)
end

"""
    σ(d, P)

Return the support vector of a convex hull in a given direction.

INPUT:

- ``d``  -- direction
- ``ch`` -- the convex hull of two sets
"""
function σ(d::Union{Vector{Float64}, SparseVector{Float64,Int64}}, ch::ConvexHull)::Vector{Float64}
    σ1 = σ(d, ch.X)
    σ2 = σ(d, ch.Y)
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
