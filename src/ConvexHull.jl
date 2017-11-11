include("convex_hull_algorithms.jl")

export ConvexHull,    # convex hull lazy type
       CH,            # type alias
       convex_hull,   # explicit convex hull of a list of points in the plane
       convex_hull!

"""
    ConvexHull <: LazySet

Type that represents the convex hull of the union of two convex sets.

### Fields

- `X` -- a convex set
- `Y` -- another convex set
"""
struct ConvexHull{T1<:LazySet,T2<:LazySet} <: LazySet
    X::T1
    Y::T2
    ConvexHull{T1,T2}(X::T1, Y::T2) where {T1<:LazySet,T2<:LazySet} = dim(X) != dim(Y) ? throw(DimensionMismatch) : new(X, Y)
end
ConvexHull(X::T1, Y::T2) where {T1<:LazySet,T2<:LazySet} = ConvexHull{T1,T2}(X, Y)
# function alias
CH = ConvexHull

"""
    dim(P)

Return the ambient dimension of the convex hull of two sets.

### Input

- `ch` -- the convex hull of two sets
"""
function dim(ch::ConvexHull)::Int64
    dim(ch.X)
end

"""
    σ(d, P)

Return the support vector of a convex hull in a given direction.

### Input

- `d`  -- direction
- `ch` -- the convex hull of two sets
"""
function σ(d::AbstractVector{Float64}, ch::ConvexHull)::Vector{Float64}
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
