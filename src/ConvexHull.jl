include("convex_hull_algorithms.jl")

export ConvexHull, CH,
       convex_hull,
       convex_hull!

"""
    ConvexHull{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}

Type that represents the convex hull of the union of two convex sets.

### Fields

- `X` -- convex set
- `Y` -- convex set
"""
struct ConvexHull{N<:Real, S1<:LazySet{N}, S2<:LazySet{N}} <: LazySet{N}
    X::S1
    Y::S2

    # default constructor with dimension check
    ConvexHull{N, S1, S2}(X::S1, Y::S2) where
        {S1<:LazySet{N}, S2<:LazySet{N}} where {N<:Real} =
            dim(X) != dim(Y) ? throw(DimensionMismatch) : new(X, Y)
end
# type-less convenience constructor
ConvexHull(X::S1, Y::S2) where {S1<:LazySet{N}, S2<:LazySet{N}} where {N<:Real} =
    ConvexHull{N, S1, S2}(X, Y)

"""
    CH

Alias for `ConvexHull`.
"""
const CH = ConvexHull

"""
    ConvexHull(X, ∅)

Convex hull of a set with the empty set from the right.

### Input

- `X` -- a convex set
- `∅` -- an empty set

### Output

The given set because the empty set is neutral for the convex hull.
"""
ConvexHull(X::LazySet, ::EmptySet) = X

"""
    ConvexHull(∅, X)

Convex hull of a set with the empty set from the left.

### Input

- `∅` -- an empty set
- `X` -- a convex set

### Output

The given set because the empty set is neutral for the convex hull.
"""
ConvexHull(::EmptySet, X::LazySet) = X

# special case: pure empty set convex hull (we require the same numeric type)
(ConvexHull(∅::E, ::E)) where {E<:EmptySet} = ∅

"""
    dim(ch::ConvexHull)::Int

Return the dimension of a convex hull of two convex sets.

### Input

- `ch` -- convex hull of two convex sets

### Output

The ambient dimension of the convex hull of two convex sets.
"""
function dim(ch::ConvexHull)::Int
    return dim(ch.X)
end

"""
    σ(d::AbstractVector{<:Real}, ch::ConvexHull)::AbstractVector{<:Real}

Return the support vector of a convex hull of two convex sets in a given
direction.

### Input

- `d`  -- direction
- `ch` -- convex hull of two convex sets
"""
function σ(d::AbstractVector{<:Real}, ch::ConvexHull)::AbstractVector{<:Real}
    σ1 = σ(d, ch.X)
    σ2 = σ(d, ch.Y)
    ρ1 = dot(d, σ1)
    ρ2 = dot(d, σ2)
    return ρ1 >= ρ2 ? σ1 : σ2
end
